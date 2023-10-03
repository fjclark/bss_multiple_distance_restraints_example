"""
Set up and run the bound leg of an ABFE calculation 
with multiple distance restraints using BioSimSpace.
"""

import argparse
import BioSimSpace.Sandpit.Exscientia as BSS
import logging
import os
from subprocess import Popen, PIPE, STDOUT, run
from typing import Optional, List

# Set up logging
logger = logging.getLogger(__name__)

# Dict of allowed stage names and the corresponding perturbation type.
PERTURBATION_TYPES = {
    "restrain": "restraint",
    "discharge": "discharge_soft",
    "vanish": "vanish_soft",
    "release_restraint": "release_restraint",
}


def setup_abfe_bound_leg(
    topology_file: str,
    coordinates_file: str,
    stage_options: dict,
    trajectory_file: Optional[str] = None,
    restraint_type: str = "multiple_distance",
) -> None:
    """
    Set up the bound leg of the ABFE calculation with multiple distance restraints.

    Parameters
    ----------
    topology_file : str
        The path to the topology file.
    coordinates_file : str
        The path to the coordinates file.
    stage_options : dict
        A dictionary of options for the stage, taking the form {"stage_name": {"option": value}, ...}.
        The stage names must be one of "restrain", "discharge", "vanish", or "release_restraint".
    trajectory_file : str, optional
        The path to the trajectory file to derive the restraints.
        If not supplied, the restraints will be derived from a short
        trajectory.
    restraint_type : str, optional
        The type of restraint to use. Must be one of "distance" or "multiple_distance".
    """
    # Load the system
    logger.info("Loading system...")
    system = BSS.IO.readMolecules([topology_file, coordinates_file])
    # Assume that the ligand is called "LIG"
    for i, mol in enumerate(system):
        if mol._getSireObject().name().value() == "LIG":
            lig = BSS.Align.decouple(mol)
            logger.info(f"Identified ligand as molecule {i}: {lig}")
            system.updateMolecule(i, lig)

    # If we don't have a trajectory file, run a short simulation to derive the restraints.
    if trajectory_file is None:
        logger.info(
            "No trajectory file provided. Running short simulation to select the restraints..."
        )
        protocol = BSS.Protocol.Production(
            timestep=2 * BSS.Units.Time.femtosecond,
            runtime=5 * BSS.Units.Time.nanosecond,
        )
        restraint_search = BSS.FreeEnergy.RestraintSearch(
            system,
            protocol=protocol,
            engine="somd",
            work_dir="restraint_search",
        )
        restraint_search.start()
        restraint = restraint_search.analyse(
            method="BSS", block=True, restraint_type=restraint_type
        )

    else:  # We have been supplied with a trajectory file
        logger.info("Loading trajectory file and selecting restraints...")
        # Make the output dir
        run(f"mkdir restraint_search", shell=True)
        traj = BSS.Trajectory.Trajectory(
            trajectory=trajectory_file, topology=topology_file
        )

        restraint = BSS.FreeEnergy.RestraintSearch.analyse(
            "restraint_search",
            system,
            traj,
            298 * BSS.Units.Temperature.kelvin,
            method="BSS",
            restraint_type=restraint_type,
        )

    # Save the restraint correction.
    with open("restraint_correction.txt", "w") as f:
        f.write("Restraint correction in kcal/mol:\n")
        f.write(str(restraint.getCorrection(method="numerical")))

    # Save the restraint selected
    with open("restraint.txt", "w") as f:
        f.write(restraint.toString("Somd"))
        f.write("\n")

    # Set up the stages of the bound leg.
    for stage_name, stage_options in stage_options.items():
        if stage_name not in PERTURBATION_TYPES:
            raise ValueError(
                f"Stage name {stage_name} not recognised. Must be one of {list(PERTURBATION_TYPES.keys())}"
            )

        perturbation_type = PERTURBATION_TYPES[stage_name]
        logger.info(f"Setting up {stage_name} stage...")
        protocol = BSS.Protocol.FreeEnergy(
            perturbation_type=perturbation_type,
            lam_vals=stage_options["lambda array"],
            # Runtime options are overwritten by the stage_options dictionary.
            runtime=0.5 * BSS.Units.Time.nanosecond,
        )
        stage_fe_calc = BSS.FreeEnergy.AlchemicalFreeEnergy(
            restraint.system,
            protocol,
            engine="somd",
            restraint=restraint,
            work_dir=stage_name,
            extra_options=stage_options,
        )


def run_and_analyse_abfe_bound_leg(
    work_dirs: List[str],
    scripts_dir_path: str,
) -> None:
    """
    Run and analyse the abfe bound leg calculations.

    Parameters
    ----------
    work_dirs : List[str]
        A list of the work directories for the bound leg stages.
    scripts_dir_path : str
        The path to the directory containing the slurm scripts.
    """
    # Make sure that all the required scripts are present.
    required_scripts = [
        "slurm_run.sh",
        "slurm_analysis.sh",
        "collect_results.sh",
        "collect_results.py",
    ]
    for script in required_scripts:
        if not os.path.exists(os.path.join(scripts_dir_path, script)):
            raise FileNotFoundError(
                f"Required script {script} not found in {scripts_dir_path}"
            )

    # Get the absolute paths to the scripts.
    slurm_run_script = os.path.abspath(os.path.join(scripts_dir_path, "slurm_run.sh"))
    slurm_analysis_script = os.path.abspath(
        os.path.join(scripts_dir_path, "slurm_analysis.sh")
    )
    collect_results_script = os.path.abspath(
        os.path.join(scripts_dir_path, "collect_results.sh")
    )
    collect_results_py = os.path.abspath(
        os.path.join(scripts_dir_path, "collect_results.py")
    )

    # Save the analysis job ids so that we can submit the collection script as a dependency.
    analysis_job_ids = []

    # Run the bound leg stages and submit the analysis scipts as dependencies.
    for work_dir in work_dirs:
        # Find the lambda sub-dirs
        lambda_dirs = [d for d in os.listdir(work_dir)]
        job_ids = []
        for lambda_dir in lambda_dirs:
            # Run the simulation and get the slurm job ID
            cmd = f"sbatch --chdir={os.path.join(work_dir, lambda_dir)} {slurm_run_script}"
            p = Popen(
                cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True
            )
            output = p.stdout.read()
            try:
                job_ids.append(int(output.split()[-1]))
            except ValueError:
                raise ValueError(
                    f"Could not find job ID in slurm output: {output}. Command was: {cmd}"
                    "Check that the slurm options in the submission script are correct (e.g) "
                    "the partition name."
                )

        # Submit the analysis script as a dependency of the simulations.
        cmd = f'sbatch --dependency=afterok:{":".join([str(j) for j in job_ids])} --chdir={work_dir} {slurm_analysis_script}'
        p = Popen(
            cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True
        )
        output = p.stdout.read()
        try:
            analysis_job_ids.append(int(output.split()[-1]))
        except ValueError:
            raise ValueError(
                f"Could not find job ID in slurm output: {output}. Command was: {cmd}"
                "Check that the slurm options in the submission script are correct (e.g) "
                "the partition name."
            )

    # Submit the result correction script as a dependency of the analysis scripts.
    cmd = f'sbatch --dependency=afterok:{":".join([str(j) for j in analysis_job_ids])} {collect_results_script} {collect_results_py}'
    run(cmd, shell=True)


def main() -> None:
    """Set up and run the bound leg ABFE calculation with multiple distance restraints."""

    # Parse arguments from the command line
    parser = argparse.ArgumentParser("Set up the bound leg ABFE calculation.")
    parser.add_argument(
        "-top",
        "--topology",
        type=str,
        help="The topology file.",
        default="../input/SYSTEM.prm7",
    )
    parser.add_argument(
        "-crd",
        "--coordinates",
        type=str,
        help="The coordinates file.",
        default="../input/SYSTEM.rst7",
    )
    parser.add_argument(
        "-trj",
        "--trajectory",
        type=str,
        help="The trajectory file.",
        default="../input/traj.dcd",
    )
    parser.add_argument(
        "-scr",
        "--scripts_dir_path",
        type=str,
        help="The path to the scripts directory.",
        default="../scripts/",
    )
    parser.add_argument(
        "-rt",
        "--restraint_type",
        type=str,
        help="The type of restraint to use, either 'boresch' or 'multiple_distance.'",
        default="multiple_distance",
        # Set the choices to be the keys of the restraint dictionary
        choices=["boresch", "multiple_distance"],
    )
    parser.add_argument(
        "--skip_setup",
        type=bool,
        help="Skip the setup stage and run the simulations.",
    )
    args = parser.parse_args()

    # Customise the simulation options for SOMD
    shared_extra_options = {
        # 0.1 ns per cycle, default 5 ns per simulation
        "nmoves": 25_000,
        "ncycles": 50,
        "buffered coordinates frequency": 500,
        "save coordinates": True,
        "cutoff distance": "12 * angstrom",
        "energy frequency": 100,
        "center solute": False,
        "minimise": True,
        "reaction field dielectric": 78.4,
        "hydrogen mass repartitioning factor": 3,
        "constraint": "allbonds",
    }

    stage_specific_options = {
        "restrain": {
            # From Scott input
            "lambda array": (0.000, 0.125, 0.250, 0.375, 0.500, 1.000),
            "ncycles": 50,  # 5 ns
            "timestep": "4 * femtosecond",
        },
        "discharge": {
            # From optimised MIF run
            "lambda array": (0.0, 0.291, 0.54, 0.776, 1.0),
            "ncycles": 50,  # 5 ns
            "timestep": "4 * femtosecond",
        },
        "vanish": {
            # From optimised MIF run
            "lambda array": (
                0.0,
                0.026,
                0.054,
                0.083,
                0.111,
                0.14,
                0.173,
                0.208,
                0.247,
                0.286,
                0.329,
                0.373,
                0.417,
                0.467,
                0.514,
                0.564,
                0.623,
                0.696,
                0.833,
                1.0,
            ),
            "ncycles": 80,  # 8 ns
        },
        "release_restraint": {
            # Initial guess
            "lambda array": (
                0.000,
                0.125,
                0.250,
                0.375,
                0.500,
                0.625,
                0.750,
                0.875,
                1.000,
            ),
            "ncycles": 50,  # 5 ns
            "timestep": "2 * femtosecond",
        },
    }

    # Combine the above dicts
    stage_options = {}
    for stage, options in stage_specific_options.items():
        stage_options[stage] = {**shared_extra_options, **options}

    # If using Boresch restraints, remove the release_restraint stage as this is unnecessary.
    if args.restraint_type == "boresch":
        stage_options.pop("release_restraint")

    # Set up the bound leg of the ABFE calculation with multiple distance restraints.
    if not args.skip_setup:
        setup_abfe_bound_leg(
            args.topology,
            args.coordinates,
            stage_options,
            args.trajectory,
            restraint_type=args.restraint_type,
        )

    # Run and analyse the bound leg of the ABFE calculation with multiple distance restraints using SOMD.
    run_and_analyse_abfe_bound_leg(
        work_dirs=list(stage_options.keys()), scripts_dir_path=args.scripts_dir_path
    )


if __name__ == "__main__":
    main()
