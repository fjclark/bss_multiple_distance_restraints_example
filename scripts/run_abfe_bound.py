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
    """
    # Load the system
    logger.info("Loading system...")
    system = BSS.IO.readMolecules([topology_file, coordinates_file])
    # Assume that the ligand is called "MOL"
    for i, mol in enumerate(system):
        if mol._getSireObject().name().value() == "MOL":
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
            method="BSS", block=True, restraint_type="multiple_distance"
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
            restraint_type="multiple_distance",
        )

    # Save the restraint correction.
    with open("restraint_correction.txt", "w") as f:
        f.write("Restraint correction in kcal/mol:\n")
        f.write(str(restraint.getCorrection(method="numerical")))

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
    slurm_run_script: str,
    slurm_analysis_script: str,
    collect_results_script: str,
) -> None:
    """
    Run and analyse the abfe bound leg calculations.

    Parameters
    ----------
    work_dirs : List[str]
        A list of the work directories for the bound leg stages.
    slurm_run_script : str
        The path to the slurm run script.
    slurm_analysis_script : str
        The path to the slurm analysis script.
    collect_results_script : str
        The path to the script to collect the results.
    """

    # Save the analysis job ids so that we can submit the collection script as a dependency.
    analysis_job_ids = []

    # Run the bound leg stages and submit the analysis scipts as dependencies.
    for work_dir in work_dirs:
        # Find the lambda sub-dirs
        lambda_dirs = [
            d for d in os.listdir(work_dir) if d in PERTURBATION_TYPES.keys()
        ]
        job_ids = []
        for lambda_dir in lambda_dirs:
            # Run the simulation and get the slurm job ID
            cmd = f"sbatch --chdir={os.path.join(work_dir, lambda_dir)} {slurm_run_script}"
            p = Popen(
                cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True
            )
            output = p.stdout.read()
            job_ids.append(int(output.split()[-1]))

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
    cmd = f'sbatch --dependency=afterok:{":".join([str(j) for j in analysis_job_ids])} {collect_results_script}'
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
        "-run",
        "--slurm_run_script",
        type=str,
        help="The slurm run script.",
        default=None,
    )
    parser.add_argument(
        "-ana",
        "--slurm_analysis_script",
        type=str,
        help="The slurm analysis script.",
        default="../scripts/slurm_analysis.sh",
    )
    parser.add_argument(
        "-res",
        "--collect_results_script",
        type=str,
        help="The script to collect the results.",
        default="../scripts/collect_results.py",
    )
    args = parser.parse_args()

    # Customise the simulation options for SOMD
    shared_extra_options = {
        # 0.1 ns per cycle, default 5 ns per simulation
        "nmoves": 25_000,
        "timestep": "4 * femtosecond",
        "ncycles": 50,
        "buffered coordinates frequency": 500,
        "save coordinates": True,
        "cutoff distance": "12 * angstrom",
        "energy frequency": 100,
        "center solute": True,
        "reaction field dielectric": 78.4,
        "hydrogen mass repartitioning factor": 3,
    }

    stage_specific_options = {
        "restrain": {
            # From Scott input
            "lambda array": [0.000, 0.125, 0.250, 0.375, 0.500, 1.000],
            "ncycles": 50,  # 5 ns
        },
        "discharge": {
            # From optimised MIF run
            "lambda array": [0.0, 0.291, 0.54, 0.776, 1.0],
            "ncycles": 50,  # 5 ns
        },
        "vanish": {
            # From optimised MIF run
            "lambda array": [
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
            ],
            "ncycles": 80,  # 8 ns
        },
        "release_restraint": {
            # Initial guess
            "lambda array": [0.000, 0.125, 0.250, 0.375, 0.500, 1.000],
            "ncycles": 50,  # 5 ns
        },
    }

    # Combine the above dicts
    stage_options = {}
    for stage, options in stage_specific_options.items():
        stage_options[stage] = {**shared_extra_options, **options}

    # Set up the bound leg of the ABFE calculation with multiple distance restraints.
    setup_abfe_bound_leg(
        args.topology, args.coordinates, stage_options, args.trajectory
    )

    # Run and analyse the bound leg of the ABFE calculation with multiple distance restraints using SOMD.
    run_and_analyse_abfe_bound_leg(
        work_dirs=list(stage_options.keys()),
        slurm_run_script=args.slurm_run_script,
        slurm_analysis_script=args.slurm_analysis_script,
        collect_results_script=args.collect_results_script,
    )


if __name__ == "__main__":
    main()
