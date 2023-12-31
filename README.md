# Example Scripts and Inputs to Run an ABFE Calculation with Multiple Distance Restraints

This uses BioSimSpace to prepare the calculation and select the restraints, and SOMD to run the 
ABFE calculations. Note that this currently only runs the bound leg of an ABFE calculation. SLURM
is used for job submission.

## Requirements

- Custom version of Sire available [here](https://github.com/fjclark/sire-openbiosim/tree/feature_permanent_mdr). Follow the instructions to install from source, making sure to install the BioSimSpace dependencies. Make sure to check out the correct branch before installing (feature_permanent_mdr).
- Custom version of BioSimSpace available [here](https://github.com/fjclark/biosimspace-openbiosim/tree/feature_mdr). Follow the instructions to install from source. Make sure to check out the correct branch before installing (feature_mdr).

## Customising Input/ Scripts

To allow the calculation to run for your system on your cluster:

- Edit the SLURM options in the  `slurm_run.sh` and `slurm_analysis.sh` scripts to match your cluster.
- Replace the `SYSTEM.rst7` and `SYSTEM.prm7` files in the input directory with your own equilibrated
  input files. Add a trajectory (of the fully interacting complex) file called "traj.dcd" to allow
  the restraint to be selected. **Ensure that your ligand is named LIG**.

## Running the Calculation

Create a calculation directory and move into it:
```
mkdir calculation
cd calculation
```
From the `calculation` directory, run:
```
python ../scripts/run_abfe_bound.py
```
To use Boresch restraints instead, run:
```
python ../scripts/run_abfe_bound.py -rt boresch
```
If you are running from a different relative path, or have renamed/ moved some of the scripts
or inputs, you will have to specify this. Run `python ../scripts/run_abfe_bound.py --help`
for details.
