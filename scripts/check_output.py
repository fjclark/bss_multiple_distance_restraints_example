"""Script to check that all lambda windows have completed successfully."""

import pathlib
from typing import List


def check_successful_run(slurm_file: pathlib.Path) -> bool:
    """Check that the slurm file has completed successfully."""
    with open(slurm_file, "r") as f:
        for line in f:
            if "Simulation took" in line:
                return True

    return False


def check_all_runs(
    output_dir: pathlib.Path = pathlib.Path("."), slurm_output_base: str = "slurm_"
) -> List:
    """Check that all slurm files have completed successfully."""
    failed_slurm_files = []
    succesfull_slurm_files = []
    output_dir = pathlib.Path(output_dir)
    for slurm_file in output_dir.glob(f"lambda*/{slurm_output_base}*"):
        list_to_append = (
            succesfull_slurm_files
            if check_successful_run(slurm_file)
            else failed_slurm_files
        )
        list_to_append.append(slurm_file)
    if failed_slurm_files:
        print(
            "FAILURE: Some slurm files have not completed successfully (or are still running)."
        )
    else:
        print("SUCCESS: All slurm files have completed successfully.")
    print(f"Failed slurm files: {failed_slurm_files}")
    print(f"Successful slurm files: {succesfull_slurm_files}")
    return failed_slurm_files


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", type=str, default=".")
    parser.add_argument("--slurm_output_base", type=str, default="slurm_")
    args = parser.parse_args()
    check_all_runs(pathlib.Path(args.output_dir), args.slurm_output_base)
