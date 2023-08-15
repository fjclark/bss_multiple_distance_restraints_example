"""
Script to collect the overall results for the bound leg, once
all individual simulations and analyses are complete.
"""

import os
from typing import List

RECOGNISED_STAGES = ["restrain", "discharge", "vanish", "release_restraint"]


class AnalysisError(Exception):
    pass


def read_dg_file(filename: str) -> float:
    """
    Read the dg.txt file and return the value of dG
    """
    with open(filename, "r") as f:
        return float(f.readlines()[-4].split(",")[0])  # kcal mol-1


def collect_overall_results() -> None:
    """
    Collect the overall results for the bound leg, once
    all individual simulations and analyses are complete.
    """
    # Figure out what legs have been run
    stages_run = []
    for d in os.listdir("."):
        if d in RECOGNISED_STAGES:
            stages_run.append(d)

    # Check that all legs have output data
    for stage in stages_run:
        if not os.path.exists(os.path.join(stage, "dg.txt")):
            raise AnalysisError(f"Could not find dg.txt in {stage}")

    # Collect the results
    results = {
        stage: read_dg_file(os.path.join(stage, "dg.txt")) for stage in stages_run
    }

    # Collect the restraint correction
    with open("restraint_correction.txt", "r") as f:
        corr = float(f.readlines()[-1].split(" ")[0])
        results["restraint_correction"] = corr

    # Reverse the sign of the release_restraint leg
    if "release_restraint" in results:
        results["release_restraint"] *= -1

    # Write the results to a file
    with open("overall_dg.txt", "w") as f:
        for stage, result in results.items():
            f.write(f"{stage} {result:.3f} kcal mol-1\n")
        # Add the overall result
        f.write(f"overall {sum(results.values()):.3f} kcal mol-1\n")


if __name__ == "__main__":
    collect_overall_results()
