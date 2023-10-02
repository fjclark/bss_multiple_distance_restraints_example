"""Collect overall statistics based on a set of calculations."""

from typing import List, Dict

import glob
import numpy as np
import os
import scipy.stats as stats


def get_dirs(dir_start: str) -> List[str]:
    """Glob directories from dir_start."""
    # get the directories
    dirs = glob.glob(dir_start + "*")
    # Make sure they're all directories and not files
    dirs = [directory for directory in dirs if os.path.isdir(directory)]
    return dirs


def read_results(output_dir: str) -> Dict[str, float]:
    """Read results from output_dir."""
    results = {}
    with open(output_dir + "/overall_dg.txt", "r") as f:
        for line in f:
            quantity = line.split()[0]
            value = float(line.split()[1])
            results[quantity] = value
    return results


def get_stats(dir_start: str) -> Dict[str, float]:
    """Get statistics from a set of calculations."""
    # Get the results dicts
    dirs = get_dirs(dir_start)
    results = []
    for directory in dirs:
        results.append(read_results(directory))

    # Get the statistics - specifically mean and 95 % CI
    overall_stats = {}
    for quantity in results[0]:
        overall_stats[quantity] = {}
        values = np.array([result[quantity] for result in results])
        overall_stats[quantity]["mean"] = values.mean()
        if values.std() == 0:
            overall_stats[quantity]["ci"] = 0
        else:
            overall_stats[quantity]["ci"] = (
                stats.t.interval(
                    0.95, len(values) - 1, loc=np.mean(values), scale=stats.sem(values)
                )[1]
                - values.mean()
            )

    return overall_stats


def write_stats(
    overall_stats: Dict[str, float], output_dir: str, output_name: str
) -> None:
    """Write the statistics to output_dir."""
    with open(output_dir + output_name + "_stats.txt", "w") as f:
        for quantity in overall_stats:
            f.write(
                f"{quantity} {overall_stats[quantity]['mean']:.2f} +/- {overall_stats[quantity]['ci']:.2f} kcal mol-1\n"
            )


def main():
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-d",
        "--dir_start",
        type=str,
        help="Common start to directories to get statistics from",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        help="Output directory to write statistics to",
        default="./",
    )

    args = parser.parse_args()

    overall_stats = get_stats(args.dir_start)
    write_stats(overall_stats, args.output_dir, args.dir_start)


if __name__ == "__main__":
    main()
