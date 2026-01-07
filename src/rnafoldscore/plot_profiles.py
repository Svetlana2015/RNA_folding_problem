from __future__ import annotations
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

from rnafoldscore.constants import PAIR_TYPES
from rnafoldscore.profiles import load_profiles


def main():
    p = argparse.ArgumentParser(
        prog="rnafoldscore.plot_profiles",
        description="Plot interaction profiles (score vs distance).",
    )
    p.add_argument("--profiles", required=True, help="Directory with profile txt files")
    p.add_argument("--out", required=True, help="Output directory for plots")
    p.add_argument("--max-dist", type=float, default=20.0)
    p.add_argument("--bins", type=int, default=20)
    args = p.parse_args()

    profiles = load_profiles(args.profiles, PAIR_TYPES, n_bins=args.bins)

    os.makedirs(args.out, exist_ok=True)

    bin_width = args.max_dist / args.bins
    distances = np.linspace(
        bin_width / 2,
        args.max_dist - bin_width / 2,
        args.bins,
    )

    for pt, values in profiles.items():
        plt.figure()
        plt.plot(distances, values, marker="o")
        plt.xlabel("Distance (Å)")
        plt.ylabel("Score u(r)")
        plt.title(f"Interaction profile {pt}")
        plt.grid(True)
        out_path = os.path.join(args.out, f"{pt}.png")
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close()

    print("Plots written to", args.out)


if __name__ == "__main__":
    main()
