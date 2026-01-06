import argparse
import os
import numpy as np
from collections import defaultdict

from rnafoldscore.constants import PAIR_TYPES
from rnafoldscore.pdb_utils import (
    parse_pdb_c3prime,
    iter_pairs_intrachain,
    euclidean_distance,
)
from rnafoldscore.binning import distance_to_bin


def pair_type(a, b) -> str | None:
    """
    Returns the pair type (AA, AU, ..., GG).
    Order does not matter: AU == UA.
    """
    pair = "".join(sorted([a.nuc, b.nuc]))
    if pair in PAIR_TYPES:
        return pair
    return None


def main():
    p = argparse.ArgumentParser(
        prog="rnafoldscore.train",
        description="Train distance-based interaction profiles from PDB files.",
    )
    p.add_argument("--pdb-dir", required=True, help="Directory with training PDB files")
    p.add_argument("--out-dir", required=True, help="Output directory for profiles")
    p.add_argument("--max-dist", type=float, default=20.0)
    p.add_argument("--bins", type=int, default=20)
    p.add_argument("--min-seq-sep", type=int, default=4)
    p.add_argument("--clip-max", type=float, default=10.0)
    args = p.parse_args()

    # Counters
    counts = {pt: np.zeros(args.bins, dtype=float) for pt in PAIR_TYPES}
    totals = {pt: 0.0 for pt in PAIR_TYPES}

    xx_counts = np.zeros(args.bins, dtype=float)
    xx_total = 0.0

    pdb_files = [
        os.path.join(args.pdb_dir, f)
        for f in os.listdir(args.pdb_dir)
        if f.lower().endswith(".pdb")
    ]

    if not pdb_files:
        raise RuntimeError("No PDB files found")

    for pdb in pdb_files:
        atoms = parse_pdb_c3prime(pdb)

        for a, b in iter_pairs_intrachain(atoms, min_seq_sep=args.min_seq_sep):
            d = euclidean_distance(a, b)
            bin_idx = distance_to_bin(d, args.max_dist, args.bins)
            if bin_idx is None:
                continue

            pt = pair_type(a, b)
            if pt is None:
                continue

            counts[pt][bin_idx] += 1.0
            totals[pt] += 1.0

            xx_counts[bin_idx] += 1.0
            xx_total += 1.0

    # Frequencies
    f_ref = xx_counts / xx_total

    os.makedirs(args.out_dir, exist_ok=True)

    for pt in PAIR_TYPES:
        if totals[pt] > 0:
            f_obs = counts[pt] / totals[pt]
        else:
            f_obs = np.zeros(args.bins)

        with np.errstate(divide="ignore"):
            u = -np.log((f_obs + 1e-12) / (f_ref + 1e-12))

        u = np.minimum(u, args.clip_max)

        out_path = os.path.join(args.out_dir, f"{pt}.txt")
        np.savetxt(out_path, u, fmt="%.6f")

    print("Training finished. Profiles written to", args.out_dir)


if __name__ == "__main__":
    main()
