from __future__ import annotations
import argparse

from rnafoldscore.constants import PAIR_TYPES
from rnafoldscore.pdb_utils import (
    parse_pdb_c3prime,
    iter_pairs_intrachain,
    euclidean_distance,
)
from rnafoldscore.profiles import load_profiles, score_with_linear_interpolation


def pair_type(a, b) -> str | None:
    pt = "".join(sorted([a.nuc, b.nuc]))
    return pt if pt in PAIR_TYPES else None


def main():
    p = argparse.ArgumentParser(
        prog="rnafoldscore.score_structure",
        description="Score a PDB structure using trained profiles (linear interpolation).",
    )
    p.add_argument("--pdb", required=True, help="Input PDB file to score")
    p.add_argument("--profiles", required=True, help="Directory containing trained profiles")
    p.add_argument("--max-dist", type=float, default=20.0)
    p.add_argument("--bins", type=int, default=20)
    p.add_argument("--min-seq-sep", type=int, default=4)
    args = p.parse_args()

    profiles = load_profiles(args.profiles, PAIR_TYPES, n_bins=args.bins)

    atoms = parse_pdb_c3prime(args.pdb)

    total = 0.0
    n_terms = 0

    for a, b in iter_pairs_intrachain(atoms, min_seq_sep=args.min_seq_sep):
        d = euclidean_distance(a, b)

        pt = pair_type(a, b)
        if pt is None:
            continue

        s = score_with_linear_interpolation(profiles[pt], d, max_dist=args.max_dist)

        total += s
        n_terms += 1

    print(f"Score (sum): {total:.6f}")
    print(f"Number of scored pairs: {n_terms}")


if __name__ == "__main__":
    main()
