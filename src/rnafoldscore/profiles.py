from __future__ import annotations
import os
import numpy as np


def load_profiles(
    profiles_dir: str,
    pair_types: list[str],
    n_bins: int = 20,
) -> dict[str, np.ndarray]:
    """
    Loads interaction profiles from text files: AA.txt, AU.txt, etc.

    Returns a dictionary:
        pair_type -> numpy array of length n_bins
    """
    profiles: dict[str, np.ndarray] = {}
    for pt in pair_types:
        path = os.path.join(profiles_dir, f"{pt}.txt")
        arr = np.loadtxt(path, dtype=float)
        if arr.shape[0] != n_bins:
            raise ValueError(f"{pt}: expected {n_bins} lines, got {arr.shape[0]}")
        profiles[pt] = arr
    return profiles


def score_with_linear_interpolation(
    profile: np.ndarray,
    distance: float,
    max_dist: float = 20.0,
) -> float:
    """
    Computes a score using linear interpolation.

    The profile contains n_bins values covering the range 0..max_dist
    (e.g. 20 values for 0–20 Å). Each profile[i] is assumed to represent
    distances in the interval around i..i+1 (with a typical bin width of 1 Å).

    The score is obtained by linearly interpolating between neighboring bins.

    Distances greater than max_dist are ignored and return 0.0.
    """
    n_bins = len(profile)
    if distance < 0.0 or distance > max_dist:
        return 0.0

    bin_width = max_dist / n_bins  # typically 1.0 Å
    x = distance / bin_width       # e.g. 4.3 -> 4.3

    i0 = int(np.floor(x))
    if i0 >= n_bins - 1:
        return float(profile[n_bins - 1])

    i1 = i0 + 1
    t = x - i0  # fraction between i0 and i1 (0..1)

    return float(profile[i0] * (1.0 - t) + profile[i1] * t)
