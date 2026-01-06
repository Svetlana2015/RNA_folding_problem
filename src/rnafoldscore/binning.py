from __future__ import annotations
import numpy as np
from typing import Optional


def make_bin_edges(max_dist: float = 20.0, n_bins: int = 20) -> np.ndarray:
    """
    Creates bin edges from 0 to max_dist.

    Returns an array of length n_bins + 1.
    """
    if max_dist <= 0:
        raise ValueError("max_dist must be > 0")
    if n_bins <= 0:
        raise ValueError("n_bins must be > 0")

    return np.linspace(0.0, float(max_dist), int(n_bins) + 1)


def distance_to_bin(
    distance: float,
    edges: np.ndarray,
) -> Optional[int]:
    """
    Converts a distance value to a bin index.

    Returns None if the distance is outside the bin range.
    """
    if distance < edges[0] or distance > edges[-1]:
        return None

    # np.digitize returns indices in the range 1..n_bins
    idx = int(np.digitize([distance], edges, right=True)[0]) - 1

    # Guard against falling exactly on the upper edge
    if idx == len(edges) - 1:
        idx = len(edges) - 2

    if idx < 0 or idx >= len(edges) - 1:
        return None

    return idx
