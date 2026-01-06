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


def distance_to_bin(distance: float, max_dist: float = 20.0, n_bins: int = 20) -> Optional[int]:
    """
     Converts a distance value to a bin index.

    Returns None if the distance is outside the bin range.
    """
    if distance < 0.0 or distance > max_dist:
        return None

    bin_width = max_dist / n_bins  # например 20/20 = 1.0
    idx = int(distance // bin_width)

    # distance == max_dist даст idx == n_bins, отправляем в последний бин
    if idx == n_bins:
        idx = n_bins - 1

    return idx