from __future__ import annotations
from dataclasses import dataclass
from typing import List, Optional


@dataclass(frozen=True)
class C3Atom:
    chain_id: str
    res_seq: int          # residue sequence number in PDB
    ins_code: str         # insertion code (often empty)
    nuc: str              # A/U/C/G
    x: float
    y: float
    z: float


def resname_to_nuc(resname: str) -> Optional[str]:
    """
    Converts a residue name to a nucleotide base: A/U/C/G.

    PDB files may store nucleotides as single-letter names ('A')
    or three-letter names ('ADE'), including modified residues.
    This function performs a simple "best effort" conversion.
    """
    r = resname.strip().upper()

    # Common standard names
    if r in {"A", "ADE"}:
        return "A"
    if r in {"U", "URI"}:
        return "U"
    if r in {"C", "CYT"}:
        return "C"
    if r in {"G", "GUA"}:
        return "G"

    # Modified nucleotides often start with the base letter
    if r.startswith("A"):
        return "A"
    if r.startswith("U"):
        return "U"
    if r.startswith("C"):
        return "C"
    if r.startswith("G"):
        return "G"

    return None


def parse_pdb_c3prime(pdb_path: str, atom_name: str = "C3'") -> List[C3Atom]:
    """
    Reads a PDB file and returns a list of C3' atoms for RNA residues.

    The function processes ATOM/HETATM records and parses
    fixed-width PDB columns.

    Note:
        This is a minimal parser sufficient for the current project.
        It can be extended later if more complex PDB files are needed.
    """
    atoms: List[C3Atom] = []

    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            name = line[12:16].strip()
            if name != atom_name:
                continue

            resname = line[17:20].strip()
            nuc = resname_to_nuc(resname)
            if nuc is None:
                continue

            chain_id = (line[21] or " ").strip() or "_"

            # Residue sequence number
            try:
                res_seq = int(line[22:26])
            except ValueError:
                continue

            ins_code = (line[26] or " ").strip()

            # Atomic coordinates
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            atoms.append(C3Atom(chain_id, res_seq, ins_code, nuc, x, y, z))

    return atoms

import math
from typing import Iterator


def euclidean_distance(a: C3Atom, b: C3Atom) -> float:
    """
    Computes the Euclidean distance between two C3' atoms.
    """
    dx = a.x - b.x
    dy = a.y - b.y
    dz = a.z - b.z
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def iter_pairs_intrachain(
    atoms: list[C3Atom],
    min_seq_sep: int = 4,
) -> Iterator[tuple[C3Atom, C3Atom]]:
    """
    Generates (i, j) atom pairs within the same chain (intrachain),
    such that j >= i + min_seq_sep (the i and i+4 rule).

    IMPORTANT:
    - The residue order follows the order in the file (after sorting)
    - The insertion code is taken into account during sorting
    """
    # Group atoms by chain
    chains: dict[str, list[C3Atom]] = {}
    for atom in atoms:
        chains.setdefault(atom.chain_id, []).append(atom)

    for chain_atoms in chains.values():
        # Sort residues within each chain
        chain_atoms_sorted = sorted(
            chain_atoms,
            key=lambda a: (a.res_seq, a.ins_code)
        )

        n = len(chain_atoms_sorted)
        for i in range(n):
            for j in range(i + min_seq_sep, n):
                yield chain_atoms_sorted[i], chain_atoms_sorted[j]
