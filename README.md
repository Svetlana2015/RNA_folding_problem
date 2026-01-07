# RNA Objective Function - Distance-Based Statistical Potential
## Project context

### Master 2 GENIOMHE - Bioinformatique de la structure de l’ARN

This project implements a distance-based statistical potential for RNA folding, based on inter-atomic distances between C3’ atoms extracted from experimentally determined RNA 3D structures (PDB).

The goal is to build an objective (scoring) function that estimates the Gibbs free energy of an RNA conformation and can be used to rank predicted RNA structures.

## Method overview

For a given RNA 3D structure:

- Only C3’ atoms are considered.

- Only intrachain residue pairs are used.

- Only pairs separated by at least 4 positions in the sequence (i, i+4, i+5, …).

- Distances are binned in 20 bins between 0 and 20 Å.

- Ten base-pair types are considered:

``` text
AA, AU, AC, AG,
UU, UC, UG,
CC, CG,
GG
```

###  Training (statistical potential)

For each base-pair type and each distance bin:

- **Observed frequencies** are computed from training PDB structures:

$$
f^{\mathrm{OBS}}_{ij}(r) = \frac{N_{ij}(r)}{N_{ij}}
$$

- A **reference distribution (XX)** is computed by ignoring base identities:

$$
f^{\mathrm{REF}}_{XX}(r) = \frac{N_{XX}(r)}{N_{XX}}
$$

- The **interaction potential** is defined as:

$$
u_{ij}(r) = -\log \left( \frac{f^{\mathrm{OBS}}_{ij}(r)}{f^{\mathrm{REF}}_{XX}(r)} \right)
$$

- Values are capped at a maximum of **10**.

Each interaction profile is saved as a text file with **20 values**
(one per distance bin).

## Scoring a structure

For a new RNA structure:

1. **All valid intrachain pairs** ($j \ge i + 4$) are extracted.
2. **Distances** between C3' atoms are computed.
3. For each distance, the corresponding score is obtained using
   **linear interpolation** between neighboring bins.
4. The **total score** is the **sum of all interaction scores**,
   acting as a proxy for the Gibbs free energy.

Lower scores correspond to more favorable (native-like) conformations.

## Repository structure

```bash
RNA_project/
├── src/rnafoldscore/           # Python package
│   ├── __init__.py
│   ├── constants.py            # Base-pair definitions
│   ├── pdb_utils.py            # PDB parsing, C3' extraction, distances
│   ├── binning.py              # Distance binning
│   ├── profiles.py             # Profile loading & interpolation
│   ├── train.py                # Training script
│   ├── plot_profiles.py        # Plotting script
│   └── score_structure.py      # Scoring script
│
├── data/pdb/                   # Training PDB files (not versioned)
├── examples/                   # Example PDB files
├── outputs/
│   ├── profiles/               # Trained profiles (AA.txt, AU.txt, ...)
│   └── plots/                  # Profile plots
│
├── pyproject.toml
├── README.md
└── .gitignore
```

## Installation

Python **3.10+** is required.

```bash
python -m venv .venv
.venv\Scripts\activate    # Windows
pip install -e .
```
## Usage

### 1. Training interaction profiles

``` bash
python -m rnafoldscore.train --pdb-dir data\pdb --out-dir outputs\profiles
```

This generates 10 profile files, each containing 20 distance bins.

### 2. Plotting profiles

``` bash
python -m rnafoldscore.plot_profiles --profiles outputs\profiles --out outputs\plots
```
This produces one plot per interaction type (AA.png, AU.png, …).

### 3. Scoring an RNA structure

``` bash
python -m rnafoldscore.score_structure --pdb examples\example.pdb --profiles outputs\profiles
```
Output:

- Total score (sum of interaction terms)

- Number of scored residue pairs

## Command-line interface

All scripts provide a `--help` option:

```bash
python -m rnafoldscore.train --help
python -m rnafoldscore.plot_profiles --help
python -m rnafoldscore.score_structure --help
```

### Good practices

- Modular code with reusable functions

- No code duplication

- Version control using Git

- Clear CLI interfaces

- Reproducible environment

- Documented repository

### Limitations

- This implementation is intended for educational purposes.

- Results obtained from very small datasets are not statistically meaningful.

- Realistic profiles require a sufficiently large and diverse training set of RNA structures.

### Author

Student project - Master 2 GENIOMHE
Bioinformatique de la structure de l’ARN
