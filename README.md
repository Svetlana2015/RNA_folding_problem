# RNA Objective Function — Distance-Based Statistical Potential
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

Training (statistical potential)

For each base-pair type and each distance bin:

Observed frequencies are computed from training PDB structures:

𝑓
𝑖
𝑗
𝑂
𝐵
𝑆
(
𝑟
)
=
𝑁
𝑖
𝑗
(
𝑟
)
𝑁
𝑖
𝑗
f
ij
OBS
	​

(r)=
N
ij
	​

N
ij
	​

(r)
	​


A reference distribution (XX) is computed by ignoring base identities:

𝑓
𝑋
𝑋
𝑅
𝐸
𝐹
(
𝑟
)
=
𝑁
𝑋
𝑋
(
𝑟
)
𝑁
𝑋
𝑋
f
XX
REF
	​

(r)=
N
XX
	​

N
XX
	​

(r)
	​













# RNA_folding_problem

# RNA Objective Function

Проект: статистический потенциал для РНК по расстояниям C3' (обучение профилей и скоринг PDB).

## Installation
``` bash
python -m venv .venv
.venv\Scripts\activate
pip install -e
```

## Usage

Train:

``` bash
python -m rnafoldscore.train --pdb-dir data\pdb --out-dir outputs\profiles
```

Plot:
``` bash
python -m rnafoldscore.plot_profiles --profiles outputs\profiles --out outputs\plots
```

Score:
``` bash
python -m rnafoldscore.score_structure --pdb examples\example.pdb --profiles outputs\profiles
```

