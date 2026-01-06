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

