# Study Repository

This repository demonstrates matrix-displacement analysis for simple trusses.

## Files
- `king_post_analysis.py` — analyzes a king-post roof truss.
- `pratt_truss_analysis.py` — analyzes a six-panel Pratt through-truss bridge.

## Requirements
The examples require `numpy` and `matplotlib`.
Install them using:

```bash
pip install -r requirements.txt
```

## Usage
Run either analysis script with Python 3:

```bash
python3 king_post_analysis.py
python3 pratt_truss_analysis.py
```

Each program prints joint displacements, support reactions and axial forces and
shows a plot comparing undeformed and deformed shapes.
