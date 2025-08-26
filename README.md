# Description
**Diagonalgo** is a simple tool for the identification of [**replication fountains**](https://doi.org/10.1126/science.adj7606) and other structures in Repli-HiC data. It extracts the diagonals of chromosomal contact matrices and finds peaks that persist with increasing offsets. If a peak in the signal persists, even if it is not exactly on the antidiagonal, this is indicative of the continuous interactions characterizing a fountain. Further work is needed to refine the process and to find the best parameters for analysis.

# Installation
```
git clone https://github.com/burgundy-snail/diagonalgo.git
cd diagonalgo
pip install -r requirements.txt
```

## Requirements
Made in Python 3.13.5. The minimum version should be 3.10 although I need to confirm this and make adjustments. See [requirements.txt](https://github.com/burgundy-snail/diagonalgo/blob/main/requirements.txt) for package dependencies.

# Features and usage
Diagonalgo currently is able to **extract** diagonals, **preprocess** them with denoising and interpolation, **detect peaks** in diagonal signals, and **visualize** the results.

Currently, it is not packaged. An example of usage can be found in [analysis.py](https://github.com/burgundy-snail/diagonalgo/blob/main/analysis.py) - proper examples TBA.

## Extraction
1. Create a `DiagCool` object - identical to initializing a normal Cooler.
2. Generate a dataframe using `get_aligned` for a specified chromosome and range of offsets. A minimum of at least 40kbp is reccomended; offset is by index so make sure to divide by bin width.
See documentation for other features of `DiagCool` including preliminary visualization of raw diagonals.

## Preprocessing
Preprocessing is handled as arguments of peakfinder.
- **Denoising:** `savgol` or `gaussian` or none.
- **Interpolation:** see documentation of findpeaks.

## Processing
See documentation

## Visualization
Two options: line graph with marked peaks or dotplot of distance vs offset.

# Examples
TBA :(

# Acknowledgements
Created by [Helen Zhang](https://github.com/burgundy-snail/) under the supervision of Dr. Hossein Salari of the Chen Lab at Institut Curie.

### References
TBA

### Mantainence
- Helen Zhang - will mantain when able.
