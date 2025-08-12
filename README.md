# Description
**Diagonalgo** is a simple algorithm for the identification of [**replication fountains**](https://doi.org/10.1126/science.adj7606) in Repli-HiC data. It extracts the diagonals of chromosomal contact matrices and finds peaks that persist with increasing offsets. If a peak in the signal persists, even if it is not exactly on the antidiagonal, this is indicative of the continuous interactions characterizing a fountain. Currently, there are only functions that output the detected peaks for analysis, and functions for visualization. Further work is needed to refine the process and to find the best parameters for analysis.
Made under the supervision of Dr. Hossein Salari of the Chen Lab at Institut Curie.

# Requirements
Made in Python 3. The packages [cooler](https://github.com/open2c/cooler) and [findpeaks](https://github.com/erdogant/findpeaks) and their requirements must be installed before use. (check specific dependencies later)
