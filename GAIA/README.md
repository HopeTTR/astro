# Gaia EDR3 Binary Star Catalog Creation

This repository contains code to create a catalog of binary star candidates from Gaia EDR3 data, based on the methodology from El-Badry et al. 2021.

## Files in this Repository

- `find_binaries_edr3.py`: Main script to identify binary star candidates
- `num_neighbors_edr3.py`: Script to identify crowded regions (clusters)
- `find_binaries_edr3_small.py`: Modified script to work with 100-row dataset
- `num_neighbors_edr3_small.py`: Modified script to work with 100-row dataset
- `find_binaries_edr3_500.py`: Modified script to work with 500-row dataset
- `num_neighbors_edr3_500.py`: Modified script to work with 500-row dataset
- `extract_subset.py`: Script to extract 100 rows from the full dataset
- `extract_subset_500.py`: Script to extract 500 rows from the full dataset

## Required Input Data

- `edr3_parallax_snr5_goodG.fits.gz`: Gaia EDR3 data with parallax SNR > 5
  - You can obtain this by running the ADQL query described in the El-Badry et al. 2021 paper
  - For testing, you can create smaller subsets using the extraction scripts

## Running the Full Analysis

The process requires three steps:

1. Obtain the input data file `edr3_parallax_snr5_goodG.fits.gz`
2. Run `num_neighbors_edr3.py` to generate `neighbor_counts_edr3_all.npz` (identifies crowded regions)
3. Run `find_binaries_edr3.py` to create the final binary catalog `binary_catalog.fits`

```bash
# Step 1: First run the ADQL query from the paper to get edr3_parallax_snr5_goodG.fits.gz

# Step 2: Count neighbors (identifies clusters)
python num_neighbors_edr3.py
# This generates neighbor_counts_edr3_all.npz

# Step 3: Find binary candidates
python find_binaries_edr3.py
# This generates binary_catalog.fits
```

## Running on Smaller Datasets

For testing or if you have limited computational resources, you can run the analysis on smaller subsets:

### 100-Row Dataset

```bash
# Create a small subset of the data (first 100 rows)
python extract_subset.py

# Count neighbors for the small dataset
python num_neighbors_edr3_small.py

# Find binary candidates in the small dataset
python find_binaries_edr3_small.py
```

### 500-Row Dataset

```bash
# Create a larger subset of the data (first 500 rows)
python extract_subset_500.py

# Count neighbors for the 500-row dataset
python num_neighbors_edr3_500.py

# Find binary candidates in the 500-row dataset
python find_binaries_edr3_500.py
```

## Results Summary

Here's a summary of the binary candidates found in our test datasets:

| Dataset | Stars | Binary Candidates | File |
|---------|-------|-------------------|------|
| 100 rows | 100 | 6 | binary_catalog_small.csv |
| 500 rows | 500 | 27 | binary_catalog_500.csv |
| Full dataset | ~64.4 million | ~1 million (estimated) | binary_catalog.fits |

## Output Format

The final catalog contains pairs of stars that are likely to be gravitationally bound binary systems. For each pair, the catalog includes:

- Coordinates (ra, dec) for both stars
- Proper motions for both stars
- Parallaxes for both stars
- G-band magnitudes for both stars
- Angular separation (arcsec)
- Physical separation (AU)
- Binary type (for the full dataset only)

## Notes

- The full analysis requires significant computational resources
- For the full dataset, the binary type classification is performed (WD-WD, WD-MS, MS-MS)
- For the small datasets, binary types are labeled as "UNKNOWN" due to missing bp-rp color information
- The 500-row dataset identified 27 binary candidates, which suggests a binary fraction of about 5.4% in the tested sample

## Reference

El-Badry et al. 2021, "A million binaries from Gaia eDR3: sample selection and validation of Gaia parallax uncertainties" 
