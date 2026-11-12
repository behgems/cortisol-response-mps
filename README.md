# cortisol-response-mps
Code for 'Development and cross-tissue validation of a machine learning-based methylation profile score for the cortisol response to stress' (Balfour et al., 2026).
https://doi.org/10.1016/j.pnpbp.2026.111992

## License and citation
Code is released under the MIT License (see LICENSE). The MPS model object is released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Please cite the paper above if this code or model (MPS) are used.

## Contact
Contact:
david.balfour@flinders.edu.au
sarah.cohenwoods@flinders.edu.au

## Applying the MPS to new data
MPS model object and script for applying it to new data are currently unavailable. Will be uploaded shortly.

## Development pipeline
Script for development pipeline is ./scripts/run.R.

Environment with required software can be created using:
mamba env create -f env.yml
renv::restore()

Required data:
From Chybowska et al. (2025) (https://doi.org/10.1038/s41467-025-58357-6):
./data/Predictors_long.csv

From Miller et al. (2025) (https://doi.org/10.1038/s41598-025-09299-y):
./data/methyl/champagne/beta_no_gap_probes.rda
./data/methyl/champagne/DataAcuteStressMeth.xlsx

From Houtepen et al. (2016) (https://doi.org/10.1038/ncomms10967):
./data/GSE77445/GSE77445_series_matrix.txt

From Hearn et al. (2020) (https://doi.org/10.1186/s12920-020-0670-9):
./data/GSE119078/GSE119078_series_matrix.txt.gz

From Provençal et al. (2020) (https://doi.org/10.1073/pnas.1820842116):
./data/provencal_2020/pnas_1820842116_sd01_S7.xlsx
./data/provencal_2020/pnas_1820842116_sd01_S8.xlsx