# cptpsp-tutorial-2019

### Layout:

- Model files are under `model`:
    - Voriconazole PBPK model `voriPBPK.cpp`
    - MAPK QSP model `mapkQSP.cpp`
- Digitized observed data are under `data`:
    - Observed voriconazole data for adult IV `Adult_IV.csv` and Pediatric IV `Pediatric_IV.csv` https://www.ncbi.nlm.nih.gov/pubmed/25245942
    - Virtual population data `s10vpop_pk.RDS` for MAPK QSP model https://www.ncbi.nlm.nih.gov/pubmed/28649441
    - Tissue composition data `tissue_comp_PT.csv` for tissue:plasma partition coefficient calculation via the Poulin and Theil method https://www.ncbi.nlm.nih.gov/pubmed/11782904
- Scripts for running simulations and reproducing manuscript figures are under `script`

### To reproduce manuscript figures:

1. Make the folder: `script` your working directory
2. Run `pkgSetup.R` to install required packages
3. To reproduce the plot in **Figure 1**, run `mrgsolveIntro_script.R`
3. To reproduce **Figure 3**, run `voriPBPK_script.R`
4. To reproduce **Figure 5**, run `mapkQSP_script.R`