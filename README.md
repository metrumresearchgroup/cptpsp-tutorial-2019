# cptpsp-tutorial-2019

### Layout:

- Model files are under `model`:
    - Voriconazole PBPK model `voriPBPK.cpp`
    - Extended voriconazole PBPK model `voriPBPK_ext.cpp`
    - MAPK QSP model `mapkQSP.cpp`
- Digitized observed data are under `data`:
    - Observed voriconazole data for adult IV `Adult_IV.csv`, Pediatric IV `Pediatric_IV.csv`, adult oral `Adult_PO.csv`, and pediatric oral `Pediatric_PO.csv` https://www.ncbi.nlm.nih.gov/pubmed/25245942
    - Virtual population data `s10vpop_pk.RDS` for MAPK QSP model https://www.ncbi.nlm.nih.gov/pubmed/28649441
    - Tissue composition data `tissue_comp_PT.csv` for tissue:plasma partition coefficient calculation via the Poulin and Theil method https://www.ncbi.nlm.nih.gov/pubmed/11782904
- Scripts for running simulations and reproducing manuscript figures are under `script`

### To reproduce manuscript figures:

1. Make the folder: `script` your working directory
2. Run `pkgSetup.R` to install required packages
3. To reproduce the plot in **Figure 1**, run `mrgsolveIntro_script.R`
3. To reproduce **Figure 3**, run `voriPBPK_script.R`
4. To reproduce **Figure 5**, run `mapkQSP_script.R`


### Voriconazole PBPK model extension:

To better capture oral voriconazole PK, an extension to the original model was added to include a mechanistic absorption model with an intestinal clearance component. The permeability was optimized and the intestinal clearance was scaled from hepatic clearance guided by relative enzyme expression data https://www.jstage.jst.go.jp/article/yakushi/123/5/123_5_369/_article. The extended model is saved as `voriPBPK_ext.cpp` and the script to run the oral voriconazole simulations is `voriPBPK_ext_script.R`.

