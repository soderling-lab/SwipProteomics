### `SwipProteomics/analysis`

This directory contains executable scripts which perform the data analysis.

```
┌── 1_WASH-iBioID
│   │ # Analysis of WASH iBioID proteomics
│   └── 1_BioID-analysis.R
│
├── 2_SWIP-TMT
│   │ # Analysis of SWIP TMT Proteomics
│   ├── 0_PD-data-preprocess.R
│   ├── 1_MSstatsTMT-analysis.R
│   └── 2_protein-variancePartition.R
│
├── 3_Clustering
│   │ # Clustering of the spatial proteomics network 
│   ├── 1_generate-network.R
│   ├── 2_leidenalg-clustering.py
│   └── 3_post-leidenalg.R
│
├── 4_Module-Analysis
│   │ # Module-level analysis of spatial proteomics modules
│   ├── 1_module-lmerTest-analysis.R
│   ├── 2-module-variancePartition.R
│   ├── 3_module-gsea.R
│   └── 4_complex-lmerTest-analysis.R
│
└── 5_Plotting
    │ # Generate plots
    ├── 0_generate-colors.R
    ├── 1_plot-proteins.R
    ├── 2_plot-protein-profiles.R
    ├── 3_plot-module-profiles.R
    ├── 4_network-graph.R
    ├── 5_Protein-PCA.R
    ├── 6_plot-WASHC-profile.R
    └── 7_generate-module-networks.R
```
