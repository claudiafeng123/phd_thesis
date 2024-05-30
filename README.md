Pipeline for analyzing single-cell CRISPR-based experiments.

1. Read alignment (CellRanger)
2. Donor assignment (Vireo)
3. Quality control (R)
4. Guide assignment (R, Bespoke)
5. Data integration (R, Seurat)
6. Differential expression testing (R, Linear model)
7. Validation of successful knockdown.
8. Reproducibility across guides, experiments and timepoints.
9. Power estimations
10. Mean transcriptional effect of knockdown (trans effects)
11. Similarity of transcriptional effect due to knockdown (target-target correlation).
12. Similarity of perturbation profile (downstream-downstream gene correlation).
13. Deep dive into iPSC biology
14. Variation in transcriptional effect due to knockdown
15. Heritability of transcriptional effect due to knockdown
16. Intersection between heritability of transcriptional effect due to knockdown and natural genetic variation

Required software: R
Required packages: Seurat, dplyr, tidyverse, lme4, stats, data.table, igraph, doParallel
