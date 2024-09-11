# GPlag
GPlag algorithm has written in both R (`R/GPlag.R`) and Python (`Python/GPlag.py`). Simulation and real data - Dynamic Chromatin Interactions was run in R.

# R Dependencies
`install.packages(c("plgp", "mvtnorm", "dplyr","MASS","pracma","Matrix","condmixt","ggplot2","ggsci","patchwork","wesanderson"))`

# Python Dependencies
`pip install -r Python/requirements.txt`

# Scripts
* Figure 1: Example genes which could be plotted via `Rscript R/plot_epCountsLoop_new.R`.
* Figure 2: `Rscript R/Figure2.R`. Simulated datasets are saved from R to compare with Lead-Lag in `python/lead-lag.ipynb`
* Figure 3: `Rscript R/Figure3.R`.
* Figure 4: `Rscript R/multits.R`.
* Figure 5: `Rscript R/plot_epCountsLoop_new.R` after run `Rscript R/epCountsLoop_rbf.R`.
* Figure S1: `Rscript R/Figure2.R` under "MLE of sigma2 and b".
* Figure S2: `Rscript R/mle_b.R`.
* Figure S3: `Rscript R/Figure3.R`.
* Figure S4: (A)`Rscript R/table2_rbf.R`, (B)`Rscript R/table2_matern.R`.
* Figure S5-7: `Rscript R/sample_path.R`.
* Table 2: Synthetic data clustering. `Rscript R/Figure3.R` under "Simulation 3 -clustering/ranking". Simulated datasets are saved from R to compare with soft-DTW divergence in `python/softdtw_div.ipynb`
* Table 2: Dynamic Chromatin Interactions: Derive GPlag results `Rscript R/epCountsLoop_rbf.R`. The soft-DTW divergence result can be found in `python/softdtw_div.ipynb`. Results from other baseline methods and analysis with statistical tests are included in `Rscript R/plot_epCountsLoop_new.R`. The data could be accessed at Reed, Kathleen SM, et al. Cell reports 41.5 (2022)(PMID: 36323252).