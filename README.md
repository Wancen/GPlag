# GPlag
Implementation in R and Python. Simulation and real data - Dynamic Chromatin Interactions was run in R, real data - Housing market was run in Python.

# R Dependencies
`install.packages(c("plgp", "mvtnorm", "dplyr","MASS","pracma","Matrix","condmixt","ggplot2","ggsci","patchwork","wesanderson"))`

# Python Dependencies
`pip install -r Python/requirements.txt`

# Scripts
* Figure 1: `Rscript R/Figure1.R`. Simulated datasets were saved from R to compare with Lead-Lag in `python/lead-lag.ipynb`
* Figure 2: `Rscript R/Figure2.R`
* Table 1: Synthetic data clustering. `Rscript R/Figure2.R` under "Simulation 3 -clustering/ranking". Simulated datasets were saved from R to compare with soft-DTW divergence in `python/softdtw_div.ipynb`
* Table 1: Dynamic Chromatin Interactions: Derive GPlag results `Rscript R/epCountsLoop_rbf.R`. soft-DTW divergence result is in `python/softdtw_div.ipynb`. Results from other baseline methods and analysis with statistical tests is in `Rscript R/plot_epCountsLoop_new.R`. Data is too large to upload to Github, but can be requested if needed. 
* Table 1: Housing market. Running all metropolitan areas time series and all methods with `python/inventory_price.ipynb`. The output CSV stores results of GPlag, Lead-Lag, DTW, soft-DTW, soft-DTW divergence. Then summarized the results with `Rscript R/plot_houseprice.R`
