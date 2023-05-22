# GPlag
Implementation in R and Python. Simulation and real data - Dynamic Chromatin Interactions was run in R, real data - Housing market was run in Python.

# R Dependencies
`install.packages(c("plgp", "mvtnorm", "dplyr","MASS","pracma","Matrix","condmixt","ggplot2","ggsci","patchwork","wesanderson"))`

# Python Dependencies
`pip install -r Python/requirements.txt`

# Scripts
* Figure 1: `Rscript R/Figure1.R`
* Figure 2: `Rscript R/Figure2.R`
* : Derive MTSGP results `Rscript Figure4.R`. Save the datasets from R to compare ranking results with soft-dtw and soft-dtw divergence `Python dtw_ranking_Figure4.py`
* Table 1: Synthetic data clustering. `Rscript R/Figure2.R` in section Simulation 3 -clustering/ranking
* Table 1: Dynamic Chromatin Interactions: Derive GPlag results `Rscript R/epCountsLoop_rbf.R`. Analyze the results with statistical tests `Rscript R/plot_epCountsLoop_new.R`. Data is too large to upload to Github, but can be requested if needed. 
* Table 1: Housing market. Running all metropolitan areas time series with `python/inventory_price.ipynb`. It output csv stores results of GPlag, Lead-Lag, DTW, soft-DTW, soft-DTW divergence. Then summarized the results with `Rscript R/plot_houseprice.R`
