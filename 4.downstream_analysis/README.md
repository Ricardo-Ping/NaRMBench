# Downstream Analysis
After model evaluation, the following scripts can be used for visualize the downstream analyses
## Input Data Format
The scripts below requires model output data in the following format (see example date in `example_data/` directory):
```sh
pos	label	prob	motif	depth	ratio
chr1_1013990	0	0.0103	GGACC	97	0
```
- `label`: 0 or 1 representing whether the site is a positive site; `prob`: Model-predicted probability; `ratio`: Model-predicted ratio of positive sites

## Script Descriptions
- Generate ROC and PR curves, and quantification accuracy: [Script](Fig2-3.py)
- Generate meta-gene plots, and WT vs. KO (or unmodified IVT) comparisons: [Script](Fig4.r)
- Determine the optimal classification cut-off for each tool: [Script](Suppl_Fig5_optimal_cutoff.py)
- Calculate replicates overlap and correlation, and evaluates sequencing depth bias, modification level bias, and motif bias: [Script](Fig5.py)
- Generate F1 score curve: [Script](Suppl_Fig6_all_F1.py)
- Generate radar plots for performance comparison: [Script](Fig6.r)
- Generate site overlap(heatmaps) and modification level correlation(scatter plots) across different tools: [Script](Suppl_Fig15_consistency.py)
