# Downstream Analysis Scripts
After model evaluation, the following scripts can be used for visualization and downstream analyses
## Data Input Format
Each script requires model output data in the following format:
```sh
$ {model}_m6A_with_labels_converted.txt
#pos	label	prob	motif	depth	ratio
#chr1_1013990	0	0.010309278350515464	GGACC	97	0
```
Example data are provided in the example_data/ directory.
## Script Descriptions
- [Fig2-3.py](Fig2-3.py) ：Generates ROC and PR curves, and computes quantification accuracy.
- [Fig4.r](Fig4.r) ：Analyzes the distribution of predicted modification sites across transcript regions and compares WT vs. KO (or unmodified IVT) samples.
- [Suppl_Fig5_optimal_cutoff.py](Suppl_Fig5_optimal_cutoff.py) : Determines the optimal classification cut-off for each model.
- [Fig5.py](Fig5.py) : Computes replicate overlap (Jaccard index) and correlation (Pearson’s r), and evaluates sequencing depth bias, modification level bias, and motif bias.
- [Suppl_Fig6_all_F1.py](Suppl_Fig6_all_F1.py) : Plots the F1 score curve.
- [Fig6.r](Fig6.r) ：Generates radar plots for model performance comparison.
- [Suppl_Fig15_consistency.py](Suppl_Fig15_consistency.py) : Creates a heatmap of site overlap ratios and computes Pearson’s correlation of predicted modification levels across different quantification tools.
