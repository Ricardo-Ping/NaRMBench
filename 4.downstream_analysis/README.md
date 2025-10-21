# Downstream Analysis Scripts
After model evaluation, the following scripts can be used for visualization and downstream analyses:
- [Fig2-3.py](Fig2-3.py) ：Generates ROC and PR curves, and computes quantification accuracy.
- [Fig4.r](Fig4.r) ：Analyzes the distribution of predicted modification sites across transcript regions and compares WT vs. KO (or unmodified IVT) samples.
- [Suppl_Fig5_optimal_cutoff.py](Suppl_Fig5_optimal_cutoff.py) : Determines the optimal classification cut-off for each model.
- [Fig5.py](Fig5.py) : Computes replicate overlap (Jaccard index) and correlation (Pearson’s r), and evaluates sequencing depth bias, modification level bias, and motif bias.
- [Suppl_Fig6_all_F1.py](Suppl_Fig6_all_F1.py) : Plots the F1 score curve.
- [Fig6.r](Fig6.r) ：Generates radar plots for model performance comparison.
- [Suppl_Fig12_ambiguous.r](Suppl_Fig12_ambiguous.r) : Examines the proportions of ambiguous sites and compares signal-derived features between true positive and false positive sites identified by TandemMod.
- [Suppl_Fig15_consistency.py](Suppl_Fig15_consistency.py) : Creates a heatmap of site overlap ratios and computes Pearson’s correlation of predicted modification levels across different quantification tools.
