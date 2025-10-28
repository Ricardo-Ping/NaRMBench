import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import glob

thresholds_data = pd.read_csv("best_thresholds_precision_recall.txt", sep='\t')
files = glob.glob("./*_m6A_with_labels_converted.txt")
labels=["Xron", "Dorado","SingleMod","TandemMod-retrain","m6Anet",'Epinano-retrain']
total_sites = []
verified_sites = []
precision_values = []
recall_values = []
f1_scores = []

for label in labels:
    row = thresholds_data[thresholds_data['Label'] == label]
    if not row.empty:
        precision = row['Precision'].values[0]
        recall = row['Recall'].values[0]
        f1 = row['F1 Score'].values[0]
        
        precision_values.append(precision)
        recall_values.append(recall)
        f1_scores.append(f1)


    file_path = f"./{label}_m6A_with_labels_converted.txt"
    with open(file_path, 'r') as f:
        data = [line.strip().split('\t') for line in f]

    third_column_values = [float(row[2]) for row in data]

sorted_indices = np.argsort(f1_scores)[::-1]

labels = [labels[i] for i in sorted_indices]
precision_values = [precision_values[i] for i in sorted_indices]
recall_values = [recall_values[i] for i in sorted_indices]
f1_scores = [f1_scores[i] for i in sorted_indices]
fig_width = 60 * 0.03937  
fig_height = 40 * 0.03937 
fig, ax2 = plt.subplots(figsize=(fig_width, fig_height))
index = np.arange(len(labels))
line1, = ax2.plot(index, precision_values, label='Precision', color='green', marker='^', linestyle='--', linewidth=0.8,markersize=3)
line2, = ax2.plot(index, recall_values, label='Recall', color='orange', marker='s', linestyle='--', linewidth=0.8,markersize=3)
line3, = ax2.plot(index, f1_scores, label='F1 score', color='brown', marker='o', linestyle='-', linewidth=0.8,markersize=3)
ax2.set_xticks(index)
ax2.set_xticklabels(labels, rotation=90, fontsize=4.5)  
ax2.set_ylim(0, 1)
ax2.tick_params(axis='y', labelsize=4)
handles, labels_legend = [], []
line_handles, line_labels = ax2.get_legend_handles_labels()
handles.extend(line_handles)
labels_legend.extend(line_labels)
plt.tight_layout(pad=0.1)
plt.savefig("./F1score_curve_m6A.pdf")
