import numpy as np
import math
from sklearn import metrics
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import glob
import os
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr
sns.set(style="white")

files = glob.glob("./*_m6A_with_labels_converted.txt")
colors = ["#2e3792", "#8e5aa2", "#ef1fff", "#f6c365", "#9fcc62","#6affb9"]
labels = [os.path.basename(file).split("_m6A_with_labels_converted.txt")[0] for file in files]

##only sites covered by all tools for ROC curve / PRAUC curve
##ROC curve / PRAUC curve
plt.figure(figsize=[6.4,6.4])
for idx, file in enumerate(files):
    Y, X = [], []
    with open(file, "r") as f:
        for i in f:
            ele = i.rstrip().split()
            Y.append(int(ele[1])) 
            if "PsiNanopore" in file or "ELIGOS" in file:
                p_value = float(ele[2])  
                if p_value == 0:
                    X.append(1000)  
                else:
                    X.append(-math.log10(p_value))
            else:
                X.append(float(ele[2])) 

    y = np.array(Y)
    x = np.array(X)
    
    fpr, tpr, _ = metrics.roc_curve(y, x)
    auc = metrics.auc(fpr, tpr)
    
    plt.plot(fpr, tpr, color=colors[idx % len(colors)], label=f'{labels[idx]} (AUC = {auc:.2f})')

plt.plot([0, 1], [0, 1], 'k--')
plt.title('ROC curve')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(fontsize="small", frameon=False)
plt.savefig("./ROC_m6A.pdf")

plt.figure(figsize=[6.4,6.4])
for idx, file in enumerate(files):
    Y, X = [], []
    with open(file, "r") as f:
        for i in f:
            ele = i.rstrip().split()
            Y.append(int(ele[1]))  
            
            
            if "PsiNanopore" in file or "ELIGOS" in file:
                p_value = float(ele[2])  
                
                if p_value == 0:
                    X.append(1000)  
                else:
                    X.append(-math.log10(p_value))
            else:
                X.append(float(ele[2]))  

    y = np.array(Y)
    x = np.array(X)
    
    
    precision, recall, _ = metrics.precision_recall_curve(y, x)
    pr_auc = metrics.auc(recall, precision)
    

    plt.plot(recall, precision, color=colors[idx % len(colors)], label=f'{labels[idx]} (AUC = {pr_auc:.2f})')

plt.title('PR curve')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.legend(fontsize="small", frameon=False)
plt.savefig("./PR_m6A.pdf")

##quantification accuracy 
tools = [os.path.basename(file).split("_m6A_with_labels_converted.txt")[0] for file in files]

def process_files(dataset_files):
    # Store data for each tool and group
    tool_data = {tool: {i: {'group_values': [], 'methylation_rates': []} for i in range(7)} for tool in tools}
    
    # Create color mapping dictionary
    tool_colors = {tool: color for tool, color in zip(tools, colors)}

    for file in dataset_files:
        tool_name = file.split("/")[-1].split("_m6A_with_labels_converted.txt")[0]

        if tool_name in tool_data:
            df = pd.read_csv(file, sep="\t", header=None)
            methylation_rates = df.iloc[:, 5]
            
            if tool_name in ["xPore", "m6Anet"]:
                values = df.iloc[:, 6]
            else:
                values = df.iloc[:, 2]

            valid_data = df[methylation_rates > 0]

            quantiles = 7

            if len(valid_data) >= quantiles:
                groups = pd.cut(valid_data.iloc[:, 5], bins=quantiles, labels=False)
            else:
                print(f"Not enough data points to divide into {quantiles} groups for {tool_name}, skipping.")
                groups = []

            # Ensure groups and values are aligned
            valid_data = valid_data.copy()
            valid_data['Groups'] = groups
            valid_data['Values'] = values

            # Store both the methylation rate values and the tool's values for each group
            for group_id in range(quantiles):
                group_mask = (valid_data['Groups'] == group_id)
                tool_data[tool_name][group_id]['group_values'].extend(valid_data[group_mask]['Values'])
                tool_data[tool_name][group_id]['methylation_rates'].extend(valid_data[group_mask].iloc[:, 4])

    # Create a figure with 1:1 aspect ratio and size 40mm x 40mm
    fig, ax = plt.subplots(figsize=(40 / 25.4, 40 / 25.4))  # Convert mm to inches for plt

    group_ids = np.arange(7)
    
    # Dictionary to store correlation coefficients for each tool
    tool_correlations = {}

    # Loop over each tool and plot the medians and lines connecting them across groups
    for i, tool in enumerate(tools):
        median_values = []
        methylation_medians = []
        
        # Collect median values for this tool across all groups
        for group_id in group_ids:
            group_values = tool_data[tool][group_id]['group_values']
            meth_rates = tool_data[tool][group_id]['methylation_rates']
            
            if len(group_values) > 0 and len(meth_rates) > 0:
                median = np.median(group_values)
                median_values.append(median)
                methylation_medians.append(np.median(meth_rates))
            else:
                median_values.append(np.nan)
                methylation_medians.append(np.nan)

        # Calculate Pearson correlation between the two sets of medians
        valid_indices = ~np.isnan(median_values) & ~np.isnan(methylation_medians)
        if sum(valid_indices) >= 2:  # Need at least 2 points to calculate correlation
            corr, _ = pearsonr(np.array(median_values)[valid_indices], 
                              np.array(methylation_medians)[valid_indices])
            tool_correlations[tool] = corr
        else:
            tool_correlations[tool] = np.nan

        # Plot the line connecting the medians
        ax.plot(group_ids, median_values, color=tool_colors[tool], linewidth=0.8, marker='o')

        # Plot the medians as points
        for group_id, median in zip(group_ids, median_values):
            if not np.isnan(median):
                ax.plot(group_id, median, 'o', color=tool_colors[tool], markersize=0.3)

    # Set the x-ticks to represent groups
    ax.set_xticks(group_ids)
    ax.set_xticklabels([str(i + 1) for i in group_ids], fontsize=4)

    # Set the y-axis limits to fit the methylation rate range
    ax.set_ylim(0, 1)
    ax.set_ylabel('Methylation Rate', fontsize=4)
    ax.tick_params(axis='y', labelsize=4)

    # Save the figure
    plt.tight_layout(pad=0.1)
    plt.savefig(f"./ratio_correlation_m6A.pdf")
    plt.close()

process_files(files)
