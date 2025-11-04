import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import seaborn as sns
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42 

##replicate overlap (Jaccard index)
def read_bed_files(folder_path):
    bed_data = {}
    for filename in os.listdir(folder_path):
        if filename.endswith(".txt"):
            model_name = filename.split('_')[0]
            file_path = os.path.join(folder_path, filename)
            
            data=pd.read_csv(file_path, sep='\t', header=None)
            data[['chrom', 'start']] = data[0].str.split('_', expand=True)
            data[['end']] = data[['start']]
            bed_data[model_name] = data

    return bed_data 

def get_positions(bed_data, model):
    positions = set()
    if model in bed_data: 
        positions.update(set(zip(bed_data[model]["chrom"], bed_data[model]["start"], bed_data[model]["end"])))
    return positions

def calculate_jaccard_index(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union != 0 else 0

def save_jaccard_results(results, output_path):
    # Transpose the results and save them to file
    transposed_results = pd.DataFrame(results).T
    transposed_results.columns = ["1vs2", "1vs3", "2vs3"]
    
    transposed_results.to_csv(output_path, sep="\t", header=True, index_label="Model")
def plot_heatmap(jaccard_df, output_path):
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(jaccard_df, 
                     cmap="Reds",  
                     cbar_kws={'label': 'Jaccard Index'},  
                     linewidths=0.5,  
                     linecolor='white')  
    plt.xlabel("Comparison Pairs", fontsize=12)
    plt.ylabel("Models", fontsize=12)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()
    
def main():
    folder_path_HeLa3 = "replicate1"
    folder_path_HeLa2 = "replicate2"
    folder_path_HeLa = "replicate3"

    bed_data_HeLa3 = read_bed_files(folder_path_HeLa3)
    bed_data_HeLa2 = read_bed_files(folder_path_HeLa2)
    bed_data_HeLa = read_bed_files(folder_path_HeLa)

    results = {}

    all_models = set(bed_data_HeLa3.keys()).intersection(set(bed_data_HeLa2.keys()), set(bed_data_HeLa.keys()))

    for model in all_models:
        print("Processing " + model)

        positions_HeLa = get_positions(bed_data_HeLa, model)
        positions_HeLa2 = get_positions(bed_data_HeLa2, model)
        positions_HeLa3 = get_positions(bed_data_HeLa3, model)

        jaccard_1vs2 = calculate_jaccard_index(positions_HeLa, positions_HeLa2)
        jaccard_1vs3 = calculate_jaccard_index(positions_HeLa, positions_HeLa3)
        jaccard_2vs3 = calculate_jaccard_index(positions_HeLa2, positions_HeLa3)

        results[model] = [jaccard_1vs2, jaccard_1vs3, jaccard_2vs3]

    output_path = "jaccard_results.tsv"
    save_jaccard_results(results, output_path)
    output_path_heatmap = "jaccard_heatmap.png"
    plot_heatmap(jaccard_df, output_path_heatmap)
    
if __name__ == "__main__":
    main()

##signal correlation (Pearson’s r) 
def read_bed_files(folder_path):
    bed_data = {}
    
    for filename in os.listdir(folder_path):
        if filename.endswith("_m6A_with_labels_converted.txt"):
            model_name = filename.split('_')[0]
            file_path = os.path.join(folder_path, filename)
            data = pd.read_csv(file_path, sep="\t", header=0)
            data[['chrom', 'start']] = data.iloc[:, 0].astype(str).str.split('_', expand=True)
            data[['end']] = data[['start']]
            columns = data.columns.tolist()
            columns[:3] = ["chrom", "start", "end"]
            data.columns = columns
            data["site"] = list(zip(data["chrom"], data["start"], data["end"]))
            data["prediction"] = data.iloc[:, 5]  
            bed_data[model_name] = data 
    return bed_data

def correlation(bed_data_HeLa,bed_data_HeLa2,bed_data_HeLa3,model):
    try:
        df_merged = bed_data_HeLa[model].iloc[:, -2:].merge(
            bed_data_HeLa2[model].iloc[:, -2:], on='site', how='inner', suffixes=('_1', '_2')
        ).merge(
            bed_data_HeLa3[model].iloc[:, -2:], on='site', how='inner'
        )
        colnames = df_merged.columns.tolist()
        colnames = ["site", "1", "2", "3"]
        df_merged.columns = colnames
        correlation_matrix = df_merged[['1', '2', '3']].corr()   
        return correlation_matrix
    except KeyError as e:
        print(f"KeyError: {e} - Skipping model: {model}")
        return None
def plot_comprehensive_heatmap(correlation_results, output_path="comprehensive_correlation_heatmap.pdf"):
    heatmap_data = []
    tools = []
    for tool, corr_matrix in correlation_results.items():
        if corr_matrix is not None:
            tools.append(tool)
            corr_1vs2 = corr_matrix.iloc[0, 1]  # Rep1 vs Rep2
            corr_1vs3 = corr_matrix.iloc[0, 2]  # Rep1 vs Rep3
            corr_2vs3 = corr_matrix.iloc[1, 2]  # Rep2 vs Rep3
            heatmap_data.append([corr_1vs2, corr_1vs3, corr_2vs3])
    heatmap_df = pd.DataFrame(heatmap_data, 
                             index=tools, 
                             columns=['Rep1vsRep2', 'Rep1vsRep3', 'Rep2vsRep3'])
    plt.figure(figsize=(8, 6))
    sns.heatmap(heatmap_df, 
                cmap='Reds',  
                vmin=0, 
                vmax=1,
                linewidths=0.5,
                linecolor='white',
                cbar_kws={'label': 'Correlation Coefficient'})
    plt.xlabel('Replicate Comparisons', fontsize=12)
    plt.ylabel('Tools', fontsize=12)
    plt.xticks(rotation=45)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(output_path, format='pdf', bbox_inches='tight')
    return heatmap_df
    
folder_path_HeLa3 = "replicate3"
folder_path_HeLa2 = "replicate2"
folder_path_HeLa = "replicate1"
bed_data_HeLa3 = read_bed_files(folder_path_HeLa3)
bed_data_HeLa2 = read_bed_files(folder_path_HeLa2)
bed_data_HeLa = read_bed_files(folder_path_HeLa)
tools <- ["Xron", "Dorado","SingleMod","TandemMod-retrain","m6Anet",'EpiNano-retrain']
correlation_results = {}
for model in tools: 
    corr_matrix = correlation(bed_data_HeLa, bed_data_HeLa2, bed_data_HeLa3, model)
    correlation_results[model] = corr_matrix
    comprehensive_df = plot_comprehensive_heatmap(correlation_results)


##sequencing depth bias
thresholds_df = pd.read_csv("best_thresholds_precision_recall.txt", sep="\t")
glori_files = glob.glob("./*_m6A_with_labels_converted.txt")
depth_ranges = [(10, 30), (30, 50), (50, 70), (70, 100), (100, 150), (150, 200), (200, float('inf'))]
depth_labels = ["10-30", "30-50", "50-70", "70-100", "100-150", "150-200", ">=200"]
palette =["#6affb9", "#9fcc62", "#f6c365", "#2e3792", "#8e5aa2", "#ef1fff"]
tool_colors = {file.split("/")[-1].split('_m6A_with_labels_converted.txt')[0]: color for file, color in zip(glori_files, palette)}

def calculate_relative_metrics_by_depth(data_file, thresholds_df):
    label = data_file.split("/")[-1].split('_m6A_with_labels_converted.txt')[0]
    threshold_data = thresholds_df[(thresholds_df['Label'] == label) & (thresholds_df['Dataset'] == "m6A")]
    threshold_value = threshold_data['Best Threshold'].values[0]
    precision_base = threshold_data['Precision'].values[0]
    ecall_base = threshold_data['Recall'].values[0]
    data = pd.read_csv(data_file, sep="\t", header=None, usecols=[0, 1, 2, 3,4], names=["ID", "Label", "Score", "motif","Depth"])
    results = []
    counts = []
    for depth_range, depth_label in zip(depth_ranges, depth_labels):
        depth_data = data[(data["Depth"] >= depth_range[0]) & (data["Depth"] < depth_range[1])]
        
        if label in ["Differr", "ELIGOS", "xPore", "DRUMMER", "ELIGOS_diff", "Nanocompore"]:
            tp = np.sum((depth_data["Score"] <= threshold_value) & (depth_data["Label"] == 1))
            fp = np.sum((depth_data["Score"] <= threshold_value) & (depth_data["Label"] == 0))
            fn = np.sum((depth_data["Score"] > threshold_value) & (depth_data["Label"] == 1))
        else:
            tp = np.sum((depth_data["Score"] >= threshold_value) & (depth_data["Label"] == 1))
            fp = np.sum((depth_data["Score"] >= threshold_value) & (depth_data["Label"] == 0))
            fn = np.sum((depth_data["Score"] < threshold_value) & (depth_data["Label"] == 1))

        precision_depth = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall_depth = tp / (tp + fn) if (tp + fn) > 0 else 0

        relative_precision = precision_depth / precision_base if precision_base > 0 else 0
        relative_recall = recall_depth / recall_base if recall_base > 0 else 0
        results.append([label, depth_label, relative_precision, relative_recall])
        counts.append([label, depth_label, len(depth_data)])
    return (pd.DataFrame(results, columns=["Tool", "Depth", "Relative_Precision", "Relative_Recall"]),pd.DataFrame(counts, columns=["Tool", "Depth", "Count"]))

all_depths = pd.DataFrame()
all_counts = pd.DataFrame() 
for file in glori_files:
    depth_df, count_df = calculate_relative_metrics_by_depth(file, thresholds_df)
    all_depths = pd.concat([all_depths, depth_df])
    all_counts = pd.concat([all_counts, count_df])


fig_width = 25 * 0.03937  
fig_height = 30 * 0.03937
plt.figure(figsize=(fig_width, fig_height))
for tool in all_depths['Tool'].unique():
    subset = all_depths[all_depths['Tool'] == tool]
    plt.plot(np.array(subset['Depth']), np.array(subset['Relative_Recall']), marker='o', label=tool, color=tool_colors[tool], linewidth=0.8,markersize=1.5)

plt.axhline(1, color='gray', linestyle='--', linewidth=0.8)
plt.ylim(0, 2.5)
plt.xticks(rotation=90, fontsize=4)
plt.tick_params(axis='y', labelsize=4)
plt.tight_layout(pad=0.1)
plt.savefig("./m6A_relative_recall_by_coverage.pdf")
plt.close()

all_counts["Log10_Count"] = np.log10(all_counts["Count"].replace(0, np.nan))

avg_log_counts = (
    all_counts
    .groupby("Depth")["Count"]
    .median()
    .replace(0, np.nan)
    .apply(np.log10)
    .reset_index(name="Log10_Avg_Count")
)

avg_log_counts['Depth'] = pd.Categorical(
    avg_log_counts['Depth'],
    categories=depth_labels,  
    ordered=True
)

avg_log_counts = avg_log_counts.sort_values('Depth')

plt.figure(figsize=(fig_width, fig_height))
plt.bar(avg_log_counts["Depth"], avg_log_counts["Log10_Avg_Count"],
        color='gray', width=0.6)
plt.xticks(rotation=90, fontsize=4)
plt.tick_params(axis='y', labelsize=4)
plt.tight_layout(pad=0.1)
plt.savefig("./m6A_site_count_by_coverage_bar.pdf")
plt.close()

##modification level bias
methylation_ranges = [(0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.0)]
methylation_labels = ["0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"]

def calculate_relative_metrics_by_methylation(data_file, thresholds_df):
    label = data_file.split("/")[-1].split("_m6A_with_labels_converted.txt")[0]
    threshold_data = thresholds_df[(thresholds_df['Label'] == label) & (thresholds_df['Dataset'] == "m6A")]
    threshold_value = threshold_data['Best Threshold'].values[0]
    precision_base = threshold_data['Precision'].values[0]
    recall_base = threshold_data['Recall'].values[0]

    data = pd.read_csv(data_file, sep="\t", header=None, usecols=[0, 1, 2, 3, 4,5], names=["ID", "Label", "Score","Motif", "Depth", "Methylation"])
    results = []
    counts = []
    for methylation_range, methylation_label in zip(methylation_ranges, methylation_labels):
        methylation_data = data[(data["Methylation"] > methylation_range[0]) & (data["Methylation"] <= methylation_range[1])]
        
        if label in ["Differr", "ELIGOS", "xPore", "DRUMMER", "ELIGOS_diff", "Nanocompore"]:
            tp = np.sum((methylation_data["Score"] <= threshold_value) & (methylation_data["Label"] == 1))
            fp = np.sum((methylation_data["Score"] <= threshold_value) & (methylation_data["Label"] == 0))
            fn = np.sum((methylation_data["Score"] > threshold_value) & (methylation_data["Label"] == 1))
        else:
            tp = np.sum((methylation_data["Score"] >= threshold_value) & (methylation_data["Label"] == 1))
            fp = np.sum((methylation_data["Score"] >= threshold_value) & (methylation_data["Label"] == 0))
            fn = np.sum((methylation_data["Score"] < threshold_value) & (methylation_data["Label"] == 1))

        precision_methylation = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall_methylation = tp / (tp + fn) if (tp + fn) > 0 else 0

        relative_precision = precision_methylation / precision_base if precision_base > 0 else 0
        relative_recall = recall_methylation / recall_base if recall_base > 0 else 0
        results.append([label, methylation_label, relative_precision, relative_recall])
        counts.append([label, methylation_label, len(methylation_data)])
    return (pd.DataFrame(results, columns=["Tool", "Methylation", "Relative_Precision", "Relative_Recall"]),pd.DataFrame(counts, columns=["Tool", "Methylation", "Count"]))



all_methylations = pd.DataFrame()
all_counts = pd.DataFrame() 
for file in glori_files:
    methylation_df,count_df = calculate_relative_metrics_by_methylation(file, thresholds_df)
    all_methylations = pd.concat([all_methylations, methylation_df])
    all_counts = pd.concat([all_counts, count_df])

fig_width = 25 * 0.03937  
fig_height = 30 * 0.03937
plt.figure(figsize=(fig_width, fig_height))
for tool in all_methylations['Tool'].unique():
    subset = all_methylations[all_methylations['Tool'] == tool]
    plt.plot(np.array(subset['Methylation']), np.array(subset['Relative_Recall']), marker='o', label=tool, color=tool_colors[tool], linewidth=0.8,markersize=1.5)

plt.axhline(1, color='gray', linestyle='--', linewidth=0.8)
plt.ylim(0, 2.5)
plt.xticks(rotation=90, fontsize=4)
plt.tick_params(axis='y', labelsize=4)
plt.tight_layout(pad=0.1)
plt.savefig("./m6A_relative_recall_by_methylation.pdf")
plt.close()

all_counts["Log10_Count"] = np.log10(all_counts["Count"].replace(0, np.nan))

avg_log_counts = (
    all_counts
    .groupby("Methylation")["Count"]
    .mean()
    .replace(0, np.nan)
    .apply(np.log10)
    .reset_index(name="Log10_Avg_Count")
)

plt.figure(figsize=(fig_width, fig_height))
plt.bar(avg_log_counts["Methylation"], avg_log_counts["Log10_Avg_Count"],
        color='gray', width=0.6)
plt.xticks(rotation=90, fontsize=4)
plt.tick_params(axis='y', labelsize=4)
plt.tight_layout(pad=0.1)
plt.savefig("./m6A_site_count_by_Methylation_bar.pdf")
plt.close()

##Motif bias
motif_list = ["GGACC","GGACA","AGACA","GGACT","AAACT","GAACT","GAACC","AGACC","GAACA","AGACT","AAACA","AAACC"]
tool_color_dict = dict(zip(tools, colors))

def calculate_relative_metrics_by_motif(data_file, thresholds_df):
    label = data_file.split("/")[-1].split('_m6A_with_labels_converted.txt')[0]
    threshold_data = thresholds_df[(thresholds_df['Label'] == label) & (thresholds_df['Dataset'] == 'm6A')]

    data = pd.read_csv(data_file, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["ID", "Label", "Score", "Motif"])
    results = []

    for motif in motif_list:
        motif_data = data[data["Motif"] == motif]
        
        if label in ["Differr", "ELIGOS", "xPore", "DRUMMER", "ELIGOS_diff", "Nanocompore"]:
            tp = np.sum((motif_data["Score"] <= threshold_data['Best Threshold'].values  [0]) & (motif_data["Label"] == 1))
            fp = np.sum((motif_data["Score"] <= threshold_data['Best Threshold'].values  [0]) & (motif_data["Label"] == 0))
            fn = np.sum((motif_data["Score"] > threshold_data['Best Threshold'].values  [0]) & (motif_data["Label"] == 1))
        else:
            tp = np.sum((motif_data["Score"] >= threshold_data['Best Threshold'].values  [0]) & (motif_data["Label"] == 1))
            fp = np.sum((motif_data["Score"] >= threshold_data['Best Threshold'].values  [0]) & (motif_data["Label"] == 0))
            fn = np.sum((motif_data["Score"] < threshold_data['Best Threshold'].values  [0]) & (motif_data["Label"] == 1))

        precision_motif = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall_motif = tp / (tp + fn) if (tp + fn) > 0 else 0

        relative_precision = precision_motif / threshold_data['Precision'].values[0]
        relative_recall = recall_motif / threshold_data['Recall'].values[0]

        results.append([label, motif, relative_precision, relative_recall])
    
    return pd.DataFrame(results, columns=["Tool", "Motif", "Relative_Precision", "Relative_Recall"])

all_motifs = pd.DataFrame()
for file in glori_files:
    label = file.split("/")[-1].split('_m6A_with_labels_converted.txt')[0]
    if label in tools:  
        motif_df = calculate_relative_metrics(file, thresholds_df, dataset_label)
        ll_motifs = pd.concat([all_motifs, motif_df])


all_motifs['Motif'] = pd.Categorical(all_motifs['Motif'], categories=motif_list, ordered=True)
all_motifs = all_motifs.sort_values('Motif')
fig_width = 40 * 0.03937  
fig_height = 40 * 0.03937
plt.figure(figsize=(fig_width, fig_height))
for tool in all_motifs['Tool'].unique():
    subset = all_motifs[all_motifs['Tool'] == tool]
    if not subset.empty:
        x = np.array(subset['Motif'])
        y = np.array(subset['Relative_Recall'])
        plt.plot(x, y, marker='o', label=tool, color=tool_colors[tool], linewidth=0.8,markersize=1.5) 
plt.axhline(1, color='gray', linestyle='--', linewidth=0.8)
plt.ylim(0, 2.5) 
plt.xticks(rotation=90, fontsize=4)
plt.tick_params(axis='y', labelsize=4)
plt.tight_layout(pad=0.1)
plt.savefig(f"./m6A_recall_motif.pdf", bbox_inches='tight')
plt.close()





