import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
matplotlib.use('Agg')
from scipy.stats import pearsonr
import numpy as np
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42 

## Site overlap ratios heatmap
def read_bed_file(file_path):
    data=pd.read_csv(file_path, sep='\t', header=None)
    data[['chr', 'pos']] = data[0].str.split('_', expand=True)
    return data

def calculate_overlap(df1, df2):
    merged = pd.merge(df1, df2, on=['chr', 'start', 'end'], how='inner')
    union = pd.concat([df1, df2]).drop_duplicates(subset=['chr', 'start', 'end'])
    return len(merged) / len(union)
def get_prefix(file_name):
    return file_name.split('_m6A_with_labels_converted.txt')[0]

bed_dir = './' 

bed_files = [f for f in os.listdir(bed_dir) if f.endswith('_m6A_with_labels_converted.txt')]

file_prefixes = [get_prefix(f) for f in bed_files]

bed_data = {file: read_bed_file(os.path.join(bed_dir, file)) for file in bed_files}

overlap_matrix = pd.DataFrame(index=file_prefixes, columns=file_prefixes, dtype=float)

for (file1, prefix1), (file2, prefix2) in combinations(zip(bed_files, file_prefixes), 2):
    overlap_ratio = calculate_overlap(bed_data[file1], bed_data[file2])
    overlap_matrix.loc[prefix1, prefix2] = overlap_ratio
    overlap_matrix.loc[prefix2, prefix1] = overlap_ratio

for prefix in file_prefixes:
    overlap_matrix.loc[prefix, prefix] = 1.0

plt.figure(figsize=(12, 12))
sns.heatmap(overlap_matrix, annot=False, cmap='Reds', square=True, annot_kws={'size': 25})
plt.xticks(fontsize=10, rotation=90)
plt.yticks(fontsize=10, rotation=0)
plt.savefig("Site_overlap_ratios.pdf")

#Pearson’s correlation of predicted modification levels between different quantifiable tools
def create_pairwise_plot(tool_data_dict, groundtruth_data, dataset_name):
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.stats import pearsonr
    import numpy as np

    sns.set_style("white") 

    tool_names = list(tool_data_dict.keys())
    all_tools = [dataset_name] + tool_names
    n = len(all_tools)

    fig, axes = plt.subplots(n, n, figsize=(10, 10), squeeze=False)

    for i, tool1 in enumerate(all_tools):
        for j, tool2 in enumerate(all_tools):
            ax = axes[i, j]

            ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

            if i == j:
                if tool1 == dataset_name:
                    for tool in tool_names:
                        sns.kdeplot(groundtruth_data[tool][dataset_name], ax=ax, fill=True)
                else:
                    sns.kdeplot(tool_data_dict[tool1][tool1], ax=ax, fill=True)
            elif i > j:
                if tool1 == dataset_name:
                    merged = pd.merge(groundtruth_data[tool2], tool_data_dict[tool2], on="chr_pos")
                    merged = merged[merged[dataset_name] != 0]
                    x = merged[tool2]
                    y = merged[dataset_name]
                elif tool2 == dataset_name:
                    merged = pd.merge(groundtruth_data[tool1], tool_data_dict[tool1], on="chr_pos")
                    merged = merged[merged[dataset_name] != 0]
                    x = merged[dataset_name]
                    y = merged[tool1]
                else:
                    merged = pd.merge(tool_data_dict[tool1], tool_data_dict[tool2], on="chr_pos")
                    x = merged[tool2]
                    y = merged[tool1]
                sns.kdeplot(x=x, y=y, ax=ax, cmap="Blues", fill=True)
            else:
                if tool1 == dataset_name:
                    merged = pd.merge(groundtruth_data[tool2], tool_data_dict[tool2], on="chr_pos")
                    merged = merged[merged[dataset_name] != 0]
                    x = merged[dataset_name]
                    y = merged[tool2]
                elif tool2 == dataset_name:
                    merged = pd.merge(groundtruth_data[tool1], tool_data_dict[tool1], on="chr_pos")
                    merged = merged[merged[dataset_name] != 0]
                    x = merged[dataset_name]
                    y = merged[tool1]
                else:
                    merged = pd.merge(tool_data_dict[tool1], tool_data_dict[tool2], on="chr_pos")
                    x = merged[tool2]
                    y = merged[tool1]

                valid_data = merged[[tool2, tool1]].dropna()
                valid_data = valid_data[np.isfinite(valid_data[tool2]) & np.isfinite(valid_data[tool1])]

                ax.set_axis_off()
                if len(valid_data) > 1:
                    r, _ = pearsonr(valid_data[tool2], valid_data[tool1])
                    ax.text(0.5, 0.55, f"r={r:.2f}", ha="center", va="center", fontsize=16)
                    ax.text(0.5, 0.35, f"n={len(valid_data)}", ha="center", va="center", fontsize=10)
                else:
                    ax.text(0.5, 0.5, "Not enough data", ha="center", va="center", fontsize=10)

            if j == 0:
                ax.set_ylabel("ground truth" if tool1 == dataset_name else tool1, fontsize=10)
                ax.tick_params(left=True, labelleft=True)

            if i == n - 1:
                ax.set_xlabel("ground truth" if tool2 == dataset_name else tool2, fontsize=10)
                ax.tick_params(bottom=True, labelbottom=True)

    plt.subplots_adjust(wspace=0.0, hspace=0.0)
    plt.savefig(f"./methlation_cor/{dataset_name}_methylation_correlation_plot.pdf", format="pdf", bbox_inches="tight")

tool_files = {
    "Xron": "./Xron_m6A_with_labels_converted.txt",
    "Dorado": "./Dorado_m6A_with_labels_converted.txt",
    "SingleMod": "./SingleMod_m6A_with_labels_converted.txt",
    "TandemMod-retrain": "./TandemMod-retrain_m6A_with_labels_converted.txt",
    "m6Anet": "./m6Anet_m6A_with_labels_converted.txt"
}
tool_data_dict = {}
for tool_name, file_path in tool_files.items():
    tool_data = pd.read_csv(file_path, sep="\t", header=None)
    tool_data = tool_data[tool_data.apply(lambda row: len(row.dropna()) > 3, axis=1)]
    tool_data = tool_data.iloc[:, [0, 2]].copy()      
    tool_data[['chr', 'pos']] = tool_data[0].str.split('_', expand=True)
    tool_data[tool_data.columns[0]] = pd.to_numeric(tool_data[tool_data.columns[0]], errors='coerce')
    tool_data = tool_data.rename(columns={tool_data.columns[0]: f"{tool_name}"})
    tool_data['chr_pos'] = tool_data['chr'] + '_' + tool_data['pos']
    tool_data = tool_data[['chr_pos', f"{tool_name}"]]
    tool_data_dict[tool_name] = tool_data

groundtruth_data = {}
for tool_name, file_path in tool_files.items():
    tool_data = pd.read_csv(file_path, sep="\t", header=None)
    tool_data = tool_data[tool_data.apply(lambda row: len(row.dropna()) > 3, axis=1)]
    groundtruth_tool_data = tool_data.iloc[:, [0, 5]].copy() 
    
    groundtruth_tool_data[['chr', 'pos', 'strand']] = groundtruth_tool_data[0].str.split('_', expand=True)
    groundtruth_tool_data = groundtruth_tool_data.drop(columns=[0, 'strand'])
    groundtruth_tool_data[groundtruth_tool_data.columns[0]] = pd.to_numeric(groundtruth_tool_data[groundtruth_tool_data.columns[0]], errors='coerce')
    groundtruth_tool_data = groundtruth_tool_data.rename(columns={groundtruth_tool_data.columns[0]: f"groundtruth"})
    groundtruth_tool_data['chr_pos'] = groundtruth_tool_data['chr'] + '_' + groundtruth_tool_data['pos']
    groundtruth_tool_data = groundtruth_tool_data[['chr_pos', f"groundtruth"]]
    groundtruth_tool_data = groundtruth_tool_data[groundtruth_tool_data[f"groundtruth"] != 0] 
    groundtruth_data[tool_name] = groundtruth_tool_data


create_pairwise_plot(tool_data_dict, groundtruth_data, "groundtruth")
