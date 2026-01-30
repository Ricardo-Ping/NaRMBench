import os
import pandas as pd
import numpy as np
import pandas as pd

df = pd.read_csv("/disk147t/luott/project/nanopore/map/hek293t_2/features.csv",header=None,names=[
            "transcript_ID", "position", "base_type", "coverage", "ins", "ins_len",
            "del", "del_len", "fuzzy", "mis", "misA", "misC", "misG",
            "misT", "base_qual_mean", "base_qual_STD", "base_qual_count_0"])

# Read the BED file
pu = pd.read_csv("/public/xump/202409Nanopore_tools/04_sites/m5C/m5C_TAC_293T_clean.bed",sep="\t",header=None,names=["chr", "position", "position2", "ratio","strand"])

# Reset index
df.reset_index(drop=True, inplace=True)
df_filtered = df[df['base_type'] == 'C']
df_filtered.reset_index(drop=True, inplace=True)
df_filtered['position'] = df_filtered['position'].astype(int)
df_filtered['position'] = df_filtered['position'] - 1
df_filtered['chr'] = df_filtered['transcript_ID'].str.replace(r'_[FR]$', '', regex=True)

pu.set_index(['chr', 'position'], inplace=True)

df_filtered['label'] = df_filtered.apply(
    lambda row: 1 if (row['chr'], row['position']) in pu.index else 0,
    axis=1
)
match_0_df = df_filtered[df_filtered['label'] == 0]
match_0_sampled = match_0_df.sample(n=5000, random_state=42)
df_new = pd.concat([df_filtered[df_filtered['label'] == 1], match_0_sampled])

df_new = df_new.drop(columns=['chr'])
df_new.to_csv("Hek293T_train_basal.csv", index=False, header=False)
