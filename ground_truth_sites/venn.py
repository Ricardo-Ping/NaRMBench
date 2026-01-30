import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# 1. 定义文件路径
file1_path = "/data/pingyc/projects/20250601_RNA/data/lu_m6a_result.transcript.update/ctr-1.totalm6A_5mer_trans.FDR.tsv"
file2_path = "/data/pingyc/projects/20250601_RNA/softwares/NaRMBench/ground_truth_sites/m6A_HEK293T.tsv"

def get_venn_data():
    # 2. 读取第一个文件 (注意：根据数据特征使用 sep='\s+' 兼容空格)
    print("正在读取文件 1...")
    df1 = pd.read_csv(file1_path, sep='\s+', engine='python')
    
    # 【新增逻辑】去除 transcript_id 为 -1 的行
    # 这里将 -1 视为字符串 "-1" 或数字 -1 统一处理
    initial_count = len(df1)
    df1 = df1[df1['transcript_id'].astype(str) != "-1"]
    removed_count = initial_count - len(df1)
    print(f"文件 1 已过滤 transcript_id 为 -1 的行（共移除 {removed_count} 条记录）")

    # 3. 读取第二个文件 (Ground Truth)
    print("正在读取文件 2...")
    df2 = pd.read_csv(file2_path, sep='\t')

    # 4. 提取位点唯一标识符 (染色体:位置:链)
    # 统一转换为字符串以防格式不匹配
    set1 = set(
        df1['Chr'].astype(str) + ":" + 
        df1['Sites'].astype(str) + ":" + 
        df1['Strand'].astype(str)
    )
    
    set2 = set(
        df2['chr'].astype(str) + ":" + 
        df2['pos'].astype(str) + ":" + 
        df2['strand'].astype(str)
    )

    # 5. 计算具体数值
    total1 = len(set1)
    total2 = len(set2)
    overlap = len(set1.intersection(set2))
    only1 = total1 - overlap
    only2 = total2 - overlap

    print("-" * 30)
    print(f"统计结果（已过滤 -1）：")
    print(f"文件 1 有效位点数: {total1}")
    print(f"文件 2 (Ground Truth) 总位点数: {total2}")
    print(f"重合 (Intersection) 位点数: {overlap}")
    print(f"文件 1 特有: {only1}")
    print(f"文件 2 特有: {only2}")
    print("-" * 30)

    # 6. 绘制韦恩图
    plt.figure(figsize=(10, 7))
    venn = venn2(
        [set1, set2], 
        set_labels=('Filtered Result (File 1)', 'Ground Truth (File 2)')
    )
    
    plt.title("m6A Sites Overlap Comparison", fontsize=14)
    
    # 保存图片
    output_png = "m6a_sites_venn.png"
    plt.savefig(output_png, dpi=300)
    print(f"韦恩图已保存至: {output_png}")

if __name__ == "__main__":
    try:
        get_venn_data()
    except Exception as e:
        print(f"发生错误: {e}")