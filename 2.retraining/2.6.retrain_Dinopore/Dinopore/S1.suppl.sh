#!/bin/bash

sam="IVET_mod.combined.gm2.sam"   # 输入文件
sortedBam="IVET_mod.combined.gm2.sorted.bam"   # 输出排序后的文件

# 使用 samtools 进行转换并捕捉错误输出
error_output=$(samtools view -@ 1 -S $sam -b 2>&1)

# 检查是否有报错
if echo "$error_output" | grep -q "Parse error"; then
    echo "Error occurred during samtools view/sort, extracting error line number"

    # 提取报错行号
    error_line=$(echo "$error_output" | grep -oP 'line \K[0-9]+')

    # 如果提取到错误行号
    if [ ! -z "$error_line" ]; then
        echo "Deleting error line $error_line from the SAM file"
        
        # 使用 sed 删除报错行
        sed -i "${error_line}d" $sam

        # 重新执行 samtools 处理文件
        echo "Retrying samtools view and sort after deleting error line"
        samtools view -@ 1 -S $sam -b | samtools sort -o $sortedBam - ; samtools index $sortedBam
    else
        echo "No specific line number found in the error output."
    fi
else
    # 如果没有报错，则继续正常的处理
    samtools view -@ 1 -S $sam -b | samtools sort -o $sortedBam - ; samtools index $sortedBam
fi
