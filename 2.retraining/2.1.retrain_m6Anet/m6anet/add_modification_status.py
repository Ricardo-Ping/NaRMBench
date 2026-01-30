#!/usr/bin/env python3
import pandas as pd
import argparse

def add_modification_status(input_file, output_file, modification_status):
    # 读取数据
    print(f"Loading data from {input_file}...")
    data = pd.read_csv(input_file)

    # 添加 modification_status 列，值设为输入的值
    print("Adding modification_status column...")
    data["modification_status"] = modification_status

    # 保存结果
    print(f"Saving modified data to {output_file}...")
    data.to_csv(output_file, index=False)
    print("Process completed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add modification_status column to a CSV file.")
    parser.add_argument("input_file", help="Path to the input CSV file.")
    parser.add_argument("output_file", help="Path to the output CSV file.")
    parser.add_argument("modification_status", type=int, help="Value to assign to the modification_status column.")
    args = parser.parse_args()

    add_modification_status(args.input_file, args.output_file, args.modification_status)
