#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
批量提取 consensus 结果中的 index 和 consensus_data 两列

遍历 results 目录下各子文件夹，找到形如:
    {样本名}_consensus_results_adaptive_100bp_seq{N}.csv
的文件，提取其中的 index 和 consensus_data 两列，
以与 recovered_100bp_seq{N}.csv 相同的格式（表头 index,sequence），
保存到同一子文件夹下的 consensus_100bp_seq{N}.csv

用法:
    python3 extract_consensus_20260916.py --results_dir /path/to/results
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


def process_file(input_file: Path):
    """提取单个 csv 的 index 和 consensus_data 两列，返回 (输出文件路径, 行数)"""
    # 从文件名中解析序列数，如 adaptive_100bp_seq50.csv -> 50
    match = re.search(r'adaptive_100bp_seq(\d+)\.csv$', input_file.name)
    if match is None:
        raise ValueError(f"文件名无法解析序列数: {input_file.name}")
    reads = match.group(1)

    out_file = input_file.parent / f"consensus_100bp_seq{reads}.csv"

    df = pd.read_csv(input_file)

    if 'index' not in df.columns or 'consensus_data' not in df.columns:
        raise ValueError(f"缺少 index 或 consensus_data 列: {input_file}")

    # 提取两列，改名为 index,sequence（与 recovered_100bp_seq{N}.csv 格式一致）
    df_out = df[['index', 'consensus_data']].copy()
    df_out.columns = ['index', 'sequence']

    df_out.to_csv(out_file, index=False)
    return out_file, len(df_out)


def main():
    parser = argparse.ArgumentParser(
        description='批量提取 *_consensus_results_adaptive_100bp_seq*.csv 中的 index 和 consensus_data 两列')
    parser.add_argument(
        '--results_dir',
        default='/home/liuycomputing/lby_FASTQ_data_202408/codes/clustering_20260423/results',
        help='results 根目录（其下每个子文件夹为 一个样本）')
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    if not results_dir.is_dir():
        print(f"错误：结果目录不存在 {results_dir}")
        sys.exit(1)

    subfolders = sorted(d for d in results_dir.iterdir() if d.is_dir())
    if not subfolders:
        print(f"错误：{results_dir} 下没有子文件夹")
        sys.exit(1)

    print(f"结果目录: {results_dir}")
    print(f"发现 {len(subfolders)} 个子文件夹")
    print("=" * 50)

    total_files = 0
    total_rows = 0
    failed = 0

    for subfolder in subfolders:
        input_files = sorted(
            subfolder.glob("*_consensus_results_adaptive_100bp_seq*.csv"))
        if not input_files:
            print(f"✗ 跳过 {subfolder.name}: 未找到匹配文件")
            continue

        for input_file in input_files:
            try:
                out_file, n_rows = process_file(input_file)
                print(f"✓ {subfolder.name}: {input_file.name} -> {out_file.name} ({n_rows} 行)")
                total_files += 1
                total_rows += n_rows
            except Exception as e:
                print(f"✗ {subfolder.name}: 处理 {input_file.name} 失败: {e}")
                failed += 1

    print("=" * 50)
    print(f"处理完成: 成功 {total_files} 个文件, 共 {total_rows} 行, 失败 {failed} 个")
    if failed > 0:
        sys.exit(1)


if __name__ == '__main__':
    main()
