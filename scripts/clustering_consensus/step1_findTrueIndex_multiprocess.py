import argparse
import numpy as np
import pandas as pd
from Bio import Align
import difflib
import csv
import concurrent.futures  # Use ProcessPoolExecutor instead of ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor  # New import for multi-process
import multiprocessing
import time
import os

from Bio import Align
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def read_library(file_path):
    df = pd.read_excel(file_path)
    index_list = df.iloc[:, 0]
    library = df.iloc[:, 1]
    return index_list, library

def read_fastq(file_path):
    with open(file_path, 'r') as file:
        while True:
            header = file.readline().strip()  # 读取头部信息
            sequence = file.readline().strip()  # 读取序列信息
            plus = file.readline().strip()  # 读取加号行
            quality = file.readline().strip()  # 读取质量分数

            if not quality:
                break  # 文件结束

            yield header, sequence, plus, quality

def convert_quality(quality_str):
    return [ord(char) - 33 for char in quality_str]

def find_best_match_position(read_seq, ref_seq):
    # 执行局部比对（Smith-Waterman）
    alignments = pairwise2.align.localms(str(read_seq), str(ref_seq), 2, -1, -2, -1)

    if not alignments:
        return None, None, None

    # 获取最佳比对的结果
    best_alignment = alignments[0]

    # 提取并格式化比对信息
    alignment_str = format_alignment(*best_alignment)

    # 直接从 best_alignment 提取开始和结束位置
    start_position = best_alignment[3]  # 开始位置
    end_position = best_alignment[4]  # 结束位置
    score = best_alignment[2]  # 分数

    return start_position, end_position, alignment_str, score


# 新的 process_chunk 函数，按行范围流式处理FASTQ文件
def process_chunk_by_line(fastq_path, library, index_list, q_threshold, primer1, primer2, is_reverse, chunk_id, output_folder, start_line, end_line, index_pos):
    # 在每个子进程中独立创建aligner，避免序列化问题
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0

    temp_fastq_file = f"{output_folder}temp_chunk_{chunk_id}.fq"  # 每个进程写入临时文件

    count = 0  # 计数器
    current_line = 0  # 当前行数

    try:
        with open(fastq_path, 'r') as fastq_file, open(temp_fastq_file, 'w') as output_fastq:
            # 跳到起始行
            for _ in range(start_line):
                fastq_file.readline()
                current_line += 1

            while current_line < end_line:
                # 读取一条FASTQ记录（4行）
                header = fastq_file.readline().strip()
                if not header:
                    break  # 文件结束

                sequence = fastq_file.readline().strip()
                plus = fastq_file.readline().strip()
                quality = fastq_file.readline().strip()

                current_line += 4  # 更新行数
                count += 1

                if count % 1000 == 0:  # 每处理1000条，打印进度
                    print(f"Chunk {chunk_id}: Processed {count} sequences")

                # 根据q_threshold筛选数据
                quality_scores = convert_quality(quality_str=quality)
                if np.mean(quality_scores) < q_threshold:
                    continue  # 质量不达标，跳过

                # 根据primer1和primer2确定序列的起始和结束位置，筛选数据
                if is_reverse:
                    sequence = str(Seq(sequence).reverse_complement())

                _, end_pos1, _, score1 = find_best_match_position(read_seq=sequence, ref_seq=primer1)
                start_pos2, _, _, score2 = find_best_match_position(read_seq=sequence, ref_seq=primer2)

                if start_pos2 - end_pos1 < 0 or score1 < 40 or score2 < 42:  # 引物位置不正确，引物内有错
                    continue

                # 根据index进行分类
                most_similar_index = 0
                best_score = -1

                # 根据index寻找library中与read匹配的ideal_seq
                sub_read_seq = sequence[start_pos2-10:start_pos2]
                if len(sub_read_seq) < 10:
                    continue

                # 只进行index比对，一旦找到最佳匹配就完成分类
                for index, ideal_seq in enumerate(library):
                    alignment = aligner.align(sub_read_seq, ideal_seq[index_pos[0]:index_pos[1]])
                    score = alignment.score
                    if score > best_score:
                        best_score = score
                        most_similar_index = index
                        # print(best_score, most_similar_index)
                    if score == 10.0:  # 找到完美匹配，直接跳出
                        break
                
                # print('final match', best_score, most_similar_index)

                # 删除了payload和unique段的比对逻辑
                # 一旦index分类成功就直接完成序列的分类
                if best_score < 10:  # 设置index匹配的阈值，要求完全匹配
                    continue

                # 立即写入到临时文件
                key = index_list[most_similar_index]
                new_header = f'@{key}_{header[1:]}'
                fastq_record = f"{new_header}\n{sequence}\n+\n{quality}\n"
                output_fastq.write(fastq_record)
                output_fastq.flush()  # 强制刷新缓冲区

        print(f"Chunk {chunk_id}: 完成，共处理 {count} 条序列")
    except Exception as e:
        print(f"Chunk {chunk_id}: 发生错误 - {str(e)}")
        print(f"Chunk {chunk_id}: 错误发生在处理第 {count} 条序列时")
        import traceback
        traceback.print_exc()
        raise

def write_filtered_reads_to_fastq(output_file_path, filtered_reads):
    with open(output_file_path, 'w') as fastq_file:
        for header, sequence, quality in filtered_reads:
            fastq_record = f"{header}\n{sequence}\n+\n{quality}\n"
            fastq_file.write(fastq_record)

def merge_temp_files(output_folder, num_chunks, final_output_file):
    """合并所有临时文件到一个最终文件"""
    with open(final_output_file, 'w') as outfile:
        for chunk_id in range(num_chunks):
            temp_file = f"{output_folder}temp_chunk_{chunk_id}.fq"
            if os.path.exists(temp_file):
                with open(temp_file, 'r') as infile:
                    outfile.write(infile.read())
                os.remove(temp_file)  # 删除临时文件

def main(args):
    # 创建PairwiseAligner实例
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0

    # 读取文库列表，从索引1开始
    index_list, library = read_library(args.library_file)

    # 读出一条idea_seq用来获取引物序列
    for index, ideal_seq in enumerate(library):
        primer1 = Seq(ideal_seq[:args.primer_pos[0]])
        primer2 = Seq(ideal_seq[args.primer_pos[1]:])
        continue

    is_reverse = args.is_reverse

    # 流式处理FASTQ文件，按行数分配给不同进程
    num_processes = min(args.num_processes, multiprocessing.cpu_count())  # 限制进程数为 CPU 核数

    # 首先统计文件总行数，确定每个进程处理的起始行
    total_lines = 0
    with open(args.file_folder + args.fastq_name, 'r') as f:
        for _ in f:
            total_lines += 1

    total_sequences = total_lines // 4  # FASTQ每条序列占4行
    sequences_per_process = total_sequences // num_processes

    start_time = time.time()  # 记录开始时间

    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = []
        for i in range(num_processes):
            start_line = i * sequences_per_process * 4
            if i == num_processes - 1:  # 最后一个进程处理剩余的所有序列
                end_line = total_lines
            else:
                end_line = (i + 1) * sequences_per_process * 4

            future = executor.submit(process_chunk_by_line,
                                    args.file_folder + args.fastq_name,
                                    library, index_list, args.q_threshold,
                                    primer1, primer2, is_reverse,
                                    i, args.output_folder,
                                    start_line, end_line, args.index_pos)
            futures.append(future)

        # 等待所有任务完成，但不获取结果
        concurrent.futures.wait(futures)

    end_time = time.time()  # 记录结束时间
    print(f"Total time: {end_time - start_time} seconds")  # 输出运行时间

    # 合并所有临时文件
    print("正在合并临时文件...")
    merge_temp_files(args.output_folder, num_processes, args.output_folder+args.fastq_out)
    print(f"分类成功的序列已保存到: {args.output_folder+args.fastq_out}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--library_file', type=str,
                        default='/home/liuycomputing/lby_FASTQ_data_202408/refLib-SEQ2210.xlsx')

    parser.add_argument('--file_folder', type=str,
                        default='/home/liuycomputing/lby_FASTQ_data_202408/TN3055/TN3055-2.fq/')
    parser.add_argument('--fastq_name', type=str,
                        default='TN3055-2_1.fq')
    parser.add_argument('--is_reverse', action='store_true', default=False,
                        help='是否为反向互补序列，默认为False')

    parser.add_argument('--q_threshold', type=float, default=30.0)

    parser.add_argument('--output_folder', type=str,
                        default='/home/liuycomputing/lby_FASTQ_data_202408/clustering_20260423/results/TN3055-2/')
    parser.add_argument('--fastq_out', type=str, default='classified_sequences.fq',
                        help='分类成功序列的输出fastq文件名')

    parser.add_argument('--index_pos', type=int, nargs=2, default=[110, 120])
    parser.add_argument('--primer_pos', type=int, nargs=2, default=[20, -21])

    parser.add_argument('--num_processes', type=int, default=1)  # 线程数参数

    args = parser.parse_args()

    main(args)