import argparse
import numpy as np
import pandas as pd
from Bio import Align
import difflib
import csv
import concurrent.futures  # Use ProcessPoolExecutor instead of ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor  # New import for multi-process
import multiprocessing
import pickle
import time
import os

from Bio import Align
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def read_library(file_path):  
    df = pd.read_excel(file_path)
    index_list = list(df.iloc[:, 0])
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

def build_read_dict(index_list, library):
    read_dict = {}
    lib_len = len(library[0])
    for i in index_list:
        read_dict[i] = {}
        read_dict[i]['A'] = [0] * lib_len
        read_dict[i]['T'] = [0] * lib_len
        read_dict[i]['C'] = [0] * lib_len
        read_dict[i]['G'] = [0] * lib_len
        read_dict[i]['D'] = [0] * lib_len  # 表示这个位置是空的
        read_dict[i]['I'] = [0] * lib_len  # 表示这个位置被插入
        read_dict[i]['N'] = [0] * lib_len  # 表示这个位置是低质量测序
    
    return read_dict

def build_error_dict(index_list, library):
    lib_len = len(library[0])

    error_dict = {}
    for i in index_list:
        error_dict[i] = {}
        error_dict[i]['Right'] = [0] * lib_len
        error_dict[i]['Deletion'] = [0] * lib_len
        error_dict[i]['Insertion'] = [0] * lib_len
        error_dict[i]['Substitution'] = [0] * lib_len
        error_dict[i]['LowQuality'] = [0] * lib_len

    complete_error_dict = {}
    complete_error_dict['Right'] = [0] * lib_len
    complete_error_dict['Deletion'] = [0] * lib_len
    complete_error_dict['Insertion'] = [0] * lib_len
    complete_error_dict['Substitution'] = [0] * lib_len
    complete_error_dict['LowQuality'] = [0] * lib_len

    return error_dict, complete_error_dict

def match_seq(ideal_seq, tmp_seq, sub_read_dict):
    matcher = difflib.SequenceMatcher(None, ideal_seq, tmp_seq)
    opcodes = matcher.get_opcodes()

    pre = 'equal'
    for tag, i1, i2, j1, j2 in opcodes:
        # print(f"{tag}: seq1[{i1}:{i2}] -> seq2[{j1}:{j2}] | {ideal_seq[i1:i2]} -> {tmp_seq[j1:j2]}")
        
        if tag == 'equal':  # 正确的情况，直接用ideal_seq做key就行
            if pre == 'insert':
                for i in range(i1+1,i2):
                    key = ideal_seq[i]
                    sub_read_dict[key][i] += 1
            else:
                for i in range(i1,i2):
                    key = ideal_seq[i]
                    sub_read_dict[key][i] += 1
        elif tag == 'delete':  # 如果这个位置空缺，就赋值给Deletion
            if i1 == i2:
                sub_read_dict['D'][i1] += 1
            else:
                for i in range(i1, i2):
                    sub_read_dict['D'][i] += 1
        elif tag == 'insert':
            if i1 < len(ideal_seq):
                if i1 == i2:
                    sub_read_dict['I'][i1] += 1
                else:
                    for i in range(i1,i2):
                        sub_read_dict['I'][i] += 1
        elif tag == 'replace':  # 替换的情况，要用tmp_seq的对应位置做key
            if i1 < len(ideal_seq):
                key = tmp_seq[j1]
                sub_read_dict[key][i1] += 1  # !!! 长替换只取头的位置
        
        pre = tag
    
    return sub_read_dict

# 新的 process_chunk 函数，用于处理每个数据块
def process_chunk(sequence_data_chunk, library, index_list, is_reverse, chunk_id):
    local_read_dict = build_read_dict(index_list, library)

    count = 0  # 初始化计数器
    total = len(sequence_data_chunk)  # 获取当前块的总数据量

    for sequence_data in sequence_data_chunk:
        count += 1  # 每处理一条数据，计数器加1
        if count % 1000 == 0:  # 每处理1000条，打印进度
            print(f"Chunk {chunk_id}: Processed {count}/{total} sequences")

        # 获取序列内容
        header, sequence, plus, quality = sequence_data

        # 根据primer1和primer2确定序列的起始和结束位置，筛选数据
        if is_reverse:
            sequence = str(Seq(sequence).reverse_complement())

        # 根据header区分index
        head = header.split('_')[0][1:]

        # 逐位比对，统计错误率
        index = index_list.index(head)
        local_read_dict[head] = match_seq(
            library[index], 
            sequence, 
            local_read_dict[head]
        )

    return local_read_dict

def save_dict_pickle(to_dict, file_path):
    with open(file_path, 'wb') as f:
        pickle.dump(to_dict, f)

def main(args):
    # 读取文库列表，从索引1开始
    index_list, library = read_library(args.library_file)

    if '_1.' in args.fastq_name:  # 正向
        is_reverse = False
    elif '_2.' in args.fastq_name:  # 反向互补
        is_reverse = True

    # 使用线程池并行处理FASTQ数据
    sequence_data_list = list(read_fastq(args.file_folder + args.fastq_name))

    # 根据进程数划分数据集
    num_processes = min(args.num_processes, multiprocessing.cpu_count())  # 限制进程数为 CPU 核数
    chunk_size = len(sequence_data_list) // num_processes
    sequence_chunks = [sequence_data_list[i:i + chunk_size] for i in range(0, len(sequence_data_list), chunk_size)]

    start_time = time.time()  # 记录开始时间
    global_read_dict = build_read_dict(index_list, library)

    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = [
            executor.submit(process_chunk, 
                            chunk, library, index_list, is_reverse, 
                            i)
            for i, chunk in enumerate(sequence_chunks)
        ]

        for future in concurrent.futures.as_completed(futures):
            local_read_dict = future.result()

            for key in index_list:
                for base_type in global_read_dict[key]:
                    global_read_dict[key][base_type] += np.array(local_read_dict[key][base_type])

    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    save_dict_pickle(global_read_dict, file_path=args.output_folder+args.mapped_read_dict_file_name)

    # 计算最终序列并ideal_seq比对
    bases = ['A', 'T', 'C', 'G', 'D', 'I', 'N']
    error_dict, complete_error_dict = build_error_dict(index_list, library)
    for index, ideal_seq in enumerate(library):
        key = index_list[index]
        for i, target_base in enumerate(ideal_seq):
            max_base, max_value = max(zip(bases, (global_read_dict[key][base][i] for base in bases)), key=lambda x: x[1])

            if max_value == 0:
                continue

            # print(max_base, ': ', max_value)

            if max_base == target_base:
                error_dict[key]['Right'][i] += 1
                complete_error_dict['Right'][i] += 1
            elif max_base == 'D':
                error_dict[key]['Deletion'][i] += 1
                complete_error_dict['Deletion'][i] += 1
            elif max_base == 'I':
                error_dict[key]['Insertion'][i] += 1
                complete_error_dict['Insertion'][i] += 1
            elif max_base == 'N':
                error_dict[key]['LowQuality'][i] += 1
                complete_error_dict['LowQuality'][i] += 1
            else:
                error_dict[key]['Substitution'][i] += 1
                complete_error_dict['Substitution'][i] += 1

    save_dict_pickle(error_dict, file_path=args.output_folder+args.mapped_error_file_name)
    save_dict_pickle(complete_error_dict, file_path=args.output_folder+args.complete_error_file_name)

    end_time = time.time()  # 记录结束时间
    print(f"Total time: {end_time - start_time} seconds")  # 输出运行时间


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--library_file', type=str, 
                        default='/home/liuycomputing/lby_FASTQ_data_202408/refLib-SEQ2210.xlsx')
    
    parser.add_argument('--file_folder', type=str,
                        default='/home/liuycomputing/lby_FASTQ_data_202408/process_20250326/YZX-148_6/')
    parser.add_argument('--fastq_name', type=str, default='mapped_1.fq')

    parser.add_argument('--output_folder', type=str, 
                        default='/home/liuycomputing/lby_FASTQ_data_202408/process_20250326/YZX-148_6/')
    parser.add_argument('--mapped_read_dict_file_name', type=str, default='mapped_read_dict.pkl')
    parser.add_argument('--mapped_error_file_name', type=str, default='mapped_error_dict.pkl')
    parser.add_argument('--complete_error_file_name', type=str, default='complete_error_dict.pkl')

    parser.add_argument('--num_processes', type=int, default=1)  # 线程数参数

    args = parser.parse_args()

    main(args)