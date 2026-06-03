import argparse
import pandas as pd
import numpy as np
import difflib
import pickle
import time

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
            print(f"{tag}: seq1[{i1}:{i2}] -> seq2[{j1}:{j2}] | {ideal_seq[i1:i2]} -> {tmp_seq[j1:j2]}")
            if i1 < len(ideal_seq):
                key = tmp_seq[j1]
                sub_read_dict[key][i1] += 1  # !!! 长替换只取头的位置
                    
        
        pre = tag
    
    return sub_read_dict

def save_dict_pickle(to_dict, file_path):
    with open(file_path, 'wb') as f:
        pickle.dump(to_dict, f)

def main(args):
    # 读取文库列表，从索引1开始
    index_list, library = read_library(args.library_file)

    # 建库
    read_dict = build_read_dict(index_list, library)

    if '_1.' in args.fastq_name:  # 正向
        is_reverse = False
    elif '_2.' in args.fastq_name:  # 反向互补
        is_reverse = True

    # 从FASTQ文件中读取指定质量的序列
    total_num = 0
    for header, sequence, _, _ in read_fastq(args.file_folder + args.fastq_name):
        # 小数据量debug结果
        total_num += 1
        if total_num >= 50:
            break

        # 根据primer1和primer2确定序列的起始和结束位置，筛选数据
        if is_reverse:
            sequence = str(Seq(sequence).reverse_complement())     

        # 根据header区分index
        head = header.split('_')[0][1:]

        # 逐位比对，统计错误率
        index = index_list.index(head)
        read_dict[head] = match_seq(
            library[index], 
            sequence, 
            read_dict[head]
        )
    
    # 计算最终序列并ideal_seq比对
    bases = ['A', 'T', 'C', 'G', 'D', 'I']
    error_dict, complete_error_dict = build_error_dict(index_list, library)
    for index, ideal_seq in enumerate(library):
        key = index_list[index]
        for i, target_base in enumerate(ideal_seq):
            max_base, max_value = max(zip(bases, (read_dict[key][base][i] for base in bases)), key=lambda x: x[1])

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
    
    save_dict_pickle(read_dict, file_path=args.output_folder+args.mapped_read_dict_file_name)
    save_dict_pickle(error_dict, file_path=args.output_folder+args.mapped_error_file_name)
    save_dict_pickle(complete_error_dict, file_path=args.output_folder+args.complete_error_file_name)


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

    args = parser.parse_args()

    main(args)