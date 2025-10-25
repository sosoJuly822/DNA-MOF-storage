import argparse
import pandas as pd
import numpy as np
import difflib
import pickle
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

def build_error_dict(index_list, library):
    error_dict = {}
    lib_len = len(library[0])
    for i in index_list:
        error_dict[i] = {}
        error_dict[i]['Right'] = [0] * lib_len
        error_dict[i]['Deletion'] = [0] * lib_len
        error_dict[i]['Insertion'] = [0] * lib_len
        error_dict[i]['Substitution'] = [0] * lib_len
    
    return error_dict

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

def match_seq(ideal_seq, tmp_seq, sub_error_dict):
    matcher = difflib.SequenceMatcher(None, ideal_seq, tmp_seq)
    opcodes = matcher.get_opcodes()

    pre = 'equal'
    for tag, i1, i2, j1, j2 in opcodes:
        # print(f"{tag}: seq1[{i1}:{i2}] -> seq2[{j1}:{j2}] | {ideal_seq[i1:i2]} -> {tmp_seq[j1:j2]}")
        
        if tag == 'equal':
            if pre == 'insert':
                for i in range(i1+1,i2):
                    sub_error_dict['Right'][i] += 1
            else:
                for i in range(i1,i2):
                    sub_error_dict['Right'][i] += 1
        elif tag == 'delete':
            if i1 == i2:
                sub_error_dict['Deletion'][i1] += 1
            else:
                for i in range(i1, i2):
                    sub_error_dict['Deletion'][i] += 1
        elif tag == 'insert':
            if i1 < len(ideal_seq):
                if i1 == i2:
                    sub_error_dict['Insertion'][i1] += 1
                else:
                    for i in range(i1,i2):
                        sub_error_dict['Insertion'][i] += 1
        elif tag == 'replace':
            if i1 < len(ideal_seq):
                if i1 == i2:
                    sub_error_dict['Insertion'][i1] += 1
                else:
                    for i in range(i1,i2):
                        sub_error_dict['Substitution'][i] += 1
        
        pre = tag
    
    return sub_error_dict

def write_filtered_reads_to_fastq(output_file_path, filtered_reads):
    with open(output_file_path, 'w') as fastq_file:
        for header, sequence, quality in filtered_reads:
            fastq_record = f"{header}\n{sequence}\n+\n{quality}\n"
            fastq_file.write(fastq_record)

def save_error_dict_pickle(error_dict, file_path):
    with open(file_path, 'wb') as f:
        pickle.dump(error_dict, f)

def main(args):
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    # 创建PairwiseAligner实例
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    # 设置aligner参数
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0

    # 读取文库列表
    index_list, library = read_library(args.library_file)

    # 建库
    error_dict = build_error_dict(index_list, library)

    # 读出一条idea_seq用来获取引物序列
    for index, ideal_seq in enumerate(library):
        primer1 = Seq(ideal_seq[:20])
        primer2 = Seq(ideal_seq[-21:])
        continue

    if '_R1.' in args.fastq_name:  # 正向
        is_reverse = False
    elif '_R2.' in args.fastq_name:  # 反向互补
        is_reverse = True

    # 从FASTQ文件中读取指定质量的序列
    total_num, count_num = 0, 0
    filtered_reads = []
    for header, sequence, _, quality in read_fastq(args.file_folder + args.fastq_name):
        # 小数据量debug结果
        total_num += 1
        if total_num >= 50:
            break
        
        # 根据q_threshold筛选数据
        quality_scores = convert_quality(quality_str=quality)
        if np.mean(quality_scores) < args.q_threshold:
            continue
        
        print(sequence)

        # 根据primer1和primer2确定序列的起始和结束位置，筛选数据
        if is_reverse:
            sequence = str(Seq(sequence).reverse_complement())
    
        start_pos1, end_pos1, _, score1 = find_best_match_position(read_seq=sequence, ref_seq=primer1)
        start_pos2, end_pos2, _, score2 = find_best_match_position(read_seq=sequence, ref_seq=primer2)

        if start_pos2 - end_pos1 <= 0 or score1 < 20 or score2 < 21:  # 引物位置不正确
            continue
        
        print('The first primer position is from ', start_pos1, ' to ', end_pos1)
        print('The first primer position is from ', start_pos2, ' to ', end_pos2)

        # 根据index寻找library中与read匹配的ideal_seq
        sub_read_seq = sequence[end_pos1:end_pos1+10]
        if len(sub_read_seq) < 10:  # 索引长度不够
            continue

        best_score = -1
        most_similar_string = ''
        most_similar_index = 0

        for index, ideal_seq in enumerate(library):
            alignment = aligner.align(sub_read_seq, ideal_seq[20:30])
            score = alignment.score
            if score > best_score:
                best_score = score
                most_similar_string = ideal_seq
                most_similar_index = index

                print('index score: ', best_score)
                print('most similar index: ', most_similar_index)
            if score == 10.0:
                break
        
        if best_score == len(sub_read_seq):  # 能够找到完全匹配的索引
            # 根据payload排除index测序错误
            sub_read_seq = sequence[end_pos1+10:end_pos1+100]
            alignment = aligner.align(sub_read_seq, most_similar_string[30:120])
            best_score = alignment.score

            print('payload score: ', best_score)
            print('most similar index: ', most_similar_index)

            if best_score < 70:  # payload相似度不符合条件，找相似的 unique，获得对应的索引和序列
                seq_unique = sequence[end_pos1:start_pos2]
                for index, ideal_seq in enumerate(library):
                    alignment = aligner.align(seq_unique, ideal_seq[20:120])
                    score = alignment.score
                    if score > best_score:
                        best_score = score
                        most_similar_string = ideal_seq
                        most_similar_index = index

                        print('unique score: ', best_score)
                        print('most similar index: ', most_similar_index)

        elif best_score < len(sub_read_seq):  # 找不到完全匹配的索引，找相似的 unique，获得对应的索引和序列
            seq_unique = sequence[end_pos1:start_pos2]
            for index, ideal_seq in enumerate(library):
                alignment = aligner.align(seq_unique, ideal_seq[20:120])
                score = alignment.score
                if score > best_score:
                    best_score = score
                    most_similar_string = ideal_seq
                    most_similar_index = index

                    print('unique score: ', best_score)
                    print('most similar index: ', most_similar_index)

        # 增加一个判定条件，舍弃相似度较差的序列
        if best_score < 75:
            continue

        key = index_list[most_similar_index]
        error_dict[key] = match_seq(
            most_similar_string, 
            sequence, 
            error_dict[key]
        )
        print(error_dict[key])

        new_read = (f'@{key}_{header[1:]}', sequence, quality)
        filtered_reads.append(new_read)

        count_num += 1

        if count_num % 1000 == 0:
            print("完成对比的序列数：", count_num)
    
    save_error_dict_pickle(error_dict, file_path=args.output_folder+args.error_file_name)
    write_filtered_reads_to_fastq(output_file_path=args.output_folder+args.fastq_out, 
                                  filtered_reads=filtered_reads)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--q_threshold', type=float, default=30.0)

    parser.add_argument('--library_file', type=str,
                        default='/home/liuycomputing/lby_FASTQ_data_202408/refLib-SEQ2210.xlsx')
    
    parser.add_argument('--file_folder', type=str,
                        default='/home/liuycomputing/lby_FASTQ_data_202408/YZX-148/')
    parser.add_argument('--fastq_name', type=str, default='JZ25044006-pool-1-pool-1_combined_R1.fq')

    parser.add_argument('--output_folder', type=str, 
                        default='/home/liuycomputing/lby_FASTQ_data_202408/process_20250326/YZX-148_6/')
    parser.add_argument('--fastq_out', type=str, default='test.fq')
    parser.add_argument('--error_file_name', type=str, default='test.pkl')
    
    args = parser.parse_args()

    main(args)