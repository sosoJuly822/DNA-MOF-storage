import argparse
import pandas as pd
import numpy as np
from collections import Counter, defaultdict
from Bio import Align
from Bio import pairwise2
from Bio.Seq import Seq
import csv
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
import time

def read_library(file_path):
    """读取参考序列库"""
    df = pd.read_excel(file_path)
    index_list = df.iloc[:, 0]
    library = df.iloc[:, 1]
    return index_list, library

def read_fastq(file_path):
    """读取FASTQ文件"""
    with open(file_path, 'r') as file:
        while True:
            header = file.readline().strip()  # 读取头部信息
            sequence = file.readline().strip()  # 读取序列信息
            plus = file.readline().strip()  # 读取加号行
            quality = file.readline().strip()  # 读取质量分数

            if not quality:
                break  # 文件结束

            yield header, sequence, plus, quality

def extract_data_segment(sequence, primer1, primer2):
    """从序列中提取data段（去掉primer部分）"""
    def find_best_match_position(read_seq, ref_seq):
        """执行局部比对（Smith-Waterman）"""
        alignments = pairwise2.align.localms(str(read_seq), str(ref_seq), 2, -1, -2, -1)

        if not alignments:
            return None, None, None, 0

        best_alignment = alignments[0]
        start_position = best_alignment[3]
        end_position = best_alignment[4]
        score = best_alignment[2]

        return start_position, end_position, best_alignment, score

    def find_primer_positions(seq, p1, p2, min_match_score=15):
        """在序列中找到primer的位置"""
        p1_start, p1_end, _, p1_score = find_best_match_position(seq, p1)
        if p1_start is None or p1_score < min_match_score:
            return None

        p2_start, p2_end, _, p2_score = find_best_match_position(seq, p2)
        if p2_start is None or p2_score < min_match_score:
            return None

        if p1_end >= p2_start:
            return None

        return p1_start, p1_end, p2_start, p2_end

    primer_positions = find_primer_positions(sequence, primer1, primer2)

    if primer_positions is None:
        return None

    p1_start, p1_end, p2_start, p2_end = primer_positions

    # 提取data段（primer1结束到primer2开始之间）
    data_segment = sequence[p1_end:p2_start]

    return data_segment

def adaptive_clustering(sequences, min_cluster_size=5):
    """
    自适应聚类算法：根据序列分布智能选择最佳策略

    策略选择逻辑：
    1. 如果只有一个序列，直接返回
    2. 如果所有序列都相同，直接返回
    3. 分析序列分布特征
    4. 根据特征选择合适的聚类策略
    """
    if not sequences:
        return None

    if len(sequences) == 1:
        return sequences[0]

    # 统计序列分布
    sequence_counts = Counter(sequences)
    sorted_sequences = sequence_counts.most_common()

    unique_sequences = len(sequence_counts)

    # 情况1：所有序列都相同
    if unique_sequences == 1:
        return sequences[0]

    # 情况2：序列数量很少，直接多数投票
    if len(sequences) <= min_cluster_size:
        return simple_majority_vote(sequences)

    # 分析序列分布特征
    main_sequence, main_count = sorted_sequences[0]
    main_ratio = main_count / len(sequences)

    # 计算序列间相似度分布
    def calculate_similarity(seq1, seq2):
        min_len = min(len(seq1), len(seq2))
        matches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a == b)
        return matches / min_len if min_len > 0 else 0

    # 分析所有序列与主序列的相似度
    similarities = []
    for seq, count in sorted_sequences:
        sim = calculate_similarity(main_sequence, seq)
        similarities.extend([sim] * count)

    similarities.sort()

    # 根据分布特征选择策略
    median_similarity = np.median(similarities)
    mean_similarity = np.mean(similarities)
    high_similarity_count = sum(1 for s in similarities if s >= 0.95)

    print(f"  聚类分析: {len(sequences)}条序列, {unique_sequences}个唯一序列")
    print(f"  主序列占比: {main_ratio:.1%}, 中位相似度: {median_similarity:.2%}")
    print(f"  高相似度序列比例: {high_similarity_count/len(similarities):.1%}")

    # 策略选择逻辑
    if main_ratio >= 0.7 and median_similarity >= 0.95:
        # 情况1：主序列占主导且序列高度相似
        # 使用主聚类
        cluster_members = [main_sequence] * main_count

        # 添加高度相似的序列
        for seq, count in sorted_sequences[1:]:
            sim = calculate_similarity(main_sequence, seq)
            if sim >= 0.95:
                cluster_members.extend([seq] * count)

        if len(cluster_members) >= min_cluster_size:
            print(f"  策略: 主聚类（高度相似）")
            return majority_vote_consensus(cluster_members)
        else:
            print(f"  策略: 回退到简单多数投票")
            return simple_majority_vote(sequences)

    elif main_ratio >= 0.5 and mean_similarity >= 0.90:
        # 情况2：主序列占多数，相似度较高
        # 使用宽松阈值聚类
        cluster_members = []
        for seq, count in sorted_sequences:
            sim = calculate_similarity(main_sequence, seq)
            if sim >= 0.85:  # 使用更宽松的阈值
                cluster_members.extend([seq] * count)

        if len(cluster_members) >= min_cluster_size * 2:  # 要求更严格的大小要求
            print(f"  策略: 宽松阈值聚类（相似度≥85%）")
            return majority_vote_consensus(cluster_members)
        else:
            print(f"  策略: 回退到简单多数投票")
            return simple_majority_vote(sequences)

    elif unique_sequences <= 5:
        # 情况3：唯一序列不多，可能是真实变异
        # 使用加权多数投票（基于频率）
        weighted_sequences = []
        for seq, count in sorted_sequences:
            weighted_sequences.extend([seq] * count)

        print(f"  策略: 加权多数投票（真实变异）")
        return weighted_majority_vote(weighted_sequences)

    else:
        # 情况4：序列多样性高，可能混合了不同来源
        # 尝试找到高质量聚类
        # 策略：基于频率和相似度的混合聚类
        best_cluster = []
        best_cluster_quality = 0

        # 尝试不同的前几个序列作为聚类中心
        for candidate_seq, candidate_count in sorted_sequences[:5]:
            if candidate_count < min_cluster_size:
                continue

            # 构建这个候选的聚类
            temp_cluster = [candidate_seq] * candidate_count
            for seq, count in sorted_sequences:
                if seq == candidate_seq:
                    continue
                sim = calculate_similarity(candidate_seq, seq)
                if sim >= 0.90:
                    temp_cluster.extend([seq] * count)

            # 评估聚类质量
            if len(temp_cluster) >= min_cluster_size * 2:
                # 计算聚类内的平均相似度
                cluster_similarities = []
                for i, seq1 in enumerate(temp_cluster[:min(10, len(temp_cluster))]):
                    for j, seq2 in enumerate(temp_cluster[i+1:min(10, len(temp_cluster))]):
                        sim = calculate_similarity(seq1, seq2)
                        cluster_similarities.append(sim)

                cluster_quality = np.mean(cluster_similarities) if cluster_similarities else 0

                if cluster_quality > best_cluster_quality and cluster_quality > 0.90:
                    best_cluster = temp_cluster
                    best_cluster_quality = cluster_quality

        if best_cluster and len(best_cluster) >= min_cluster_size:
            print(f"  策略: 混合聚类（高质量聚类，相似度={best_cluster_quality:.2f}）")
            return majority_vote_consensus(best_cluster)
        else:
            print(f"  策略: 回退到简单多数投票")
            return simple_majority_vote(sequences)

def simple_majority_vote(sequences):
    """简单多数投票算法"""
    if not sequences:
        return ""

    seq_length = len(sequences[0])
    for seq in sequences:
        if len(seq) != seq_length:
            pass

    consensus = []

    for i in range(seq_length):
        bases = []
        for seq in sequences:
            if i < len(seq):
                bases.append(seq[i])

        base_counts = Counter(bases)

        if base_counts:
            most_common = base_counts.most_common(1)[0][0]
            consensus.append(most_common)
        else:
            consensus.append('N')

    return ''.join(consensus)

def weighted_majority_vote(sequences):
    """加权多数投票算法（考虑序列频率）"""
    return simple_majority_vote(sequences)

def majority_vote_consensus(sequences):
    """多数投票法生成共识序列"""
    if not sequences:
        return ""

    seq_length = len(sequences[0])
    for seq in sequences:
        if len(seq) != seq_length:
            pass  # 使用第一个序列的长度，其他序列如果不同则截断或填充

    consensus = []

    for i in range(seq_length):
        bases = []
        for seq in sequences:
            if i < len(seq):
                bases.append(seq[i])

        # 统计每个碱基的出现次数
        base_counts = Counter(bases)

        # 选择出现次数最多的碱基
        if base_counts:
            most_common = base_counts.most_common(1)[0][0]
            consensus.append(most_common)
        else:
            consensus.append('N')

    return ''.join(consensus)

def process_chunk(fastq_path, search_keys, key_to_index_dict, primer1, primer2, reference_dict, max_sequences_per_index, chunk_id, min_cluster_size):
    """
    处理一个chunk的key，使用自适应聚类算法
    """
    print(f"Chunk {chunk_id}: Starting with {len(search_keys)} keys")

    # 建立序列收集字典（使用index作为键）
    sequences_dict = {}
    local_search_keys = list(search_keys)

    # print('payload of ref_sequence:', reference_dict[0][20:-21])

    # 读取FASTQ文件并收集序列
    for header, sequence, plus, quality in read_fastq(fastq_path):
        if header.startswith('@'):
            try:
                # 从header中提取key，格式：@key_rest_of_header
                header_parts = header[1:].split('_', 1)
                key = header_parts[0]  # 第一部分是key

                if key not in local_search_keys:
                    continue

                # 将key映射到index
                index = key_to_index_dict[key]

                # if index == 0:
                #     print(sequence)

                if index not in sequences_dict:
                    sequences_dict[index] = []

                sequences_dict[index].append(sequence)

                if len(sequences_dict[index]) >= max_sequences_per_index:
                    local_search_keys.remove(key)
                    print(f"Chunk {chunk_id}: Key {key} (Index {index}) completed with {len(sequences_dict[index])} sequences")

            except (IndexError, ValueError, KeyError):
                continue

        if not local_search_keys:
            print(f"Chunk {chunk_id}: All target keys completed, stopping early...")
            break

    print(f"Chunk {chunk_id}: Collection finished, collected {len(sequences_dict)} indexes")

    # 分析收集到的序列（使用自适应聚类）
    results = []
    for index, sequences in sequences_dict.items():
        if index not in reference_dict:
            continue

        reference_seq = reference_dict[index]

        result = {
            'index': index,
            'num_sequences': len(sequences),
            'success': False,
            'error': None
        }

        try:
            # 提取所有序列的data段
            data_segments = []
            for seq in sequences:
                data_segment = extract_data_segment(seq, primer1, primer2)
                if data_segment:
                    data_segments.append(data_segment)

            if not data_segments:
                result['error'] = 'No valid data segments extracted'
                results.append(result)
                continue

            # 使用自适应聚类构建共识序列
            consensus_data = adaptive_clustering(data_segments, min_cluster_size=min_cluster_size)

            if not consensus_data:
                result['error'] = 'Failed to build consensus from adaptive clustering'
                results.append(result)
                continue

            # 提取参考序列的data段（固定90bp: [20:110]）
            reference_data = reference_seq[20:110]

            # 固定取90bp进行比较
            target_length = 90
            consensus_data_90 = consensus_data[:target_length]
            reference_data_90 = reference_data[:target_length]

            # 计算正确率
            correct_bases = sum(1 for a, b in zip(consensus_data_90, reference_data_90) if a == b)
            accuracy = correct_bases / target_length

            # 只有完全匹配（90/90）才算成功
            is_perfect_match = (correct_bases == target_length)

            result.update({
                'consensus_data': consensus_data_90+reference_seq[110:120],  # 包含data段和后续部分
                'reference_data': reference_data_90+reference_seq[110:120],
                'accuracy': accuracy,
                'correct_bases': correct_bases,
                'total_bases': target_length,
                'success': is_perfect_match
            })

            if not is_perfect_match:
                result['error'] = f'Not perfect match: {correct_bases}/{target_length} correct ({accuracy:.1%})'

            results.append(result)

        except Exception as e:
            result['error'] = str(e)
            results.append(result)

    successful_count = sum(1 for r in results if r['success'])
    print(f"Chunk {chunk_id}: Completed {successful_count}/{len(results)} successful analyses")

    return results

def main(args):
    print("=" * 60)
    print("Adaptive Clustering Consensus Sequence Comparison Analysis")
    print("=" * 60)

    # 读取参考序列库
    print(f"\nReading reference library: {args.library_file}")
    index_list, library = read_library(args.library_file)
    print(f"Total reference sequences: {len(library)}")

    # 建立参考序列字典
    reference_dict = {idx: str(seq) for idx, seq in enumerate(library)}

    # 建立key到index的映射（用于从FASTQ header中提取key）
    key_to_index_dict = {key: idx for idx, key in enumerate(index_list)}

    # 从library中提取primer
    print(f"Extracting primers from library (positions: {args.primer_pos[0]}, {args.primer_pos[1]})")
    for index, ideal_seq in enumerate(library):
        primer1 = Seq(ideal_seq[:args.primer_pos[0]])
        primer2 = Seq(ideal_seq[args.primer_pos[1]:])
        continue  # 只读取第一条
    print(f"Primer1: {primer1}")
    print(f"Primer2: {primer2}")

    # 建立待收集的key列表
    search_keys = list(index_list)
    print(f"\nTarget keys to process: {len(search_keys)}")
    print(f"Max sequences per index: {args.max_sequences_per_index}")
    print(f"Min cluster size: {args.min_cluster_size}")

    # 确定进程数
    num_processes = min(args.num_processes, multiprocessing.cpu_count())
    print(f"Using {num_processes} processes")

    # 将key列表分成多个chunk
    chunk_size = len(search_keys) // num_processes
    chunks = []
    for i in range(num_processes):
        start_idx = i * chunk_size
        if i == num_processes - 1:
            end_idx = len(search_keys)
        else:
            end_idx = (i + 1) * chunk_size

        chunk_keys = search_keys[start_idx:end_idx]
        chunks.append(chunk_keys)
        print(f"Chunk {i}: {len(chunk_keys)} keys (keys {start_idx}-{end_idx-1})")

    start_time = time.time()

    # 使用多进程处理
    print(f"\nStarting multi-process processing with adaptive clustering...")

    all_results = []
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = []
        for i, chunk_keys in enumerate(chunks):
            future = executor.submit(
                process_chunk,
                args.fastq_file,
                chunk_keys,
                key_to_index_dict,
                primer1,
                primer2,
                reference_dict,
                args.max_sequences_per_index,
                i,
                args.min_cluster_size
            )
            futures.append(future)

        # 收集所有结果
        for future in concurrent.futures.as_completed(futures):
            chunk_results = future.result()
            all_results.extend(chunk_results)

    end_time = time.time()

    # 输出结果
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total indexes analyzed: {len(all_results)}")
    print(f"Successful analyses: {sum(1 for r in all_results if r['success'])}")
    print(f"Failed analyses: {sum(1 for r in all_results if not r['success'])}")
    print(f"Total processing time: {end_time - start_time:.2f} seconds")

    # 详细输出
    if args.verbose:
        print("\n" * 60)
        print("DETAILED RESULTS (First 20)")
        print("=" * 60)

        for result in all_results[:20]:
            if result['success']:
                print(f"\nIndex {result['index']}:")
                print(f"  ✓ SUCCESS - Perfect match!")
                print(f"  Sequences used: {result['num_sequences']}")
                print(f"  Accuracy: {result['accuracy']:.1%} ({result['correct_bases']}/{result['total_bases']})")
                print(f"  Consensus data: {result['consensus_data']}")
                print(f"  Reference data: {result['reference_data']}")
            else:
                print(f"\nIndex {result['index']}:")
                print(f"  ✗ FAILED - {result['error']}")
                print(f"  Sequences used: {result['num_sequences']}")
                if 'accuracy' in result:
                    print(f"  Accuracy: {result['accuracy']:.1%} ({result['correct_bases']}/{result['total_bases']})")
                    print(f"  Consensus data: {result['consensus_data']}")
                    print(f"  Reference data: {result['reference_data']}")

        # 统计聚合结果分布
        strategy_counts = {}
        error_details = {}

        for result in all_results:
            if not result['success']:
                error = result.get('error', 'unknown')
                if '策略' in error:
                    strategy = error.split(':')[0].strip()
                    strategy_counts[strategy] = strategy_counts.get(strategy, 0) + 1
                error_details[error] = error_details.get(error, 0) + 1

        if strategy_counts:
            print("\n" * 60)
            print("STRATEGY USAGE DISTRIBUTION (Failed Cases)")
            print("=" * 60)
            for strategy, count in sorted(strategy_counts.items()):
                print(f"  {strategy}: {count} cases")

    # 保存结果到CSV
    if args.output_csv:
        print(f"\nSaving results to: {args.output_csv}")

        with open(args.output_csv, 'w', newline='') as csvfile:
            fieldnames = ['index', 'num_sequences', 'success', 'accuracy',
                         'correct_bases', 'total_bases', 'consensus_data',
                         'reference_data', 'error']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            for result in all_results:
                row = {
                    'index': result['index'],
                    'num_sequences': result['num_sequences'],
                    'success': result['success'],
                    'consensus_data': result.get('consensus_data', ''),
                    'reference_data': result.get('reference_data', ''),
                    'error': result.get('error', '')
                }

                if 'accuracy' in result:
                    row.update({
                        'accuracy': f"{result['accuracy']:.4f}",
                        'correct_bases': result['correct_bases'],
                        'total_bases': result['total_bases']
                    })
                else:
                    row.update({
                        'accuracy': '',
                        'correct_bases': '',
                        'total_bases': ''
                    })

                writer.writerow(row)

    # 统计accuracy分布
    if len(all_results) > 0:
        accuracies = [r['accuracy'] for r in all_results if 'accuracy' in r]
        print(f"\n" + "=" * 60)
        print("ACCURACY STATISTICS")
        print("=" * 60)
        print(f"Mean accuracy: {np.mean(accuracies):.2%}")
        print(f"Median accuracy: {np.median(accuracies):.2%}")
        print(f"Std accuracy: {np.std(accuracies):.2%}")
        print(f"Min accuracy: {np.min(accuracies):.2%}")
        print(f"Max accuracy: {np.max(accuracies):.2%}")

        # Accuracy分布
        perfect_match = sum(1 for acc in accuracies if acc == 1.0)
        high_accuracy = sum(1 for acc in accuracies if 0.95 <= acc < 1.0)
        medium_accuracy = sum(1 for acc in accuracies if 0.80 <= acc < 0.95)
        low_accuracy = sum(1 for acc in accuracies if acc < 0.80)

        total_analyses = len(all_results)
        print(f"\nAccuracy distribution:")
        print(f"  100% (Perfect match): {perfect_match} ({perfect_match/total_analyses:.1%})")
        print(f"  95-100%: {high_accuracy} ({high_accuracy/total_analyses:.1%})")
        print(f"   80-95%: {medium_accuracy} ({medium_accuracy/total_analyses:.1%})")
        print(f"  < 80%: {low_accuracy} ({low_accuracy/total_analyses:.1%})")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Adaptive clustering consensus sequence comparison analysis')

    parser.add_argument('--library_file', type=str,
                        default='/home/liuycomputing/lby_FASTQ_data_202408/refLib-SEQ2210.xlsx',
                        help='Path to reference library Excel file')
    parser.add_argument('--fastq_file', type=str, 
                        default='/home/liuycomputing/lby_FASTQ_data_202408/clustering_20260423/results/TN3132-1/TN3132-1_classified_sequences.fq',
                        help='Path to classified FASTQ file')
    parser.add_argument('--primer_pos', type=int, nargs=2, default=[20, -21],
                        help='Positions to extract primers from reference sequences (default: [20, -21])')
    parser.add_argument('--max_sequences_per_index', type=int, default=50,
                        help='Maximum number of sequences to collect per index (default: 50)')
    parser.add_argument('--min_cluster_size', type=int, default=5,
                        help='Minimum cluster size for consensus building (default: 5)')
    parser.add_argument('--num_processes', type=int, default=1,
                        help='Number of processes to use (default: 4)')
    parser.add_argument('--output_csv', type=str, default=None,
                        help='Path to output CSV file')
    parser.add_argument('--verbose', action='store_true',
                        help='Print detailed results for each index')

    args = parser.parse_args()

    main(args)