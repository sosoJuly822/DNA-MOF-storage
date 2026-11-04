#!/bin/bash
# 批量提取 results 各子文件夹中 *_consensus_results_adaptive_100bp_seq*.csv 的
# index 和 consensus_data 两列，生成 consensus_100bp_seq{N}.csv 保存到对应子目录
# 格式与 recovered_100bp_seq{N}.csv 一致（表头 index,sequence）

# 定义路径
WORK_DIR="/home/liuycomputing/lby_FASTQ_data_202408/codes/clustering_20260423"
RESULT_DIR="${WORK_DIR}/results"

echo "========================================="
echo "批量提取 consensus 两列 (index, consensus_data)"
echo "========================================="
echo "结果目录: ${RESULT_DIR}"
echo "开始时间: $(date)"
echo "========================================="

# 调用 python 脚本批量处理
python3 ${WORK_DIR}/extract_consensus.py \
    --results_dir ${RESULT_DIR}

EXIT_CODE=$?

echo "========================================="
echo "结束时间: $(date)"
if [ ${EXIT_CODE} -eq 0 ]; then
    echo "处理完成！生成的文件:"
    ls -la ${RESULT_DIR}/*/consensus_100bp_seq*.csv
else
    echo "处理失败，退出码: ${EXIT_CODE}"
fi
echo "========================================="

exit ${EXIT_CODE}
