#!/bin/bash
# 批量处理results目录下各子文件夹中的.fq文件，测试不同序列数的自适应聚类效果

# 定义路径
BASE_DIR="/home/liuycomputing/lby_FASTQ_data_202408"
WORK_DIR="/home/liuycomputing/lby_FASTQ_data_202408/clustering_20260423"
RESULT_DIR="/home/liuycomputing/lby_FASTQ_data_202408/clustering_20260423/results"
LIBRARY_FILE="/home/liuycomputing/lby_FASTQ_data_202408/refLib-SEQ2210.xlsx"

# 检查必要文件是否存在
if [ ! -f "${LIBRARY_FILE}" ]; then
    echo "错误：library文件不存在 ${LIBRARY_FILE}"
    exit 1
fi

# 定义不同的序列数测试
declare -a SEQUENCE_NUMS=(50 100 200)

# 参数设置
MIN_CLUSTER_SIZE=5           # 最小聚类大小
NUM_PROCESSES=16             # 进程数
VERBOSE=false                # 是否输出详细信息
PRIMER_POS="20 -21"          # primer位置（从library中提取）

echo "========================================="
echo "自适应聚类算法批量处理（多序列数测试）"
echo "========================================="
echo "结果目录: ${RESULT_DIR}"
echo "Library文件: ${LIBRARY_FILE}"
echo "Primer位置: ${PRIMER_POS}（从library中提取）"
echo "测试序列数: ${SEQUENCE_NUMS[@]}"
echo "最小聚类大小: ${MIN_CLUSTER_SIZE}"
echo "进程数: ${NUM_PROCESSES}"
echo "========================================="
echo ""

# 统计样本数量
TOTAL_SAMPLES=$(find ${RESULT_DIR} -maxdepth 1 -type d | wc -l)
TOTAL_SAMPLES=$((TOTAL_SAMPLES - 1))  # 减去results目录本身

echo "发现 ${TOTAL_SAMPLES} 个样本文件夹"
echo ""

# 遍历不同的序列数设置
for MAX_SEQUENCES in "${SEQUENCE_NUMS[@]}"; do
    echo ""
    echo "========================================="
    echo "测试序列数: ${MAX_SEQUENCES}"
    echo "========================================="
    echo "测试开始时间: $(date)"
    echo ""

    # 为每个样本提交作业
    for SUBFOLDER in ${RESULT_DIR}/*/; do
        # 去掉末尾的斜杠
        SUBFOLDER=${SUBFOLDER%/}
        SAMPLE_NAME=$(basename ${SUBFOLDER})

        # 查找子文件夹中的.fq文件
        FQ_FILE=$(find ${SUBFOLDER} -maxdepth 1 -name "*.fq" | head -1)

        if [ -z "${FQ_FILE}" ]; then
            echo "✗ 跳过 ${SAMPLE_NAME}: 没有找到.fq文件"
            continue
        fi

        # 设置输出路径（包含序列数信息）
        OUTPUT_CSV="${SUBFOLDER}/${SAMPLE_NAME}_consensus_results_adaptive_seq${MAX_SEQUENCES}.csv"

        echo "========================================="
        echo "测试样本: ${SAMPLE_NAME} (序列数: ${MAX_SEQUENCES})"
        echo "========================================="
        echo "输入文件: ${FQ_FILE}"
        echo "输出文件: ${OUTPUT_CSV}"

        # 生成SLURM脚本
        cat > ${WORK_DIR}/job_adaptive_${SAMPLE_NAME}_seq${MAX_SEQUENCES}.slurm << EOF
#!/bin/bash
#SBATCH -J adaptive_${SAMPLE_NAME}_seq${MAX_SEQUENCES}
#SBATCH -o ${WORK_DIR}/adaptive_${SAMPLE_NAME}_seq${MAX_SEQUENCES}_%j.out
#SBATCH -e ${WORK_DIR}/adaptive_${SAMPLE_NAME}_seq${MAX_SEQUENCES}_%j.err
#SBATCH -p cpu,gpu_1,gpu_2,normal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=${NUM_PROCESSES}
#SBATCH -t 24:00:00

# 设置环境变量
module purge
ulimit -s unlimited
ulimit -l unlimited

echo "========================================="
echo "自适应聚类算法处理"
echo "========================================="
echo "作业ID: \$SLURM_JOB_ID"
echo "样本名称: ${SAMPLE_NAME}"
echo "输入文件: ${FQ_FILE}"
echo "输出文件: ${OUTPUT_CSV}"
echo "每Index最大序列数: ${MAX_SEQUENCES}"
echo "最小聚类大小: ${MIN_CLUSTER_SIZE}"
echo "进程数: ${NUM_PROCESSES}"
echo "Primer位置: ${PRIMER_POS}"
echo "========================================="
echo "开始时间: \$(date)"

cd ${WORK_DIR}

# 运行自适应聚类版本分析
python step2_consensus_comparison_adaptive.py \\
    --library_file ${LIBRARY_FILE} \\
    --fastq_file ${FQ_FILE} \\
    --primer_pos ${PRIMER_POS} \\
    --max_sequences_per_index ${MAX_SEQUENCES} \\
    --min_cluster_size ${MIN_CLUSTER_SIZE} \\
    --num_processes ${NUM_PROCESSES} \\
    --output_csv ${OUTPUT_CSV} \\
    $( [ "${VERBOSE}" = true ] && echo "--verbose" )

echo "完成时间: \$(date)"

# 检查输出文件
if [ -f "${OUTPUT_CSV}" ]; then
    TOTAL_LINES=\$(wc -l < "${OUTPUT_CSV}")
    SUCCESS_COUNT=\$(grep -c ",True," "${OUTPUT_CSV}")
    FAILED_COUNT=\$(grep -c ",False," "${OUTPUT_CSV}")

    # 计算成功率
    TOTAL_ANALYSES=\$((SUCCESS_COUNT + FAILED_COUNT))
    if [ \$TOTAL_ANALYSES -gt 0 ]; then
        SUCCESS_RATE=\$(echo "scale=2; \${SUCCESS_COUNT} * 100 / \${TOTAL_ANALYSES}" | bc)
    else
        SUCCESS_RATE="0.00"
    fi

    echo "========================================="
    echo "处理结果: ${SAMPLE_NAME} (序列数: ${MAX_SEQUENCES})"
    echo "========================================="
    echo "输出文件: ${OUTPUT_CSV}"
    echo "总行数: \${TOTAL_LINES} (包含表头)"
    echo "成功: \${SUCCESS_COUNT} (完全匹配)"
    echo "失败: \${FAILED_COUNT} (不完全匹配)"
    echo "成功率: \${SUCCESS_RATE}%"
    echo "========================================="

    # 统计accuracy分布
    if [ \$FAILED_COUNT -gt 0 ]; then
        echo "失败case的accuracy分布:"
        # 提取失败case的accuracy
        awk -F',' 'NR>1 && $4=="False" {print $5}' "${OUTPUT_CSV}" | sort -n | uniq -c | sort -rn | head -10 | while read count accuracy; do
            echo "  accuracy=\${accuracy}: \${count}个"
        done
    fi
else
    echo "警告：输出文件未生成 ${OUTPUT_CSV}"
fi

echo "========================================="
echo "作业完成: ${SAMPLE_NAME} (序列数: ${MAX_SEQUENCES})"
echo "========================================="
EOF

        # 提交作业
        echo "正在提交作业: adaptive_${SAMPLE_NAME}_seq${MAX_SEQUENCES}..."
        sbatch ${WORK_DIR}/job_adaptive_${SAMPLE_NAME}_seq${MAX_SEQUENCES}.slurm

        if [ $? -eq 0 ]; then
            echo "✓ 作业提交成功: ${SAMPLE_NAME} (序列数: ${MAX_SEQUENCES})"
        else
            echo "✗ 作业提交失败: ${SAMPLE_NAME} (序列数: ${MAX_SEQUENCES})"
        fi

        sleep 2  # 等待2秒，避免提交过快
    done
done

echo ""
echo "========================================="
echo "所有作业已提交！"
echo "========================================="
echo "测试设置："
echo "  序列数: ${SEQUENCE_NUMS[@]}"
echo "  样本数: ${TOTAL_SAMPLES}"
echo "  总作业数: $((${#SEQUENCE_NUMS[@]} * ${TOTAL_SAMPLES}))"
echo "========================================="
echo "查看作业状态: squeue -u \$USER"
echo "查看作业日志: ls -la ${WORK_DIR}/adaptive_*_seq*.out"
echo "查看结果文件: ls -la ${RESULT_DIR}/*/*_consensus_results_adaptive_seq*.csv"
echo "========================================="

# 显示当前提交的作业
sleep 3
echo ""
echo "当前作业队列:"
squeue -u $USER

echo ""
echo "========================================="
echo "提示：测试完成后，可以对比不同序列数的效果"
echo "========================================="
echo "1. 对比不同序列数的自适应聚类结果："
echo "   ls -la ${RESULT_DIR}/*/*_consensus_results_adaptive_seq*.csv"
echo ""
echo "2. 对比同一样本不同序列数的结果："
echo "   ls -la ${RESULT_DIR}/TN3055-2/*_consensus_results_adaptive_seq*.csv"
echo ""
echo "3. 分析成功率差异："
echo "   # 可以编写脚本比较不同序列数配置的成功率"
echo "========================================="
