#!/bin/bash
# 批量提交所有样本到SLURM计算节点

# 定义路径
BASE_DIR="/home/liuycomputing/lby_FASTQ_data_202408"
WORK_DIR="/home/liuycomputing/lby_FASTQ_data_202408/clustering_20260423"
RESULT_DIR="/home/liuycomputing/lby_FASTQ_data_202408/clustering_20260423/results"
LIBRARY_FILE="/home/liuycomputing/lby_FASTQ_data_202408/refLib-SEQ2210.xlsx"

# 创建结果目录
mkdir -p ${RESULT_DIR}

# 定义要处理的文件列表（格式：相对路径|样本名称）
declare -a FILES=(
    "TN3132/TN3132-1.fastq/TN3132-1_1.fastq|TN3132-1"
    "TN3132/TN3132-2.fastq/TN3132-2_1.fastq|TN3132-2"
    "TN3132/TN3132-3.fastq/TN3132-3_1.fastq|TN3132-3"
    "TN3132/TN3132-4.fastq/TN3132-4_1.fastq|TN3132-4"
    "TN3196/TN3196-1.fastq/TN3196-1_1.fastq|TN3196-1"
    "TN3196/TN3196-2.fastq/TN3196-2_1.fastq|TN3196-2"
    "TN3196/TN3196-3.fastq/TN3196-3_1.fastq|TN3196-3"
    "TN3196/TN3196-4.fastq/TN3196-4_1.fastq|TN3196-4"
    "TN3196/TN3196-5.fastq/TN3196-5_1.fastq|TN3196-5"
    "TN3196/TN3196-6.fastq/TN3196-6_1.fastq|TN3196-6"
    "TN3055/TN3055-2.fq/TN3055-2_1.fq|TN3055-2"
    "TN3055/TN3055-9.fq/TN3055-9_1.fq|TN3055-9"
)

# 为每个文件生成并提交一个独立的SLURM作业
for FILE_INFO in "${FILES[@]}"; do
    INPUT_PATH=$(echo ${FILE_INFO} | cut -d'|' -f1)
    SAMPLE_NAME=$(echo ${FILE_INFO} | cut -d'|' -f2)

    # 拼接完整路径
    FULL_INPUT_PATH="${BASE_DIR}/${INPUT_PATH}"
    INPUT_DIR=$(dirname ${FULL_INPUT_PATH})
    INPUT_FILE=$(basename ${FULL_INPUT_PATH})
    OUTPUT_FOLDER="${RESULT_DIR}/${SAMPLE_NAME}/"

    # 创建输出目录
    mkdir -p ${OUTPUT_FOLDER}

    # 生成SLURM脚本
    cat > ${WORK_DIR}/job_${SAMPLE_NAME}.slurm << EOF
#!/bin/bash
#SBATCH -J ${SAMPLE_NAME}
#SBATCH -o ${WORK_DIR}/${SAMPLE_NAME}_%j.out
#SBATCH -e ${WORK_DIR}/${SAMPLE_NAME}_%j.err
#SBATCH -p gpu_1,gpu_2,cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH -t 168:00:00

# 设置环境变量
module purge
ulimit -s unlimited
ulimit -l unlimited

# 设置OpenMP的线程数
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK

# 设置临时目录（使用计算节点的本地存储）
TMP_DIR=/tmp/fastq_process_${SAMPLE_NAME}_\${SLURM_JOB_ID}
mkdir -p \${TMP_DIR}

echo "========================================="
echo "作业信息"
echo "========================================="
echo "作业ID: \$SLURM_JOB_ID"
echo "节点名称: \$SLURM_NODELIST"
echo "CPU数: \$SLURM_CPUS_PER_TASK"
echo "样本名称: ${SAMPLE_NAME}"
echo "输入文件: ${FULL_INPUT_PATH}"
echo "输出目录: ${OUTPUT_FOLDER}"
echo "临时目录: \${TMP_DIR}"
echo "========================================="

# 复制必要文件到临时目录
cp ${WORK_DIR}/step1_findTrueIndex_multiprocess.py \${TMP_DIR}/
cp ${LIBRARY_FILE} \${TMP_DIR}/

cd \${TMP_DIR}

echo "开始处理: ${SAMPLE_NAME}"
echo "开始时间: \$(date)"

# 运行处理脚本
srun python -u step1_findTrueIndex_multiprocess.py \\
    --library_file \${TMP_DIR}/refLib-SEQ2210.xlsx \\
    --file_folder "${INPUT_DIR}/" \\
    --fastq_name "${INPUT_FILE}" \\
    --q_threshold 30.0 \\
    --output_folder "${OUTPUT_FOLDER}" \\
    --fastq_out "${SAMPLE_NAME}_classified_sequences.fq" \\
    --index_pos 110 120 \\
    --primer_pos 20 -21 \\
    --num_processes 24

echo "完成时间: \$(date)"

# 统计结果
if [ -f "${OUTPUT_FOLDER}${SAMPLE_NAME}_classified_sequences.fq" ]; then
    SEQUENCE_COUNT=\$(grep -c "^@" "${OUTPUT_FOLDER}${SAMPLE_NAME}_classified_sequences.fq")
    echo "分类成功序列数: \${SEQUENCE_COUNT}"
else
    echo "警告：输出文件不存在"
fi

# 清理临时目录
cd ${WORK_DIR}
rm -rf \${TMP_DIR}

echo "========================================="
echo "作业完成: ${SAMPLE_NAME}"
echo "========================================="
EOF

    # 提交作业
    echo "正在提交作业: ${SAMPLE_NAME}..."
    sbatch ${WORK_DIR}/job_${SAMPLE_NAME}.slurm

    if [ $? -eq 0 ]; then
        echo "✓ 作业提交成功: ${SAMPLE_NAME}"
    else
        echo "✗ 作业提交失败: ${SAMPLE_NAME}"
    fi

    sleep 2  # 等待2秒，避免提交过快
done

echo ""
echo "========================================="
echo "所有作业已提交！"
echo "========================================="
echo "查看作业状态: squeue -u \$USER"
echo "查看作业日志: ls -la ${WORK_DIR}/*.out"
echo "查看错误日志: ls -la ${WORK_DIR}/*.err"
echo "========================================="

# 显示当前提交的作业
sleep 3
echo ""
echo "当前作业队列:"
squeue -u $USER
