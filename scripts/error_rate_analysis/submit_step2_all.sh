#!/bin/bash
# 批量提交step2作业脚本
# 处理results目录下各子文件夹中的mapped_1.fq文件

BASE_DIR="/home/liuycomputing/lby_FASTQ_data_202408"
WORK_DIR="/home/liuycomputing/lby_FASTQ_data_202408/new_process_20260414"
RESULT_DIR="/home/liuycomputing/lby_FASTQ_data_202408/new_process_20260414/results"

# 定义要处理的样本列表
declare -a SAMPLES=(
    "TN3055-2"
    "TN3055-9"
    "TN3132-1"
    "TN3132-2"
    "TN3132-3"
    "TN3132-4"
    "TN3196-1"
    "TN3196-2"
    "TN3196-3"
    "TN3196-4"
    "TN3196-5"
    "TN3196-6"
)

# 提交作业
for SAMPLE_NAME in "${SAMPLES[@]}"; do
    INPUT_FOLDER="${RESULT_DIR}/${SAMPLE_NAME}/"
    OUTPUT_FOLDER="${RESULT_DIR}/${SAMPLE_NAME}/"
    INPUT_FILE="mapped_1.fq"

    # 检查输入文件是否存在
    if [ ! -f "${INPUT_FOLDER}${INPUT_FILE}" ]; then
        echo "Warning: ${INPUT_FOLDER}${INPUT_FILE} does not exist, skipping..."
        continue
    fi

    # 为每个样本生成并提交一个独立的SLURM脚本
    cat > ${WORK_DIR}/job_step2_${SAMPLE_NAME}.slurm << EOF
#!/bin/bash
#SBATCH -J step2_${SAMPLE_NAME}
#SBATCH -o step2_${SAMPLE_NAME}_%j.out
#SBATCH -e step2_${SAMPLE_NAME}_%j.err
#SBATCH -p gpu_1
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

# 设置临时目录
TMP_DIR=/tmp/step2_process_${SAMPLE_NAME}
mkdir -p \${TMP_DIR}

# 复制文件到计算节点
cp ${WORK_DIR}/step2_mapped_match_multiprocess.py \${TMP_DIR}/
cp ${BASE_DIR}/refLib-SEQ2210.xlsx \${TMP_DIR}/

cd \${TMP_DIR}

echo "Processing step2: ${SAMPLE_NAME}"
echo "Input: ${INPUT_FOLDER}${INPUT_FILE}"
echo "Output: ${OUTPUT_FOLDER}"

srun python -u step2_mapped_match_multiprocess.py \\
    --library_file \${TMP_DIR}/refLib-SEQ2210.xlsx \\
    --file_folder "${INPUT_FOLDER}" \\
    --fastq_name "${INPUT_FILE}" \\
    --output_folder "${OUTPUT_FOLDER}" \\
    --mapped_read_dict_file_name "mapped_read_dict.pkl" \\
    --mapped_error_file_name "mapped_error_dict.pkl" \\
    --complete_error_file_name "complete_error_dict.pkl" \\
    --num_processes 24

# 清理临时目录
rm -rf \${TMP_DIR}
EOF

    # 提交作业
    echo "Submitting step2 job for ${SAMPLE_NAME}..."
    sbatch ${WORK_DIR}/job_step2_${SAMPLE_NAME}.slurm
    sleep 1  # 稍微等待一下，避免提交过快
done

echo "========================================="
echo "All step2 jobs submitted!"
echo "Use 'squeue -u \$USER' to check job status"
echo "========================================="
