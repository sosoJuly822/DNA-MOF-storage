#!/bin/bash
# 批量提交独立作业脚本（备选方案）
# 如果不想使用数组作业，可以使用这个脚本为每个文件提交独立的SLURM作业

BASE_DIR="/home/liuycomputing/lby_FASTQ_data_202408"
WORK_DIR="/home/liuycomputing/lby_FASTQ_data_202408/new_process_20260414"
RESULT_DIR="/home/liuycomputing/lby_FASTQ_data_202408/new_process_20260414/results"

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

# 提交作业
for FILE_INFO in "${FILES[@]}"; do
    INPUT_PATH=$(echo ${FILE_INFO} | cut -d'|' -f1)
    SAMPLE_NAME=$(echo ${FILE_INFO} | cut -d'|' -f2)

    # 为每个文件生成并提交一个独立的SLURM脚本
    cat > ${WORK_DIR}/job_${SAMPLE_NAME}.slurm << EOF
#!/bin/bash
#SBATCH -J ${SAMPLE_NAME}
#SBATCH -o ${SAMPLE_NAME}_%j.out
#SBATCH -e ${SAMPLE_NAME}_%j.err
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
TMP_DIR=/tmp/fastq_process_${SAMPLE_NAME}
mkdir -p \${TMP_DIR}

# 复制文件到计算节点
cp ${WORK_DIR}/step1_findTrueIndexAndMatch_multiprocess.py \${TMP_DIR}/
cp ${BASE_DIR}/refLib-SEQ2210.xlsx \${TMP_DIR}/

# 拼接完整路径
FULL_INPUT_PATH="${BASE_DIR}/${INPUT_PATH}"
INPUT_DIR=\$(dirname \${FULL_INPUT_PATH})
INPUT_FILE=\$(basename \${FULL_INPUT_PATH})
OUTPUT_FOLDER="${RESULT_DIR}/${SAMPLE_NAME}/"

mkdir -p \${OUTPUT_FOLDER}

cd \${TMP_DIR}

echo "Processing: ${SAMPLE_NAME}"
echo "Input: \${FULL_INPUT_PATH}"

srun python -u step1_findTrueIndexAndMatch_multiprocess.py \\
    --library_file \${TMP_DIR}/refLib-SEQ2210.xlsx \\
    --file_folder \${INPUT_DIR}/ \\
    --fastq_name "\${INPUT_FILE}" \\
    --q_threshold 30.0 \\
    --output_folder \${OUTPUT_FOLDER} \\
    --fastq_out "mapped_1.fq" \\
    --error_file_name "whole_error_dict_1.pkl" \\
    --index_length 10 \\
    --index_pos 110 120 \\
    --primer_pos 20 -21 \\
    --payload_threshold 70.0 \\
    --unique_threshold 75.0 \\
    --num_processes 24

# 清理临时目录
rm -rf \${TMP_DIR}
EOF

    # 提交作业
    echo "Submitting job for ${SAMPLE_NAME}..."
    sbatch ${WORK_DIR}/job_${SAMPLE_NAME}.slurm
    sleep 1  # 稍微等待一下，避免提交过快
done

echo "========================================="
echo "All jobs submitted!"
echo "Use 'squeue -u \$USER' to check job status"
echo "========================================="
