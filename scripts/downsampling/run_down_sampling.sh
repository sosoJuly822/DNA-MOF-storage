#!/bin/bash

# TN3055 数据下采样脚本
# 下采样方式参考: /home/liuycomputing/lqt_sequencing/codes/others/run_downsampling_20251104.sh
# 对 TN3055-2_1.fq 和 TN3055-9_1.fq 按目标深度下采样，每个深度 10 次重复（不同随机种子）

# 设置工作路径
DATA_DIR="/home/liuycomputing/lby_FASTQ_data_202408/data"
OUT_DIR="${DATA_DIR}/down_sampling_data_20260915"
SAMPLES=("TN3055-2" "TN3055-9")

# 参考库大小（refLib-SEQ2210.xlsx 的 ndata 表：10800 条 x 141 bp）
SEQ_COUNT=10800
SEQ_LEN=141
REFERENCE_SIZE=$((SEQ_COUNT * SEQ_LEN))   # 1522800 bp

# 实验参数
TARGET_DEPTHS=(5 10 20 50 100 150 200)
REPLICATES=10

# 创建日志目录
mkdir -p ${OUT_DIR}/logs

# 记录实验参数
echo "实验参数:" > ${OUT_DIR}/logs/experiment_params.txt
echo "样本: ${SAMPLES[@]}" >> ${OUT_DIR}/logs/experiment_params.txt
echo "目标深度: ${TARGET_DEPTHS[@]}" >> ${OUT_DIR}/logs/experiment_params.txt
echo "重复次数: $REPLICATES" >> ${OUT_DIR}/logs/experiment_params.txt
echo "参考库: refLib-SEQ2210.xlsx (ndata 表, ${SEQ_COUNT} 条 x ${SEQ_LEN} bp)" >> ${OUT_DIR}/logs/experiment_params.txt
echo "参考库大小: ${REFERENCE_SIZE} bp" >> ${OUT_DIR}/logs/experiment_params.txt
echo "随机种子基准: 42 (同一深度系列各重复使用相同种子, 保证不同深度间嵌套)" >> ${OUT_DIR}/logs/experiment_params.txt

echo "开始下采样..."
echo "输出目录: $OUT_DIR"
echo "参考库大小: ${REFERENCE_SIZE} bp"
echo ""

# 单样本下采样函数
downsample_sample() {
    local sample=$1
    local original_file="${DATA_DIR}/TN3055/${sample}.fq/${sample}_1.fq"

    if [ ! -f "$original_file" ]; then
        echo "[${sample}] 错误: 找不到输入文件 ${original_file}"
        return 1
    fi

    # 获取原始数据信息（seqtk size 输出: reads<TAB>bases）
    local size_out=$(seqtk size $original_file)
    local original_reads=$(echo "$size_out" | cut -f1)
    local original_bases=$(echo "$size_out" | cut -f2)
    local original_depth=$(echo "scale=2; $original_bases / $REFERENCE_SIZE" | bc)

    echo "[${sample}] 原始数据信息:"
    echo "  输入文件: $original_file"
    echo "  Reads数: $original_reads"
    echo "  总碱基数: $original_bases bp"
    echo "  原始测序深度: $original_depth x"
    echo ""

    # 记录原始数据信息
    echo "输入文件: $original_file" > ${OUT_DIR}/logs/${sample}_original_info.txt
    echo "Reads数: $original_reads" >> ${OUT_DIR}/logs/${sample}_original_info.txt
    echo "总碱基数: $original_bases bp" >> ${OUT_DIR}/logs/${sample}_original_info.txt
    echo "原始测序深度: $original_depth x" >> ${OUT_DIR}/logs/${sample}_original_info.txt
    echo "参考库大小: ${REFERENCE_SIZE} bp" >> ${OUT_DIR}/logs/${sample}_original_info.txt

    # 日志 CSV 表头
    echo "目标深度,重复编号,下采样比例,实际深度,实际reads数,实际碱基数" > ${OUT_DIR}/logs/subsampling_log_${sample}.csv

    for target_depth in "${TARGET_DEPTHS[@]}"
    do
        # 计算下采样比例
        local sampling_fraction=$(echo "scale=8; $target_depth / $original_depth" | bc)

        # 跳过比例大于1的情况（超过原始深度）
        if (( $(echo "$sampling_fraction > 1" | bc -l) )); then
            echo "${target_depth},NA,超过原始深度,-,-,-" >> ${OUT_DIR}/logs/subsampling_log_${sample}.csv
            echo "[${sample}] 跳过 ${target_depth}x (超过原始深度)"
            continue
        fi

        echo "[${sample}] 处理深度 ${target_depth}x (比例: $sampling_fraction)"

        for ((rep=1; rep<=$REPLICATES; rep++))
        do
            # 为每次重复使用不同的随机种子
            local SEED=$((42 + rep))

            # 输出目录: OUT_DIR/样本/深度，文件名: 样本_深度x_重复N.fq
            local out_dir="${OUT_DIR}/${sample}/depth_${target_depth}x"
            mkdir -p "$out_dir"
            local output_file="${out_dir}/${sample}_depth${target_depth}x_rep${rep}.fq"

            # 执行下采样
            seqtk sample -s$SEED $original_file $sampling_fraction > $output_file

            # 验证实际结果
            local actual_size=$(seqtk size $output_file)
            local actual_reads=$(echo "$actual_size" | cut -f1)
            local actual_bases=$(echo "$actual_size" | cut -f2)
            local actual_depth=$(echo "scale=2; $actual_bases / $REFERENCE_SIZE" | bc)

            # 记录日志
            echo "${target_depth},${rep},${sampling_fraction},${actual_depth},${actual_reads},${actual_bases}" >> ${OUT_DIR}/logs/subsampling_log_${sample}.csv

            echo "  [${sample}] -> 重复${rep}: ${actual_depth}x (${actual_reads} reads)"
        done
        echo ""
    done

    echo "[${sample}] 下采样完成！结果保存在: ${OUT_DIR}/${sample}"
}

# 两个样本并行处理，各自独立写日志
for sample in "${SAMPLES[@]}"
do
    downsample_sample "$sample" &
done
wait

echo "全部下采样完成！"
echo "结果保存在: $OUT_DIR"
echo "日志目录: ${OUT_DIR}/logs"
