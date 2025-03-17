#!/bin/bash

# 检查参数
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <sequences.fasta> <reads_R1.fq.gz> <reads_R2.fq.gz> <output_prefix>"
    exit 1
fi

# 参数
sequences_fasta="$1"
reads_R1="$2"
reads_R2="$3"
output_prefix="$4"

# 临时文件
matched_ids="${output_prefix}_matched_ids.txt"

# 使用 seqkit grep 提取匹配的 R2 reads 并直接提取ID
zcat "$reads_R2" | seqkit grep -j 12 -R 1:8 -s -P -f "$sequences_fasta" | seqkit seq -j 12 -n -i | sed 's/ 1:N:0:.*//' > "$matched_ids"
echo "matched_ids: finished"

# 使用 seqkit grep 提取配对的 R1 和 R2 reads 并保存为 fq.gz 格式
seqkit grep -j 12 -f "$matched_ids" "$reads_R1" | pigz > "${output_prefix}_R1.fq.gz"
echo "matched_R1: finished"
seqkit grep -j 12 -f "$matched_ids" "$reads_R2" | pigz > "${output_prefix}_R2.fq.gz"
echo "matched_R2: finished"

# 删除过程文件
rm "$matched_ids"

echo "Extracted paired reads saved to ${output_prefix}_R1.fq.gz and ${output_prefix}_R2.fq.gz"


