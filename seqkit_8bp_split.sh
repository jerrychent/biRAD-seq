# usage：bash seqkit_8bp_split_without_index.sh <sequences.fasta> <reads_R1.fq.gz> <reads_R2.fq.gz> <output_prefix>
#!/bin/bash

# 参数检查
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <sequences.fasta> <reads_R1.fq.gz> <reads_R2.fq.gz> <output_prefix>"
    exit 1
fi

sequences_fasta="$1"
reads_R1="$2"
reads_R2="$3"
output_prefix="$4"

# 工具检查
if ! command -v seqkit &> /dev/null; then
    echo "Error: seqkit not found. Please install seqkit and try again."
    exit 2
fi

if ! command -v pigz &> /dev/null && ! command -v gzip &> /dev/null; then
    echo "Error: Neither pigz nor gzip found. Please install one of them."
    exit 2
fi

# 使用 pigz 或 gzip
compressor="pigz"
if ! command -v pigz &> /dev/null; then
    compressor="gzip"
fi

# 输入检查
if [ ! -s "$reads_R1" ] || [ ! -s "$reads_R2" ] || [ ! -s "$sequences_fasta" ]; then
    echo "Error: One or more input files are missing or empty."
    exit 3
fi

# 临时文件设置
tmp_dir=$(mktemp -d)
trap "rm -rf $tmp_dir" EXIT
matched_ids="$tmp_dir/matched_ids.txt"

# 提取 R1 匹配 Reads
zcat "$reads_R1" | seqkit grep -j 12 -R 1:8 -P -s -f "$sequences_fasta" | $compressor > "${output_prefix}_R1.fq.gz"
echo "Matched R1 reads saved to ${output_prefix}_R1.fq.gz."

# 提取 R1 的匹配 ID转换为 R2 的匹配 ID
seqkit seq -j 12 -n -i "${output_prefix}_R1.fq.gz" | sed 's/\/1.*/\/2/g' > "$matched_ids"
if [ ! -s "$matched_ids" ]; then
    echo "No matched IDs found. Exiting."
    exit 4
fi
echo "Matched IDs extracted."

# 提取 R2 匹配 Reads
seqkit grep -j 12 -f "$matched_ids" "$reads_R2" | $compressor > "${output_prefix}_R2.fq.gz"
if [ $? -ne 0 ] || [ ! -s "${output_prefix}_R2.fq.gz" ]; then
    echo "Error: Failed to extract R2 reads."
    exit 5
fi

echo "Matched R2 reads saved to ${output_prefix}_R2.fq.gz."

