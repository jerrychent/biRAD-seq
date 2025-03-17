#!/bin/bash

# 检查是否提供了足够的参数
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_fasta> <output_fasta>"
    exit 1
fi

# 输入和输出文件
input_fasta="$1"
output_fasta="$2"

# 使用 awk 过滤掉 contig 序列，只保留以 "chr" 开头的序列
awk '
    BEGIN {output = 0}
    /^>/ {
        if ($0 ~ /^>chr/) {
            output = 1
            print $0
        } else {
            output = 0
        }
    }
    output == 1 && $0 !~ /^>/ {print $0}
' "$input_fasta" > "$output_fasta"

echo "过滤完成，结果保存在 $output_fasta"

