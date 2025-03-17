#!/bin/bash

# 检查是否提供了足够的参数
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_vcf> <output_vcf>"
    exit 1
fi

# 输入和输出文件
input_vcf="$1"
output_vcf="$2"

# 临时文件
temp_header="temp_header.vcf"
temp_filtered="temp_filtered.vcf"

# 提取 VCF 文件头部信息
grep '^#' "$input_vcf" > "$temp_header"

# 处理 VCF 文件
grep -v '^#' "$input_vcf" | awk '
{
    filter_flag=0
    for (i=10; i<=NF; i++) {
        if ($i ~ /\.\/\./) { # 未定型的位点
            filter_flag=1
            break
        }
        if ($i ~ /0\/0/) { # 纯合参考基因型
            filter_flag=1
            break
        }
        if ($i ~ /\|/) { # 具有相位信息的基因型
            filter_flag=1
            break
        }
    }
    if (filter_flag == 0 && $1 !~ /^#/ && $1 ~ /^Chr/) {
        print $0
    }
}' > "$temp_filtered"

# 合并头部和过滤后的内容
cat "$temp_header" "$temp_filtered" > "$output_vcf"

# 删除临时文件
rm "$temp_header" "$temp_filtered"

echo "过滤完成，结果保存在 $output_vcf"

