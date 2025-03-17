#!/bin/bash

# 检查参数是否正确
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <parents_file> <A1_file> <output_file>"
  exit 1
fi

# 读取输入文件
parents_file=$1
A1_file=$2
output_file=$3

# 计算比值并输出结果，忽略大小写并跳过第一行（假设是header）
awk 'FNR==1 {next} FNR==NR {parents[tolower($1),$2,$3] = $6; next} (tolower($1),$2,$3) in parents && parents[tolower($1),$2,$3] != 0 {ratio = $6 / parents[tolower($1),$2,$3]; print $1, $2, $3, ratio}' OFS="\t" "$parents_file" "$A1_file" > "$output_file"

echo "Results saved to $output_file"
