#!/bin/bash

# 检查输入参数
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.bed output.bed"
    exit 1
fi

input=$1
output=$2

# 处理BED文件并生成新的区间
awk '{
    start = $2
    first_end = start + 10
    print $1 "\t" start "\t" first_end
}' "$input" > "$output"