#!/bin/bash

# 根目录（包含Lib1~Lib19）
base_dir="/public/agis/changyuxiao_group/chenhaotian/mywork/2024.11.14-rice_RIL/data_sequencing_test/RIL-ELF6/lib_data_split/rename"

# 初始样本编号
sample_id=1

# 遍历每个子目录
for lib_dir in "$base_dir"/Lib*; do
    # 确保是目录
    if [ -d "$lib_dir" ]; then
        echo "Processing directory: $lib_dir"

        # 获取当前目录下所有文件前缀，按字母顺序排序
        prefixes=$(ls "$lib_dir"/*_R1.fq.gz | sed 's/.*\///' | sed 's/_R1\.fq\.gz//' | sort)

        # 遍历每个前缀，按顺序重命名
        for prefix in $prefixes; do
            # 构造新前缀，修改为 RIL-前缀格式
            new_prefix="RIL-$(printf "%d" $sample_id)"  # 例如 RIL-1, RIL-2

            # 获取文件路径
            r1_file="$lib_dir/${prefix}_R1.fq.gz"
            r2_file="$lib_dir/${prefix}_R2.fq.gz"

            # 新文件路径
            new_r1_file="$lib_dir/${new_prefix}_R1.fq.gz"
            new_r2_file="$lib_dir/${new_prefix}_R2.fq.gz"

            # 重命名
            mv "$r1_file" "$new_r1_file"
            mv "$r2_file" "$new_r2_file"

            echo "Renamed: $prefix -> $new_prefix"

            # 更新样本编号
            ((sample_id++))
        done
    fi
done

echo "Renaming completed."

