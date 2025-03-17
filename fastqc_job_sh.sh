#!/bin/bash

# 输入数据所在目录
input_dir="./02_Cleandata"
# 输出目录，用于保存生成的job.sh脚本
output_dir="./02_Cleandata/fastqc/jobs"

# 创建输出目录（如果不存在）
mkdir -p $output_dir

# 遍历input_dir目录下所有.fastq.gz文件
for fastq_file in $input_dir/*.fq.gz; do
    # 获取文件名（不带路径）
    file_name=$(basename $fastq_file)
    # 去掉扩展名，获取基本文件名
    base_name=${file_name%.fq.gz}
    
    # 定义job脚本的路径
    job_file="$output_dir/job_$base_name.sh"

    # 生成job.sh脚本
    echo "#!/bin/bash" > $job_file
    echo "fastqc -o /public/agis/changyuxiao_group/chenhaotian/2024.9.12-wheat_AK58/02_Cleandata/fastqc -t 2 /public/agis/changyuxiao_group/chenhaotian/2024.9.12-wheat_AK58/02_Cleandata/$file_name" >> $job_file  # fastqc命令

    # 确保job.sh脚本可执行
    chmod +x $job_file

    echo "生成作业脚本: $job_file"
done

echo "所有job.sh脚本生成完毕！"
