# 依据overhang筛选fq文件，排除样本污染
# usage:python fq_filter_overhang.py sample_tag_list.txt --processes 4

import argparse
import gzip
import re
import os
from multiprocessing import Pool

# 将简并碱基转换为正则表达式模式
def degenerate_to_regex(tag_sequence):
    degenerate_dict = {
        'A': 'A',
        'T': 'T',
        'C': 'C',
        'G': 'G',
        'N': '[ATCG]'  # N可以匹配A、T、C、G
    }
    regex_pattern = ''.join(degenerate_dict.get(base, base) for base in tag_sequence)
    return re.compile(f"^{regex_pattern}")  # 匹配以此正则模式开头的序列

def filter_reads_per_sample(input_file, output_prefix, tag_sequence):
    # 创建子目录 fq_filter_overhang（如果不存在）
    output_dir = 'fq_filter_overhang'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_R1_file = os.path.join(output_dir, output_prefix + '_R1.fq.gz')
    output_R2_file = os.path.join(output_dir, output_prefix + '_R2.fq.gz')

    tag_regex = degenerate_to_regex(tag_sequence)  # 将tag_sequence转换为正则表达式

    with gzip.open(input_file + '_R1.fq.gz', 'rt') as f_R1, \
         gzip.open(input_file + '_R2.fq.gz', 'rt') as f_R2, \
         gzip.open(output_R1_file, 'wt') as f_out_R1, \
         gzip.open(output_R2_file, 'wt') as f_out_R2:

        while True:
            try:
                # 读取四行，即一个read的信息
                read_info_R1 = [next(f_R1) for _ in range(4)]
                read_info_R2 = [next(f_R2) for _ in range(4)]
            except StopIteration:
                break  # 读取到文件末尾，退出循环

            # 检查R1和R2的序列是否匹配tag_regex
            if tag_regex.match(read_info_R1[1]) and tag_regex.match(read_info_R2[1]):
                # 写入到输出文件中
                f_out_R1.write(''.join(read_info_R1))
                f_out_R2.write(''.join(read_info_R2))

def main():
    parser = argparse.ArgumentParser(description='Filter reads from multiple sample FASTQ files starting with specified tag sequences.')
    parser.add_argument('sample_tag_list', help='List of sample file names and corresponding tag sequences (tab separated)')
    parser.add_argument('--processes', type=int, default=1, help='Number of processes to use (default: 1)')

    args = parser.parse_args()

    # 解析样本文件和tag序列
    sample_files = []
    tag_sequences = {}
    with open(args.sample_tag_list) as f:
        for line in f:
            line = line.strip()  # 去除行末尾的空白字符
            if not line:
                continue  # 跳过空行
            sample_file, sequence = line.split('\t')
            sample_files.append(sample_file)
            tag_sequences[sample_file] = sequence

    # 并行处理多个样本
    with Pool(args.processes) as pool:
        pool.starmap(filter_reads_per_sample, [(sample, sample, tag_sequences[sample]) for sample in sample_files])

if __name__ == "__main__":
    main()
