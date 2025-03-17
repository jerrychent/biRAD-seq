# usage: python split_overhang.py -i input_R1.fq.gz -o output_R1.fq.gz -t ATCG -p 4
# split_test_multi_sample.py
# 多样本多线程同时拆分数据
import argparse
import gzip
import re
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

def filter_reads_per_sample(input_R1_file, input_R2_file, output_prefix, tag_sequence):
    output_R1_file = output_prefix + '_R1.fq.gz'
    output_R2_file = output_prefix + '_R2.fq.gz'

    tag_regex = degenerate_to_regex(tag_sequence)  # 将tag_sequence转换为正则表达式

    with gzip.open(input_R1_file, 'rt') as f_R1, gzip.open(input_R2_file, 'rt') as f_R2, \
         gzip.open(output_R1_file, 'wt') as f_out_R1, gzip.open(output_R2_file, 'wt') as f_out_R2:

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
    parser = argparse.ArgumentParser(description='Filter reads starting with specified tag sequences (including degenerate bases) from input FASTQ files.')
    parser.add_argument('input_R1', help='Input R1 FASTQ file (gzip compressed)')
    parser.add_argument('input_R2', help='Input R2 FASTQ file (gzip compressed)')
    parser.add_argument('output_list', help='List of output file prefixes and tag sequences (tab separated)')
    parser.add_argument('--processes', type=int, default=1, help='Number of processes to use (default: 1)')

    args = parser.parse_args()

    # 解析输出文件前缀和tag序列
    output_prefixes = []
    tag_sequences = {}
    with open(args.output_list) as f:
        for line in f:
            line = line.strip()  # 去除行末尾的空白字符
            if not line:
                continue  # 跳过空行
            prefix, sequence = line.split('\t')
            output_prefixes.append(prefix)
            tag_sequences[prefix] = sequence

    # 并行处理多个样本
    with Pool(args.processes) as pool:
        pool.starmap(filter_reads_per_sample, [(args.input_R1, args.input_R2, prefix, tag_sequences[prefix]) for prefix in output_prefixes])

if __name__ == "__main__":
    main()

