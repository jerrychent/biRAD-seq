import argparse
import gzip
from multiprocessing import Pool

def filter_reads_per_sample(input_R1_file, input_R2_file, output_prefix, tag_sequence):
    output_R1_file = output_prefix + '_R1.fq.gz'
    output_R2_file = output_prefix + '_R2.fq.gz'

    with gzip.open(input_R1_file, 'rt') as f_R1, gzip.open(input_R2_file, 'rt') as f_R2, \
         gzip.open(output_R1_file, 'wt') as f_out_R1, gzip.open(output_R2_file, 'wt') as f_out_R2:

        while True:
            # 读取四行，即一个 read 的信息
            read_info_R1 = [next(f_R1) for _ in range(4)]
            read_info_R2 = [next(f_R2) for _ in range(4)]

            # 如果读取到文件末尾，则退出循环
            if not read_info_R1[0] or not read_info_R2[0]:
                break

            # 检查 R1 和 R2 的序列是否以 tag_sequence 开头
            if read_info_R1[1].startswith(tag_sequence) and read_info_R2[1].startswith(tag_sequence):
                # 写入到输出文件中
                f_out_R1.write(''.join(read_info_R1))
                f_out_R2.write(''.join(read_info_R2))

def main():
    parser = argparse.ArgumentParser(description='Filter reads starting with specified tag sequences from input FASTQ files.')
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