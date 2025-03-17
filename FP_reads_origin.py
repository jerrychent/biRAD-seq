#!/usr/bin/python
import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description='Classify false positive reads to their original samples.')
    parser.add_argument('-fp', '--false_positive', required=True, help='Path to the file containing false positive read names')
    parser.add_argument('-fq', '--fastq_list', required=True, help='Path to the file containing a list of FASTQ files for each sample')
    parser.add_argument('-o', '--output_prefix', required=True, help='Output prefix for the result files')
    return parser.parse_args()
    

def read_false_positive_names(fp_file):
    false_positives = set()
    with open(fp_file, 'r') as f:
        for line in f:
            false_positives.add(line.strip())
    return false_positives
    
def read_fastq_list(fq_list_file):
    fastq_list = {}
    with open(fq_list_file, 'r') as f:
        for line in f:
            sample, fastq_file = line.strip().split()
            if sample not in fastq_list:
                fastq_list[sample] = []
            fastq_list[sample].append(fastq_file)
    return fastq_list
    
def classify_false_positives(fp_names, fastq_list, output_prefix):
    for sample, fastq_files in fastq_list.items():
        sample_fp_reads = set()
        for fastq_file in fastq_files:
            # 解压缩 .fq.gz 文件并以文本模式打开
            with gzip.open(fastq_file, 'rt') as f:
                for line in f:
                    if line.startswith('@'):
                        read_name = line.strip()[1:].split()[0]
                        for fp_name in fp_names:
                            if fp_name in read_name:
                                sample_fp_reads.add(read_name)
                                break  # 如果找到匹配的假阳性 read 名称，则跳出内层循环
        output_file = f'{output_prefix}_{sample}_FP_reads.txt'
        with open(output_file, 'w') as f:
            for read_name in sample_fp_reads:
                f.write(read_name + '\n')
    
def main():
    args = get_args()

    # 读取假阳性 read 名称和 FASTQ 文件列表
    false_positive_names = read_false_positive_names(args.false_positive)
    fastq_list = read_fastq_list(args.fastq_list)

    # 对假阳性 read 进行分类
    classify_false_positives(false_positive_names, fastq_list, args.output_prefix)

    print('False positive reads have been classified successfully.')

if __name__ == "__main__":
    main()

