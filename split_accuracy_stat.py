#!/usr/bin/python
import pysam
import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Compare read names between FASTQ and BAM files.')
    parser.add_argument('-fq1', '--fastq1', required=True, help='Path to the first FASTQ file')
    parser.add_argument('-fq2', '--fastq2', required=True, help='Path to the second FASTQ file')
    parser.add_argument('-b', '--bam', required=True, help='Path to BAM file')
    parser.add_argument('-o', '--output', required=True, help='Path to output file for false positive reads')
    return parser.parse_args()

def filter_false_positives(bamfile, fq1_names, fq2_names, output_file):
    false_positives = set()
    with pysam.AlignmentFile(bamfile, 'rb') as bamfile:
        for read in bamfile.fetch():
            read_name = read.query_name
            if read_name not in fq1_names and read_name not in fq2_names:
                false_positives.add(read_name)
    with open(output_file, 'w') as f:
        for read in sorted(false_positives):
            f.write(read + '\n')

def main():
    args = get_args()

    # 存储read名称的集合
    read_names_fq1 = set()
    read_names_fq2 = set()

    # 读取第一个FASTQ文件的read名称
    with pysam.FastxFile(args.fastq1) as fastx_file:
        for entry in fastx_file:
            read_name = entry.name[0:].split("/")[0]  # 去除开头的@和结尾的/1
            read_names_fq1.add(read_name)

    # 读取第二个FASTQ文件的read名称
    with pysam.FastxFile(args.fastq2) as fastx_file:
        for entry in fastx_file:
            read_name = entry.name[0:].split("/")[0]  # 去除开头的@和结尾的/1
            read_names_fq2.add(read_name)

    # 合并并去重 read 名称
    all_read_names = read_names_fq1.union(read_names_fq2)

    # 获取假阳性 read
    filter_false_positives(args.bam, read_names_fq1, read_names_fq2, args.output)

    # 打开BAM文件
    with pysam.AlignmentFile(args.bam, 'rb') as bamfile:
        # 初始化计数器
        count_fq1_reads = len(read_names_fq1)
        count_fq2_reads = len(read_names_fq2)
        count_fq_reads = 0
        total_bam_reads = 0

        # 遍历BAM文件中的reads
        for read in bamfile.fetch():
            # 获取read的名称
            read_name = read.query_name
            total_bam_reads += 1

            # 判断read来自哪个文件
            if read_name in all_read_names:
                count_fq_reads += 1

    # 计算比例
    if total_bam_reads > 0:
        ratio = count_fq_reads / total_bam_reads
    else:
        ratio = 0

    # 打印结果
    print(f'Reads from {args.fastq1}: {count_fq1_reads}')
    print(f'Reads from {args.fastq2}: {count_fq2_reads}')
    print(f'Total reads from BAM file: {total_bam_reads}')
    print(f'Reads from merged FASTQ files: {count_fq_reads}')
    print(f'Ratio of reads from FASTQ files in BAM file: {ratio:.4f}')

if __name__ == "__main__":
    main()
