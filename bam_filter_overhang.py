# 依据overhang筛选bam文件，排除样本污染
# usage:python bam_filter_overhang.py sample_tag_list.txt --processes 4
import pysam
import re
import argparse
import os


# 将简并碱基转换为正则表达式
def degenerate_to_regex(tag_sequence, match_start=True):
    degenerate_dict = {
        'A': 'A',
        'T': 'T',
        'C': 'C',
        'G': 'G',
        'N': '[ATCG]'  # N可以匹配A、T、C、G
    }
    regex_pattern = ''.join(degenerate_dict.get(base, base) for base in tag_sequence)
    if match_start:
        return re.compile(f"^{regex_pattern}")  # 匹配以此模式开头
    else:
        return re.compile(f"{regex_pattern}$")  # 匹配以此模式结尾


# 生成反向互补序列
def reverse_complement(seq):
    complement = str.maketrans("ATCGatcg", "TAGCtagc")
    return seq.translate(complement)[::-1]


# 筛选 BAM 文件
def filter_bam_by_tag(input_bam, output_bam, tag_start, tag_end):
    # 转换标签序列为正则表达式
    tag_start_regex = degenerate_to_regex(tag_start, match_start=True)
    tag_end_regex = degenerate_to_regex(tag_end, match_start=False)

    # 打开输入和输出 BAM 文件
    with pysam.AlignmentFile(input_bam, "rb") as in_bam, \
         pysam.AlignmentFile(output_bam, "wb", header=in_bam.header) as out_bam:

        for read in in_bam:
            if read.is_unmapped:
                continue  # 跳过未比对的 reads

            query_sequence = read.query_sequence
            if not query_sequence:
                continue  # 跳过没有序列信息的 reads

            # 检查插入片段长度为正数的 reads
            if read.template_length > 0:  # 正片段长度
                if tag_start_regex.match(query_sequence):  # 检查以指定序列开头
                    out_bam.write(read)

            # 检查插入片段长度为负数的 reads
            elif read.template_length < 0:  # 负片段长度
                if read.is_reverse:
                    # 反向链：直接检查原始序列的结尾
                    if tag_end_regex.search(query_sequence):
                        out_bam.write(read)
                else:
                    # 正向链：对序列取反向互补后检查结尾
                    rc_sequence = reverse_complement(query_sequence)
                    if tag_end_regex.search(rc_sequence):
                        out_bam.write(read)


# 主程序
def main():
    parser = argparse.ArgumentParser(description="Filter BAM reads by specific overhang tags.")
    parser.add_argument("sample_tag_list", help="Sample tag list file (tab-delimited with start and end tags)")
    parser.add_argument("--processes", type=int, default=1, help="Number of parallel processes (default: 1)")

    args = parser.parse_args()

    # 创建子目录 "filtered"（如果不存在）
    output_dir = "filtered"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 解析样本文件和标签序列
    sample_files = []
    tags = {}
    with open(args.sample_tag_list) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            sample_file, tag_start, tag_end = line.split("\t")
            sample_files.append(sample_file)
            tags[sample_file] = (tag_start, tag_end)

    # 逐个样本处理
    for sample_file in sample_files:
        input_bam = f"{sample_file}.split.bam"
        output_bam = os.path.join(output_dir, f"{sample_file}.sorted.bam")  # 输出路径在子目录中
        tag_start, tag_end = tags[sample_file]

        print(f"Processing {sample_file}...")
        filter_bam_by_tag(input_bam, output_bam, tag_start, tag_end)
        print(f"Filtered BAM saved to {output_bam}")


if __name__ == "__main__":
    main()

