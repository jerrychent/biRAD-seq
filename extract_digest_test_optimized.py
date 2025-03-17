# usage: python extract_digest_test_optimized.py <input.bam> <fragments.bed> <output.fq> [num_threads]
import pysam
import pandas as pd
import os
from concurrent.futures import ThreadPoolExecutor
import sys

def load_bed_by_chromosome(fragments_bed):
    """
    按染色体分组加载 BED 文件，并合并连续片段。
    """
    # 加载 BED 文件
    bed_df = pd.read_csv(fragments_bed, sep="\t", header=None, names=["chrom", "start", "end"])
    
    # 合并同一染色体上连续的片段
    merged_intervals = {}
    for chrom, group in bed_df.groupby("chrom"):
        group = group.sort_values("start")
        intervals = []
        current_start, current_end = group.iloc[0]["start"], group.iloc[0]["end"]
        for _, row in group.iterrows():
            if row["start"] <= current_end:  # 如果片段连续
                current_end = max(current_end, row["end"])
            else:
                intervals.append((current_start, current_end))
                current_start, current_end = row["start"], row["end"]
        intervals.append((current_start, current_end))
        merged_intervals[chrom] = intervals
    return merged_intervals

def process_chromosome(args):
    """
    处理单条染色体上的所有片段，提取子序列。
    """
    bam_file, chrom, intervals, high_quality_string = args
    results = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for frag_start, frag_end in intervals:
            fragment_length = frag_end - frag_start
            for read in bam.fetch(chrom, frag_start, frag_end):
                if read.is_unmapped or read.query_sequence is None:
                    continue
                read_start = read.reference_start
                read_end = read_start + read.query_alignment_length
                query_seq = read.query_sequence
                if frag_start >= read_start and frag_end <= read_end:
                    sub_start = frag_start - read_start
                    sub_end = sub_start + fragment_length
                    sub_seq = query_seq[sub_start:sub_end]
                    if len(sub_seq) == fragment_length:
                        results.append(f"@{read.query_name}\n{sub_seq}\n+\n{high_quality_string[:fragment_length]}\n")
    return results

def extract_subreads_with_threads(bam_file, fragments_bed, output_fq, num_threads=8):
    """
    使用多线程从 BAM 文件中提取子序列并写入 FASTQ 文件。
    """
    # 加载并合并 BED 文件片段
    grouped_intervals = load_bed_by_chromosome(fragments_bed)
    
    # 预生成高质量分值字符串
    max_length = max(end - start for intervals in grouped_intervals.values() for start, end in intervals)
    high_quality_string = "I" * max_length  # 高质量分值（ASCII 40）

    # 分配任务
    tasks = []
    for chrom, intervals in grouped_intervals.items():
        tasks.append((bam_file, chrom, intervals, high_quality_string))

    # 使用多线程处理
    results = []
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        for result in executor.map(process_chromosome, tasks):
            results.extend(result)

    # 写入结果到输出文件
    with open(output_fq, "w") as out_fq:
        out_fq.writelines(results)

if __name__ == "__main__":
    # 获取命令行参数
    bam_file = sys.argv[1]
    fragments_bed = sys.argv[2]
    output_fq = sys.argv[3]
    num_threads = int(sys.argv[4]) if len(sys.argv) > 4 else 8  # 默认使用 8 个线程

    # 检查 BAM 文件是否索引
    if not os.path.exists(bam_file + ".bai") and not os.path.exists(bam_file.rsplit(".", 1)[0] + ".csi"):
        raise FileNotFoundError("BAM 文件未索引，请运行 'samtools index' 生成索引文件。")

    # 提取子序列
    extract_subreads_with_threads(bam_file, fragments_bed, output_fq, num_threads)
