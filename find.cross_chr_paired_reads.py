# usage: python find.cross_chr_reads.py input.bam 20 output.txt
import sys
import os
import re
import pysam

def is_strand(flag):
    """根据 FLAG 字段确定比对链方向"""
    return "-" if (flag & 16) != 0 else "+"

# 检查输入参数
if len(sys.argv) != 4:
    print("Usage: python script.py <input.bam> <run_site> <output.txt>")
    sys.exit(1)

inputFile = sys.argv[1]
run_site = int(sys.argv[2])
outputFile = sys.argv[3]

# 检查输入文件是否存在
if not os.path.isfile(inputFile):
    print(f"Error: File '{inputFile}' not found.")
    sys.exit(1)

# 打开 BAM 文件
bamfile = pysam.AlignmentFile(inputFile, "rb")
sample = os.path.basename(inputFile).split(".")[0]

# 用于存储输出的列表
output_lines = []

for read in bamfile.fetch(until_eof=True):
    # 跳过未比对的 reads 或没有配对信息的 reads
    if read.is_unmapped or not read.is_paired:
        continue

    # 筛除插入片段大于 1000bp 的 reads
    if abs(read.template_length) > 1000:
        continue

    # 配对 read 信息
    mate_chrom = read.next_reference_name  # 配对 read 的染色体名称

    # 配对 read 位于主要比对染色体上或在其次要比对所在的染色体上
    if mate_chrom != "=" and mate_chrom != read.reference_name:
        continue

    # 跳过未比对的 reads 或没有 SA 标签的 reads
    if not read.has_tag("SA"):
        continue

    # 主要比对信息
    main_chrom = read.reference_name
    main_pos = read.reference_start + 1  # 转为1基坐标
    cigar_main = read.cigarstring
    strand_main = is_strand(read.flag)

    # 次要比对信息
    sa_tag = read.get_tag("SA")
    sa_fields = sa_tag.split(";")

    for sa_entry in sa_fields:
        if not sa_entry.strip():
            continue
        sa_info = sa_entry.split(",")
        sa_chrom = sa_info[0]
        sa_pos = int(sa_info[1])
        sa_strand = sa_info[2]
        cigar_sa = sa_info[3]

        # 筛选条件：主要比对在不含“D2”字符的染色体，次要比对在含“D2”字符的染色体
        if "D2" not in main_chrom and "D2" in sa_chrom:
            # 检查配对 read 的次要比对
            mate_secondary_chroms = set()
            if read.mate_is_unmapped:
                continue
            mate = bamfile.mate(read)
            if mate.has_tag("SA"):
                mate_sa_tag = mate.get_tag("SA")
                mate_sa_fields = mate_sa_tag.split(";")
                for mate_sa_entry in mate_sa_fields:
                    if not mate_sa_entry.strip():
                        continue
                    mate_sa_info = mate_sa_entry.split(",")
                    mate_sa_chrom = mate_sa_info[0]
                    mate_secondary_chroms.add(mate_sa_chrom)

            # 如果配对 read 存在次要比对，要求其次要比对染色体与当前 read 的次要比对染色体一致
            if mate_secondary_chroms and sa_chrom not in mate_secondary_chroms:
                continue

            # 解析 CIGAR 字段
            seg_main = re.findall(r"(\d+)([A-Z])", cigar_main)
            seg_sa = re.findall(r"(\d+)([A-Z])", cigar_sa)

            # 构建比对片段长度字典
            dic_main = {}
            for value, key in seg_main:
                if key not in dic_main:
                    dic_main[key] = int(value)
            if "H" in dic_main:
                dic_main["S"] = dic_main.pop("H")

            dic_sa = {}
            for value, key in seg_sa:
                if key not in dic_sa:
                    dic_sa[key] = int(value)
            if "H" in dic_sa:
                dic_sa["S"] = dic_sa.pop("H")

            # 计算 read 长度
            if "M" in dic_main and "S" in dic_main and "M" in dic_sa and "S" in dic_sa:
                read_len = dic_main["M"] + dic_main["S"]
                sml = dic_main["M"] - dic_sa["S"]
                big = dic_main["S"] - dic_sa["M"]

                # 筛选条件：长度差异在范围内
                if (-run_site < sml < run_site) and (-run_site < big < run_site):
                    # 筛选条件：主要比对和次要比对的 CIGAR 末尾字符必须不同
                    if (cigar_main[-1] == "M" and cigar_sa[-1] == "S") or (cigar_main[-1] == "S" and cigar_sa[-1] == "M"):
                        # 根据 CIGAR 值末尾字符生成输出行
                        if cigar_sa[-1] == "M":
                            if cigar_main[-1] == "M":
                                # 这种情况不会发生，因为末尾字符必须不同
                                continue
                            else:
                                # 次要比对以 M 结尾，主要比对以 S 结尾
                                line = f"{sample}\t{main_chrom}\t{dic_main['M']}\t1\t{dic_main['M']}\t" \
                                    f"{main_pos}\t{main_pos + dic_main['M'] - 1}\t{strand_main}\t" \
                                    f"{sa_chrom}\t{dic_sa['M']}\t{dic_sa['S'] + 1}\t{read_len}\t" \
                                    f"{sa_pos}\t{sa_pos + dic_sa['M'] - 1}\t{sa_strand}\t" \
                                    f"{read.query_name}\t{read.query_sequence}\n"
                        else:
                            if cigar_main[-1] == "M":
                                # 主要比对以 M 结尾，次要比对以 S 结尾
                                line = f"{sample}\t{main_chrom}\t{dic_main['M']}\t{dic_main['S'] + 1}\t{read_len}\t" \
                                    f"{main_pos}\t{main_pos + dic_main['M'] - 1}\t{strand_main}\t" \
                                    f"{sa_chrom}\t{dic_sa['M']}\t1\t{dic_sa['M']}\t" \
                                    f"{sa_pos}\t{sa_pos + dic_sa['M'] - 1}\t{sa_strand}\t" \
                                    f"{read.query_name}\t{read.query_sequence}\n"
                            else:
                                # 这种情况不会发生，因为末尾字符必须不同
                                continue

                        # 保存符合条件的输出行
                        output_lines.append((main_chrom, main_pos, line))

bamfile.close()

# 对输出数据进行排序
# 按照第二列（染色体名称）升序排列，再按照第六列（主要比对开始位置）升序排列
output_lines.sort(key=lambda x: (x[0], x[1]))

# 写入排序后的结果到文件
with open(outputFile, "w") as out:
    out.write("Sample\tChr\tspan at  AK58\tstart at Read\tend at Read\tstart at AK58\tend at AK58\t"
              "strand direction\tSecondary Chromosome\tspan at JJM\tstart at Read\tend at Read\tstart at JJM\tend at JJM\t"
              "strand direction\tReadID\tread to AK58\n")
    for _, _, line in output_lines:
        out.write(line)
