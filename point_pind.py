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

# 输出文件
with open(outputFile, "w") as out:
    out.write("Sample\tChr\tspan at  AK58\tstart at Read\tend at Read\tstart at AK58\tend at AK58\t"
              "strand direction\tspan at JJM\tstart at Read\tend at Read\tstart at JJM\tend at JJM\t"
              "strand direction\tReadID\tread to AK58\n")

    for read in bamfile.fetch(until_eof=True):
        # 跳过未比对的reads
        if read.is_unmapped or not read.has_tag("SA"):
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
            if "D2" not in sa_chrom and "D2" in main_chrom:
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
                    
                    if (-run_site < sml < run_site) and (-run_site < big < run_site):
                        if cigar_sa[-1] == "M":
                            if cigar_main[-1] == "M":
                                out.write(f"{sample}\t{sa_chrom}\t{dic_sa['M']}\t{dic_sa['S'] + 1}\t{read_len}\t"
                                          f"{sa_pos}\t{sa_pos + dic_sa['M'] - 1}\t{sa_strand}\t"
                                          f"{dic_main['M']}\t{dic_main['S'] + 1}\t{read_len}\t"
                                          f"{main_pos}\t{main_pos + dic_main['M'] - 1}\t{strand_main}\t"
                                          f"{read.query_name}\t{read.query_sequence}\n")
                            else:
                                out.write(f"{sample}\t{sa_chrom}\t{dic_sa['M']}\t{dic_sa['S'] + 1}\t{read_len}\t"
                                          f"{sa_pos}\t{sa_pos + dic_sa['M'] - 1}\t{sa_strand}\t"
                                          f"{dic_main['M']}\t1\t{dic_main['M']}\t"
                                          f"{main_pos}\t{main_pos + dic_main['M'] - 1}\t{strand_main}\t"
                                          f"{read.query_name}\t{read.query_sequence}\n")
                        else:
                            if cigar_main[-1] == "M":
                                out.write(f"{sample}\t{sa_chrom}\t{dic_sa['M']}\t1\t{dic_sa['M']}\t"
                                          f"{sa_pos}\t{sa_pos + dic_sa['M'] - 1}\t{sa_strand}\t"
                                          f"{dic_main['M']}\t{dic_main['S'] + 1}\t{read_len}\t"
                                          f"{main_pos}\t{main_pos + dic_main['M'] - 1}\t{strand_main}\t"
                                          f"{read.query_name}\t{read.query_sequence}\n")
                            else:
                                out.write(f"{sample}\t{sa_chrom}\t{dic_sa['M']}\t1\t{dic_sa['M']}\t"
                                          f"{sa_pos}\t{sa_pos + dic_sa['M'] - 1}\t{sa_strand}\t"
                                          f"{dic_main['M']}\t1\t{dic_main['M']}\t"
                                          f"{main_pos}\t{main_pos + dic_main['M'] - 1}\t{strand_main}\t"
                                          f"{read.query_name}\t{read.query_sequence}\n")
                    
bamfile.close()
