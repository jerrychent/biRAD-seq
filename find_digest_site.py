# python find_digest_site.py -g reference_genome.fa -e TTAA -min 250 -max 350 -o fragments.bed
import argparse
import re
from Bio import SeqIO

def create_regex(enzyme_seq):
    """将含有简并碱基的酶切位点转为正则表达式"""
    iupac_dict = {
        "A": "A", "T": "T", "C": "C", "G": "G",
        "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
        "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
        "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"
    }
    return "".join(iupac_dict[base] for base in enzyme_seq.upper())

def find_re_sites(genome_file, enzyme_regex):
    """在基因组文件中寻找酶切位点"""
    positions = {}
    compiled_regex = re.compile(enzyme_regex, re.IGNORECASE)  # 正则表达式忽略大小写

    with open(genome_file, "r") as genome:
        for record in SeqIO.parse(genome, "fasta"):
            seq = str(record.seq)
            chr_name = record.id
            positions[chr_name] = []
            
            # 搜索正则表达式匹配的位点
            for match in compiled_regex.finditer(seq):
                positions[chr_name].append(match.start() + 1)  # 使用1-based位置
    return positions

def filter_fragments(re_sites, min_len, max_len):
    """筛选符合长度要求的片段"""
    fragments = []
    for chr_name, sites in re_sites.items():
        for i in range(len(sites) - 1):
            start = sites[i]
            end = sites[i + 1]
            fragment_len = end - start
            if min_len <= fragment_len <= max_len:
                fragments.append((chr_name, start, end))
    return fragments

def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="模拟酶切建库生成简化基因组测序片段")
    parser.add_argument("-g", "--genome", required=True, help="参考基因组文件（FASTA 格式）")
    parser.add_argument("-e", "--enzyme", required=True, help="酶切识别位点序列（支持简并碱基）")
    parser.add_argument("-min", "--min_length", type=int, required=True, help="片段最小长度")
    parser.add_argument("-max", "--max_length", type=int, required=True, help="片段最大长度")
    parser.add_argument("-o", "--output", default="fragments.bed", help="输出 BED 文件名称")
    
    args = parser.parse_args()
    
    # 转换酶切位点为正则表达式
    enzyme_regex = create_regex(args.enzyme)
    
    # 读取参考基因组并寻找酶切位点
    print("正在寻找酶切位点...")
    re_sites = find_re_sites(args.genome, enzyme_regex)
    
    # 筛选符合条件的片段
    print(f"筛选片段，长度范围为 {args.min_length} - {args.max_length} bp...")
    filtered_fragments = filter_fragments(re_sites, args.min_length, args.max_length)
    
    # 保存片段到 BED 文件
    with open(args.output, "w") as bed_file:
        for chr_name, start, end in filtered_fragments:
            bed_file.write(f"{chr_name}\t{start}\t{end}\n")
    
    print(f"完成！共找到 {len(filtered_fragments)} 个符合条件的片段，结果已保存到 {args.output}")

if __name__ == "__main__":
    main()
