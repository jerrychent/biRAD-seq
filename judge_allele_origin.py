# 功能：判断子代等位基因来源
# usage: python judge_allele_origin.py <parent_file> <offspring_file> <output_file>
import sys

def determine_allele_origin(parent_file, offspring_file, output_file):
    # 打开双亲和子代文件
    with open(parent_file, "r") as pf, open(offspring_file, "r") as of, open(output_file, "w") as out:
        # 读取文件头
        parent_header = pf.readline().strip().split("\t")
        offspring_header = of.readline().strip().split("\t")
        
        # 检查文件格式
        if parent_header[:4] != ["#CHROM", "POS", "REF", "ALT"]:
            raise ValueError("双亲文件格式错误，必须包含 #CHROM, POS, REF, ALT 列")
        if offspring_header[:4] != ["#CHROM", "POS", "REF", "ALT"]:
            raise ValueError("子代文件格式错误，必须包含 #CHROM, POS, REF, ALT 列")
        
        # 写入输出文件头
        out.write("#CHROM\tPOS\tREF\tALT\tGENOTYPE\tSOURCE\tQUAL\n")
        
        # 创建字典存储双亲基因型信息
        parent_data = {}
        for line in pf:
            cols = line.strip().split("\t")
            chrom, pos, ref, alt, n22, nip, qual = cols[:7]
            parent_data[(chrom, pos)] = {"REF": ref, "ALT": alt, "N22": n22, "Nip": nip}
        
        # 遍历子代基因型信息
        for line in of:
            cols = line.strip().split("\t")
            chrom, pos, ref, alt, genotype, qual = cols[:6]
            
            # 检查位点是否在双亲文件中
            key = (chrom, pos)
            if key not in parent_data:
                source = "NA"  # 如果双亲中没有对应位点，标记为 NA
            else:
                # 获取双亲基因型信息
                parent_info = parent_data[key]
                n22_allele = parent_info["N22"]
                nip_allele = parent_info["Nip"]
                
                # 判断子代的等位基因来源
                alleles = genotype.split("/")
                if len(alleles) != 2:
                    source = "NA"  # 如果基因型格式不正确，标记为 NA
                else:
                    allele_sources = []
                    for allele in alleles:
                        if allele == ref:
                            if n22_allele == ref:
                                allele_sources.append("A")
                            elif nip_allele == ref:
                                allele_sources.append("B")
                            else:
                                allele_sources.append("NA")
                        elif allele == alt:
                            if n22_allele == alt:
                                allele_sources.append("A")
                            elif nip_allele == alt:
                                allele_sources.append("B")
                            else:
                                allele_sources.append("NA")
                        else:
                            allele_sources.append("NA")  # 如果等位基因与 REF/ALT 不匹配，标记为 NA
                    
                    # 根据等位基因来源判断子代来源
                    if allele_sources[0] == allele_sources[1]:
                        source = allele_sources[0]
                    elif "A" in allele_sources and "B" in allele_sources:
                        source = "H"  # 杂合
                    else:
                        source = "NA"
            
            # 写入输出文件
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{genotype}\t{source}\t{qual}\n")

    print(f"子代等位基因来源判断完成，结果已保存到文件: {output_file}")

if __name__ == "__main__":
    # 检查命令行参数
    if len(sys.argv) != 4:
        print("用法: python determine_allele_origin.py <双亲文件> <子代文件> <输出文件>")
        sys.exit(1)

    # 获取输入文件路径
    parent_file = sys.argv[1]
    offspring_file = sys.argv[2]
    output_file = sys.argv[3]

    # 执行判断
    determine_allele_origin(parent_file, offspring_file, output_file)
