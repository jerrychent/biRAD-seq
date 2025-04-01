# 从vcf文件中提取基因型信息
# usage: python extract_single_sample_genotypes.py <输入VCF文件> <输出文件> <样本名称>
import sys
from cyvcf2 import VCF

def extract_single_sample_genotypes(input_vcf, output_file, sample_name):
    # 打开 VCF 文件
    vcf = VCF(input_vcf)

    # 获取样本名称
    samples = vcf.samples
    if sample_name not in samples:
        raise ValueError(f"VCF 文件中未找到样本 {sample_name}，请检查样本名称是否正确！")

    # 确定样本的索引
    sample_index = samples.index(sample_name)

    # 打开输出文件
    with open(output_file, "w") as out:
        # 写入文件头
        out.write("#CHROM\tPOS\tREF\tALT\tGENOTYPE\tQUAL\n")
        
        # 遍历每个变异位点
        for variant in vcf:
            chrom = variant.CHROM
            pos = variant.POS
            ref = variant.REF
            alt = variant.ALT[0] if variant.ALT else "."  # 如果没有 ALT，则设为 "."
            qual = variant.QUAL
            
            # 获取样本的基因型
            sample_genotype = variant.genotypes[sample_index][:2]  # 前两位是基因型
            
            # 基因型转为等位基因
            genotype_alleles = []
            for allele in sample_genotype:
                if allele == 0:
                    genotype_alleles.append(ref)
                elif allele == 1:
                    genotype_alleles.append(alt)
                else:
                    genotype_alleles.append(".")  # 缺失基因型
            
            genotype_str = "/".join(genotype_alleles)
            
            # 写入输出文件
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{genotype_str}\t{qual}\n")

    print(f"样本 {sample_name} 的基因型信息已提取到文件: {output_file}")

if __name__ == "__main__":
    # 检查命令行参数
    if len(sys.argv) != 4:
        print("用法: python extract_single_sample_genotypes.py <输入VCF文件> <输出文件> <样本名称>")
        sys.exit(1)

    # 获取输入和输出文件路径，以及样本名称
    input_vcf = sys.argv[1]
    output_file = sys.argv[2]
    sample_name = sys.argv[3]

    # 提取单个样本的基因型信息
    extract_single_sample_genotypes(input_vcf, output_file, sample_name)
