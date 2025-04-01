# 从vcf文件中提取双亲的基因型信息
# usage: python extract_parents_genotypes.py <input_vcf> <output_file> 
import sys
from cyvcf2 import VCF

def extract_genotypes(input_vcf, output_file):
    # 打开 VCF 文件
    vcf = VCF(input_vcf)

    # 获取样本名称
    samples = vcf.samples
    if "N22" not in samples or "Nip" not in samples:
        raise ValueError("VCF 文件中未找到 N22 或 Nip 样本，请检查样本名称是否正确！")

    # 确定 N22 和 Nip 的索引
    n22_index = samples.index("N22")
    nip_index = samples.index("Nip")

    # 打开输出文件
    with open(output_file, "w") as out:
        # 写入文件头
        out.write("#CHROM\tPOS\tREF\tALT\tN22\tNip\tQUAL\n")
        
        # 遍历每个变异位点
        for variant in vcf:
            chrom = variant.CHROM
            pos = variant.POS
            ref = variant.REF
            alt = variant.ALT[0]  # 假设只有一个 ALT 等位基因
            qual = variant.QUAL
            
            # 获取 N22 和 Nip 的基因型
            n22_genotype = variant.genotypes[n22_index][:2]  # 前两位是基因型
            nip_genotype = variant.genotypes[nip_index][:2]
            
            # 基因型转为等位基因
            n22_allele = ref if n22_genotype[0] == 0 else alt
            nip_allele = ref if nip_genotype[0] == 0 else alt
            
            # 写入输出文件
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{n22_allele}\t{nip_allele}\t{qual}\n")

    print(f"基因型信息已提取到文件: {output_file}")

if __name__ == "__main__":
    # 检查命令行参数
    if len(sys.argv) != 3:
        print("用法: python extract_genotypes.py <输入VCF文件> <输出文件>")
        sys.exit(1)

    # 获取输入和输出文件路径
    input_vcf = sys.argv[1]
    output_file = sys.argv[2]

    # 提取基因型信息
    extract_genotypes(input_vcf, output_file)
