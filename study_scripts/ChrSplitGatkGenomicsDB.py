import numpy as np
import pandas as pd
import sys

list_name_txt = sys.argv[1]
chr_list_txt = sys.argv[2]

GATK = "/public/software/env01/bin/gatk"
url_ref = "/public/agis/changyuxiao_group/zhuzhenghang/LH-TAIL-seq-20230111/reference/Zea_mays.AGPv4.dna_sm.toplevel.fa"

chr_list = np.array(pd.read_csv(chr_list_txt,header = None,index_col= None))
list_name = np.array(pd.read_csv(list_name_txt,header = None,index_col= None))
for chr,total_lenth in chr_list:
    chr = str(chr)
    # 如果chr的长度大于10000000，将分成多个文件进行合并
    # 通过调整1000000的值增加进程数
    if total_lenth > 10000000:
        fold = int(np.floor(total_lenth/10000000))
        lenth_list = np.zeros(fold+1)
        for  i in np.arange(fold):
            lenth_list[i] = int(10000000*(i+1))
        lenth_list[fold] = int(total_lenth)

        # 写入shell脚本
        start_po = 1
        for chr_po in lenth_list:
            ContigName = "contig" + chr + "-" + str(int(chr_po))
            shell_txt = ContigName +"_05.sh"
            with open(shell_txt,"ab") as file:
                chr_range = chr + ":" + str(int(start_po)) + "-" + str(int(chr_po)) 

                txt1 = '#!/bin/bash\n' 
                txt2 = GATK + " GenomicsDBImport \\\n"
                txt3 = "-R " + url_ref + " \\" +'\n'

                file.write(txt1.encode()) 
                file.write(txt2.encode())
                file.write(txt3.encode())

                for sample in list_name:
                    txt4 = "-V ./gvcf/" + str(sample[0]) + ".g.vcf.gz " + "\\" +'\n'
                    file.write(txt4.encode()) 

                txt5 = "--batch-size 100" + " \\" +'\n'
                txt6 = "--genomicsdb-workspace-path ./gvcf/" + ContigName + " \\" +'\n'
                txt7 = "-L " + chr_range
                file.write(txt5.encode()) 
                file.write(txt6.encode())
                file.write(txt7.encode())
                start_po = chr_po + 1


        start_po = 1
        for chr_po in lenth_list:
            ContigName = "contig" + chr + "-" + str(int(chr_po))
            shell_txt = ContigName +"_06.sh"
            with open(shell_txt,"ab") as file:
                chr_range = chr + ":" + str(int(start_po)) + "-" + str(int(chr_po)) 

                txt1 = '#!/bin/bash\n' 
                txt2 = GATK + " GenotypeGVCFs \\\n"
                txt3 = "-R " + url_ref + " \\" +'\n'
                txt4 = "-V " + "gendb://gvcf/" + ContigName + " \\" +'\n'
                txt5 = "-O " + "./gvcf/" + ContigName + ".vcf.gz"

                file.write(txt1.encode()) 
                file.write(txt2.encode())
                file.write(txt3.encode())
                file.write(txt4.encode())
                file.write(txt5.encode())


    else:
        shell_txt = chr +"_05.sh"
        with open(shell_txt,"ab") as file:
            txt1 = '#!/bin/bash\n' 
            txt2 = GATK + " GenomicsDBImport \\\n" 
            txt3 = "-R " + url_ref + " \\" +'\n'

            file.write(txt1.encode()) 
            file.write(txt2.encode())
            file.write(txt3.encode())
            for sample in list_name:
                txt4 = "-V ./gvcf/" + str(sample[0]) + ".g.vcf.gz " + "\\" +'\n'
                file.write(txt4.encode()) 

            txt5 = "--batch-size 100" + " \\" +'\n'
            txt6 = "--genomicsdb-workspace-path ./gvcf/" + chr + " \\" + '\n'
            txt7 = "-L " + chr
            file.write(txt5.encode()) 
            file.write(txt6.encode())
            file.write(txt7.encode())


        shell_txt = chr +"_06.sh"
        with open(shell_txt,"ab") as file:
            chr_range = chr + ":" + str(int(start_po)) + "-" + str(int(chr_po)) 

            txt1 = '#!/bin/bash\n' 
            txt2 = GATK + " GenotypeGVCFs \\\n"
            txt3 = "-R " + url_ref + " \\" +'\n'
            txt4 = "-V " + "gendb://gvcf/" + chr + " \\" + '\n'
            txt5 = "-O " + "./gvcf/" + chr + ".vcf.gz"

            file.write(txt1.encode()) 
            file.write(txt2.encode())
            file.write(txt3.encode())
            file.write(txt4.encode())
            file.write(txt5.encode())