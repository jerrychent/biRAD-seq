# 筛选与亲本位点相同的 SNP。
# example ： python3 new_zd_filter.py qd_file.path zd_file.path out.file
import sys

print("CHR", "POS", "REF", "ALT", "GENOTYPE", "genotype1", "genotype2", sep="\t")


def genotype():
    qd_file = sys.argv[1]
    zd_file = sys.argv[2]
    qf = open(qd_file, mode='r')
    zf = open(zd_file, mode='r')
    dic = {}
    for qd_line in qf:
        if qd_line.startswith("#"):
            continue
        qd_sam = qd_line.strip("\r\n").split("\t")
        qds = (qd_sam[0], qd_sam[1])
        kv = {qds: qd_line}
        dic.update(kv)
    for zd_line in zf:
        if zd_line.startswith("#"):
            continue
        zd_sam = zd_line.strip("\r\n").split("\t")
        z_sam10 = zd_sam[9].split(":")[0]
        zds = (zd_sam[0], zd_sam[1])
        if zds in dic:
            qd_line_n = dic[zds]
            qd_sam_n = qd_line_n.strip("\r\n").split("\t")
            q_sam10 = qd_sam_n[9].split(":")[0]
            q_sam11 = qd_sam_n[10].split(":")[0]
            if zd_sam[4] == qd_sam_n[4]:
                if z_sam10 == q_sam10:
                    print(zd_sam[0], zd_sam[1], zd_sam[3], zd_sam[4], zd_sam[9].split(":")[0], "0", "p1", sep="\t")
                elif z_sam10 == q_sam11:
                    print(zd_sam[0], zd_sam[1], zd_sam[3], zd_sam[4], zd_sam[9].split(":")[0], "2", "P2", sep="\t")
                else:
                    print(zd_sam[0], zd_sam[1], zd_sam[3], zd_sam[4], zd_sam[9].split(":")[0], "1", "H", sep="\t")

    qf.close()
    zf.close()


if __name__ == "__main__":
    genotype()

