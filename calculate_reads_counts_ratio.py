# usage:python calculate_ratio.py coverage_before.txt coverage_after.txt

import sys

def read_coverage(file):
    coverage_dict = {}
    with open(file, 'r') as f:
        for line in f:
            fields = line.strip().split()
            chrom = fields[0]
            start = fields[1]
            end = fields[2]
            reads_count = int(fields[3])
            key = (chrom, start, end)
            coverage_dict[key] = reads_count
    return coverage_dict

def calculate_ratio(before_file, after_file):
    before_coverage = read_coverage(before_file)
    after_coverage = read_coverage(after_file)

    print("chrom\tstart\tend\tbefore_reads\tafter_reads\tratio")
    for key in before_coverage:
        before_reads = before_coverage[key]
        after_reads = after_coverage.get(key, 0)  # 如果没有匹配的窗口，设为0
        if before_reads > 0:  # 避免除以0的情况
            ratio = after_reads / before_reads
        else:
            ratio = 'NA'
        print(f"{key[0]}\t{key[1]}\t{key[2]}\t{before_reads}\t{after_reads}\t{ratio}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <coverage_before.txt> <coverage_after.txt>")
        sys.exit(1)

    before_file = sys.argv[1]
    after_file = sys.argv[2]

    calculate_ratio(before_file, after_file)

