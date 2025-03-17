import gzip
import sys

def read_fastq(file):
    with gzip.open(file, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            sequence = f.readline().strip()
            f.readline()  # Plus line
            f.readline()  # Quality line
            yield header, sequence

def count_matching_read_pairs(sub_r1_file, sub_r2_file, r1_file, r2_file):
    sub_r1_reads = {header: seq for header, seq in read_fastq(sub_r1_file)}
    sub_r2_reads = {header: seq for header, seq in read_fastq(sub_r2_file)}

    matching_pairs = 0
    total_pairs = 0

    r1_iter = read_fastq(r1_file)
    r2_iter = read_fastq(r2_file)

    for (header_r1, seq_r1), (header_r2, seq_r2) in zip(r1_iter, r2_iter):
        total_pairs += 1
        if (header_r1 in sub_r1_reads and sub_r1_reads[header_r1] == seq_r1) and \
           (header_r2 in sub_r2_reads and sub_r2_reads[header_r2] == seq_r2):
            matching_pairs += 1

    return matching_pairs, total_pairs

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <1.sub_R1.fq.gz> <1.sub_R2.fq.gz> <1_R1.fq.gz> <1_R2.fq.gz>")
        sys.exit(1)

    sub_r1_file = sys.argv[1]
    sub_r2_file = sys.argv[2]
    r1_file = sys.argv[3]
    r2_file = sys.argv[4]

    matching_pairs, total_pairs = count_matching_read_pairs(sub_r1_file, sub_r2_file, r1_file, r2_file)

    matching_percentage = (matching_pairs / total_pairs) * 100 if total_pairs > 0 else 0

    print(f'Matching read pairs: {matching_pairs} out of {total_pairs} ({matching_percentage:.2f}%)')
