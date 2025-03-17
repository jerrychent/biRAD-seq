## usage: python 20241009_split_fastq.py --fq_r1 R1.fastq.gz --fq_r2 R2.fastq.gz --barcodes barcodes.txt --output_r1 R1_filtered.fastq.gz --output_r2 R2_filtered.fastq.gz --threads 8

import gzip
import argparse
from multiprocessing import Pool, cpu_count

def read_fq(file):
    with gzip.open(file, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            yield header, seq, plus, qual

def process_r1(fq_r1, barcodes):
    r1_dict = {}
    for header, seq, plus, qual in read_fq(fq_r1):
        for barcode in barcodes:
            if seq.startswith(barcode):
                read_id = header.split()[0]
                r1_dict[read_id] = (header, seq, plus, qual)
                break
    return r1_dict

def process_r2(fq_r2, r1_dict):
    r2_reads = []
    for header, seq, plus, qual in read_fq(fq_r2):
        read_id = header.split()[0]
        if read_id in r1_dict:
            r2_reads.append((header, seq, plus, qual))
    return r2_reads

def write_fq(output_file, reads):
    with gzip.open(output_file, 'wt') as f:
        for read in reads:
            f.write('\n'.join(read) + '\n')

def main(fq_r1, fq_r2, barcode_file, output_r1, output_r2, threads):
    with open(barcode_file, 'r') as f:
        barcodes = [line.strip() for line in f]

    # Step 1: Process R1 reads
    print("Processing R1...")
    pool = Pool(threads)
    chunk_size = 10**5  # Adjust this chunk size for better memory management
    r1_chunks = [fq_r1] * threads
    r1_dict_list = pool.starmap(process_r1, [(fq_r1, barcodes) for fq_r1 in r1_chunks])
    pool.close()
    pool.join()

    # Merge dictionaries from all processes
    r1_dict = {}
    for r in r1_dict_list:
        r1_dict.update(r)

    print(f"Total matching R1 reads: {len(r1_dict)}")

    # Step 2: Process R2 reads
    print("Processing R2...")
    r2_reads = process_r2(fq_r2, r1_dict)

    # Step 3: Write R1 and R2 outputs
    print("Writing output files...")
    write_fq(output_r1, r1_dict.values())
    write_fq(output_r2, r2_reads)

    print("Processing complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split R1 reads based on barcodes and extract corresponding R2 reads.")
    parser.add_argument('--fq_r1', type=str, required=True, help="Input R1 fastq.gz file")
    parser.add_argument('--fq_r2', type=str, required=True, help="Input R2 fastq.gz file")
    parser.add_argument('--barcodes', type=str, required=True, help="File containing barcodes (one per line)")
    parser.add_argument('--output_r1', type=str, required=True, help="Output file for filtered R1 reads")
    parser.add_argument('--output_r2', type=str, required=True, help="Output file for corresponding R2 reads")
    parser.add_argument('--threads', type=int, default=cpu_count(), help="Number of threads to use")
    args = parser.parse_args()

    main(args.fq_r1, args.fq_r2, args.barcodes, args.output_r1, args.output_r2, args.threads)
