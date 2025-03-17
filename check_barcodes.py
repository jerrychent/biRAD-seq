# usage:python check_common_barcodes.py --file1 sample1_barcode.txt --file2 sample2_barcode.txt


import argparse
from collections import defaultdict

def check_common_barcodes(file1, file2):
    # 一次性读取两个文件内容，转换为大写并去除空行
    with open(file1, 'r') as f1:
        barcodes1 = [line.strip().upper() for line in f1 if line.strip()]
        
    with open(file2, 'r') as f2:
        barcodes2 = {line.strip().upper() for line in f2 if line.strip()}  # 使用集合来加快查找速度

    # 统计 file1 中与 file2 中共有的 barcodes，使用字典来保留重复
    common_barcodes = defaultdict(int)
    
    for barcode in barcodes1:
        if barcode in barcodes2:
            common_barcodes[barcode] += 1
    
    print(f"Total barcodes in {file1}: {len(barcodes1)}")
    print(f"Total barcodes in {file2}: {len(barcodes2)}")
    print(f"Common barcodes found: {len(common_barcodes)}")
    
    if common_barcodes:
        print(f"Common barcodes between {file1} and {file2}:")
        for barcode, count in common_barcodes.items():
            print(f"{barcode}: found {count} times in {file1}")
    else:
        print(f"No common barcodes found between {file1} and {file2}.")
    
    return common_barcodes

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check for common barcodes between two barcode files.")
    parser.add_argument('--file1', type=str, required=True, help="First barcode.txt file")
    parser.add_argument('--file2', type=str, required=True, help="Second barcode.txt file")
    args = parser.parse_args()

    check_common_barcodes(args.file1, args.file2)




