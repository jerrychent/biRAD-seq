# usage: python remove_common_barcodes.py --file1 barcode1.txt --file2 barcode2.txt --output1 unique_barcode1.txt --output2 unique_barcode2.txt
import argparse

def remove_common_barcodes(file1, file2, output1, output2):
    # 读取文件内容并去掉空行，全部转换为大写
    with open(file1, 'r') as f1:
        barcodes1 = {line.strip().upper() for line in f1 if line.strip()}
        
    with open(file2, 'r') as f2:
        barcodes2 = {line.strip().upper() for line in f2 if line.strip()}

    # 找到两个文件中的公共 barcodes
    common_barcodes = barcodes1 & barcodes2
    
    print(f"Total barcodes in {file1}: {len(barcodes1)}")
    print(f"Total barcodes in {file2}: {len(barcodes2)}")
    print(f"Common barcodes: {len(common_barcodes)}")

    # 从各自的集合中删除公共的 barcodes
    unique_barcodes1 = barcodes1 - common_barcodes
    unique_barcodes2 = barcodes2 - common_barcodes

    # 输出结果到新的文件
    with open(output1, 'w') as out1:
        for barcode in sorted(unique_barcodes1):  # 排序后写入文件
            out1.write(barcode + '\n')

    with open(output2, 'w') as out2:
        for barcode in sorted(unique_barcodes2):  # 排序后写入文件
            out2.write(barcode + '\n')

    print(f"Unique barcodes saved to {output1} and {output2}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove common barcodes between two barcode files.")
    parser.add_argument('--file1', type=str, required=True, help="First barcode.txt file")
    parser.add_argument('--file2', type=str, required=True, help="Second barcode.txt file")
    parser.add_argument('--output1', type=str, required=True, help="Output file for unique barcodes from file1")
    parser.add_argument('--output2', type=str, required=True, help="Output file for unique barcodes from file2")
    
    args = parser.parse_args()

    remove_common_barcodes(args.file1, args.file2, args.output1, args.output2)