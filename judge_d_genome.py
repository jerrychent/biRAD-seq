import pandas as pd

# 读取数据文件
data = pd.read_csv("reads_counts.txt", sep="\t", header=None)

# 判断第四列的数值，并在第五列输出相应的标识符
data[4] = data[3].apply(lambda x: "D2G" if x > 10 else "DG")

# 保存结果到新文件
data.to_csv("classified_reads_counts.txt", sep="\t", header=False, index=False)

print("文件已处理完成，输出保存为 classified_reads_counts.txt")
