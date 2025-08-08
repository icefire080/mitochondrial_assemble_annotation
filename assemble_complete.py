from Bio import SeqIO
from Bio.Seq import reverse_complement
import argparse
import subprocess
import os

# 设置命令行参数
parser = argparse.ArgumentParser(description="Assemble mitochondrial genome from contigs.")
parser.add_argument("-r", "--reference", required=True, help="Reference sequence file (FASTA format).")
parser.add_argument("-c", "--contigs", required=True, help="Contigs file (FASTA format).")
parser.add_argument("-o", "--output", required=True, help="Output file for the assembled genome (FASTA format).")
args = parser.parse_args()

# 输入文件
reference_file = args.reference
contigs_file = args.contigs
output_file = args.output

#从输出文件名中提取序列名称（去掉文件扩展名）
sequence_name = os.path.splitext(os.path.basename(output_file))[0]

# 临时文件
delta_file = "output.delta"
coords_file = "output.coords"

# 运行 nucmer 进行比对
print("Running nucmer to align contigs to the reference...")
nucmer_command = f"nucmer --maxmatch -p output {reference_file} {contigs_file}"
subprocess.run(nucmer_command, shell=True, check=True)

# 运行 show-coords 提取比对坐标
print("Running show-coords to extract alignment coordinates...")
show_coords_command = f"show-coords -r -c -l -T {delta_file} > {coords_file}"
subprocess.run(show_coords_command, shell=True, check=True)

# 读取参考序列
reference = next(SeqIO.parse(reference_file, "fasta")).seq
reference_length = len(reference)

# 读取 contigs
contigs = {rec.id: rec.seq for rec in SeqIO.parse(contigs_file, "fasta")}

# 解析比对坐标
with open(coords_file, "r") as f:
    lines = f.readlines()[5:]  # 跳过表头
    coords = []
    for line in lines:
        parts = line.strip().split()
        ref_start = int(parts[0])  # 参考序列起始位置
        ref_end = int(parts[1])    # 参考序列结束位置
        query_start = int(parts[2])  # 查询序列起始位置
        query_end = int(parts[3])    # 查询序列结束位置
        identity = float(parts[6])   #比对相似性百分比
        contig_id = parts[-1]        # contig ID

        # 判断比对方向
        if query_start < query_end:
            strand = "+"  # 正义链
        else:
            strand = "-"  # 负义链

        # 如果比对相似性低于90%，跳过该contig
        if identity < 95.0:
            print(f"Contig {contig_id} skipped due to low identity ({identity}%).")
            continue
        
        coords.append((ref_start, ref_end, contig_id, strand))

# 按参考序列位置排序
coords.sort()

# 处理环形参考序列
# 如果 contig 跨越起点和终点，将其拆分为两部分
new_coords = []
for start, end, contig_id, strand in coords:
    if start > end:  # 跨越起点和终点
        new_coords.append((start, reference_length, contig_id, strand))
        new_coords.append((1, end, contig_id, strand))
    else:
        new_coords.append((start, end, contig_id, strand))
coords = sorted(new_coords)

# 构建完整基因组
full_genome = ["N"] * reference_length  # 初始化全基因组为 N
for start, end, contig_id, strand in coords:
    contig_seq = str(contigs[contig_id])
    contig_length = len(contig_seq)

    # 如果是负义链，将 contig 序列反向互补
    if strand == "-":
        contig_seq = str(reverse_complement(contigs[contig_id]))

    # 处理长度超出的 contigs
    if end - start + 1 > reference_length:  # contig 长度超出参考序列
        print(f"Contig {contig_id} exceeds reference length. Trimming overlap...")
        # 计算超出部分
        overlap_length = (end - start + 1) - reference_length
        # 截断 contig 序列，去除首尾重叠区域
        contig_seq = contig_seq[overlap_length // 2 : contig_length - overlap_length // 2]
        contig_length = len(contig_seq)

    # 检查是否有重叠
    overlap_length = 0
    for pos in range(start - 1, end):  # 转换为 0-based 索引
        if full_genome[pos % reference_length] != "N":
            overlap_length += 1

    # 如果有重叠，选择质量较高的序列（这里假设当前 contig 质量更高）
    if overlap_length > 0:
        print(f"Overlap detected between contigs at reference position {start}-{end}")

    # 将 contig 序列填充到全基因组中
    for i in range(contig_length):
        ref_pos = (start - 1 + i) % reference_length  # 处理环形基因组
        if full_genome[ref_pos] == "N" or full_genome[ref_pos] == contig_seq[i]:
            full_genome[ref_pos] = contig_seq[i]
        else:
            # 如果有冲突，保留当前 contig 的碱基（或根据需要选择其他策略）
            full_genome[ref_pos] = contig_seq[i]

# 输出完整基因组
full_genome = "".join(full_genome)
with open(output_file, "w") as f:
    f.write(f">{sequence_name}\n{full_genome}\n")

print(f"Full genome assembled and saved to {output_file}")