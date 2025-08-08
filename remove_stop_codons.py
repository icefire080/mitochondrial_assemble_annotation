#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
脚本功能：去除FASTA文件中每个序列的终止密码子(TAA/TAG/TGA)
支持多行FASTA格式(每行80个碱基的标准格式)
"""

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def remove_stop_codons(input_file, output_file, line_width=80):
    """
    去除FASTA序列中的终止密码子并保存
    
    参数:
        input_file: 输入FASTA文件路径
        output_file: 输出FASTA文件路径
        line_width: 输出文件中每行的碱基数(默认80)
    """
    try:
        # 读取输入文件
        records = list(SeqIO.parse(input_file, "fasta"))
        processed_records = []
        
        for record in records:
            # 去除序列中的空格和换行符，合并为连续序列
            seq_str = str(record.seq).replace("\n", "").replace(" ", "")
            seq = Seq(seq_str)
            
            # 检查并去除终止密码子
            if len(seq) >= 3:
                last_codon = str(seq[-3:]).upper()
                if last_codon in ["TAA", "TAG", "AGA", "AGG"]:
                    seq = seq[:-3]
            
            # 创建新的SeqRecord对象，保留原始ID和描述
            new_record = SeqRecord(
                seq,
                id=record.id,
                description=record.description
            )
            processed_records.append(new_record)
        
        # 写入输出文件，保持每行line_width个碱基
        with open(output_file, "w") as f:
            for record in processed_records:
                f.write(f">{record.id} {record.description}\n")
                for i in range(0, len(record.seq), line_width):
                    line = str(record.seq[i:i+line_width])
                    f.write(f"{line}\n")
        
        print(f"成功处理并保存到 {output_file}")
        print(f"共处理 {len(processed_records)} 条序列")
        
    except Exception as e:
        print(f"处理过程中出错: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    # 命令行使用说明
    if len(sys.argv) != 3:
        print("用法: python remove_stop_codons.py 输入文件.fasta 输出文件.fasta")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    print(f"开始处理文件: {input_file}")
    remove_stop_codons(input_file, output_file)