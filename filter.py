#!/usr/bin/env python3
"""
FASTA序列过滤脚本
过滤包含过多N和"-"字符的序列
"""

import argparse
import sys
from Bio import SeqIO
import os

def count_n_and_gaps(sequence):
    """
    计算序列中N和"-"的数量
    """
    sequence = str(sequence).upper()
    n_count = sequence.count('N')
    gap_count = sequence.count('-')
    return n_count, gap_count

def filter_sequences(input_file, output_file, max_n_count=None, max_gap_count=None, 
                    max_n_ratio=None, max_gap_ratio=None, min_length=0):
    """
    过滤FASTA序列
    
    参数:
    - input_file: 输入FASTA文件路径
    - output_file: 输出FASTA文件路径
    - max_n_count: 最大N数量
    - max_gap_count: 最大gap数量
    - max_n_ratio: 最大N比例 (0-1)
    - max_gap_ratio: 最大gap比例 (0-1)
    - min_length: 最小序列长度
    """
    
    kept_count = 0
    filtered_count = 0
    
    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            sequence = str(record.seq)
            seq_length = len(sequence)
            
            # 检查最小长度
            if seq_length < min_length:
                filtered_count += 1
                continue
            
            n_count, gap_count = count_n_and_gaps(sequence)
            
            # 检查是否需要过滤
            should_filter = False
            
            # 检查N数量
            if max_n_count is not None and n_count > max_n_count:
                should_filter = True
            
            # 检查gap数量
            if max_gap_count is not None and gap_count > max_gap_count:
                should_filter = True
            
            # 检查N比例
            if max_n_ratio is not None and seq_length > 0:
                n_ratio = n_count / seq_length
                if n_ratio > max_n_ratio:
                    should_filter = True
            
            # 检查gap比例
            if max_gap_ratio is not None and seq_length > 0:
                gap_ratio = gap_count / seq_length
                if gap_ratio > max_gap_ratio:
                    should_filter = True
            
            if should_filter:
                filtered_count += 1
                print(f"过滤序列: {record.id} (长度: {seq_length}, N: {n_count}, gaps: {gap_count})")
            else:
                SeqIO.write(record, out_handle, "fasta")
                kept_count += 1
    
    return kept_count, filtered_count

def main():
    parser = argparse.ArgumentParser(description="过滤FASTA文件中包含过多N和gap的序列")
    
    # 必需参数
    parser.add_argument("input", help="输入FASTA文件路径")
    parser.add_argument("output", help="输出FASTA文件路径")
    
    # 过滤条件参数
    parser.add_argument("--max-n", type=int, help="最大N数量")
    parser.add_argument("--max-gaps", type=int, help="最大gap(-)数量")
    parser.add_argument("--max-n-ratio", type=float, help="最大N比例 (0-1)")
    parser.add_argument("--max-gap-ratio", type=float, help="最大gap比例 (0-1)")
    parser.add_argument("--min-length", type=int, default=0, help="最小序列长度")
    
    # 其他选项
    parser.add_argument("--stats", action="store_true", help="显示统计信息")
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.input):
        print(f"错误: 输入文件 {args.input} 不存在")
        sys.exit(1)
    
    # 检查是否至少指定了一个过滤条件
    if not any([args.max_n, args.max_gaps, args.max_n_ratio, args.max_gap_ratio, args.min_length]):
        print("警告: 未指定任何过滤条件，将复制所有序列")
    
    # 验证比例参数
    if args.max_n_ratio is not None and (args.max_n_ratio < 0 or args.max_n_ratio > 1):
        print("错误: N比例应该在0-1之间")
        sys.exit(1)
    
    if args.max_gap_ratio is not None and (args.max_gap_ratio < 0 or args.max_gap_ratio > 1):
        print("错误: gap比例应该在0-1之间")
        sys.exit(1)
    
    print("开始过滤序列...")
    print(f"输入文件: {args.input}")
    print(f"输出文件: {args.output}")
    
    if args.max_n:
        print(f"最大N数量: {args.max_n}")
    if args.max_gaps:
        print(f"最大gap数量: {args.max_gaps}")
    if args.max_n_ratio:
        print(f"最大N比例: {args.max_n_ratio}")
    if args.max_gap_ratio:
        print(f"最大gap比例: {args.max_gap_ratio}")
    if args.min_length:
        print(f"最小序列长度: {args.min_length}")
    
    try:
        kept_count, filtered_count = filter_sequences(
            args.input, args.output,
            max_n_count=args.max_n,
            max_gap_count=args.max_gaps,
            max_n_ratio=args.max_n_ratio,
            max_gap_ratio=args.max_gap_ratio,
            min_length=args.min_length
        )
        
        print("\n过滤完成!")
        print(f"保留序列数: {kept_count}")
        print(f"过滤序列数: {filtered_count}")
        print(f"总序列数: {kept_count + filtered_count}")
        
        if kept_count + filtered_count > 0:
            print(f"保留比例: {kept_count/(kept_count + filtered_count)*100:.2f}%")
        
    except Exception as e:
        print(f"错误: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()