#!/bin/env python 

import pysam
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse

def calculate_insert_distances(bam_file_path):
    """
    计算并统计BAM文件中reads对的insert距离。
    
    参数:
    bam_file_path (str): BAM文件路径。
    
    返回:
    list: 包含所有有效reads对insert距离的列表。
    """
    insert_sizes = []
    interchromosomal_count = 0
    total_pair = 0
    # 打开BAM文件
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        for read in bam_file.fetch(until_eof=True):
            # 确保read是配对的并且两端都已成功比对
            if not read.is_paired or read.mate_is_unmapped or read.is_secondary:
                continue
            if read.is_read1:
                total_pair += 1
                if read.reference_name != read.next_reference_name:
                    interchromosomal_count += 1
                    continue
                # 计算insert size
                insert_size = abs(read.next_reference_start - (read.reference_start + read.query_length))  # 注意加read.query_length
                # print(read.next_reference_start,read.reference_start)
                # 注意insert size的正负，根据实际情况调整逻辑
                if insert_size > 500:
                    insert_sizes.append(insert_size)
                
                    
    
    return insert_sizes,total_pair, interchromosomal_count

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process BAM file and output result.")
    parser.add_argument("bam", help="Path to the BAM file.", type=str)
    parser.add_argument("-o", "--output", help="Output file name.The format is pdf", type=str)
    args = parser.parse_args()
    if args.output:
        output = args.output
    else:
        output = "./hic_distance_histplot.pdf"
    bam_path = args.bam
    insert_sizes, total_pair, interchr = calculate_insert_distances(bam_path)
    print("Total pairs: ", total_pair)  
    print("Inter chrom: ", interchr)
    print("500-1000: ", len([i for i in insert_sizes if (i<1000)]))
    print("1k-10k: ", len([i for i in insert_sizes if (i>=1000) and (i<10000)]))
    print("10k-20k: ", len([i for i in insert_sizes if (i>=10000) and (i<20000)]))
    print(">20k:", len([i for i in insert_sizes if i>=20000]))
    #plot distance histprint(">20k:", len(insert_sizes>20000))
    bin_edges = np.arange(0, 2000000, 1000)
    plt.figure(figsize=(6, 4))
    plt.hist(insert_sizes, bins=bin_edges, edgecolor='black', alpha=0.7)
    plt.xlabel('Insertion Distance (bp)', fontweight = "bold")
    plt.ylabel('Frequency', fontweight = "bold")
    
    plt.savefig(output, format='pdf')
