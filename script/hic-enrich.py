#!/usr/bin/env python
import pandas as pd
import argparse

# 参数
def process_arguments():
	parser = argparse.ArgumentParser(description='Process BED and matrix files with an optional count.')

    # 添加必需的参数
	parser.add_argument('bed', help='Path to the BED file')
	parser.add_argument('matrix', help='Path to the matrix file')
    # 添加可选的参数
	parser.add_argument('-c','--count', type=int, default=None, help='Number of select genomes')
	return parser.parse_args()

       
        
if __name__ == "__main__":
    args = process_arguments()
    count = 30
    if args.count:
        count = args.count 
    bed = pd.read_csv(args.bed, header = None, sep = "\t")
    matrix = pd.read_csv(args.matrix, header = None, sep = "\t")
    
    id2name = {}
    for index, row in bed.iterrows():
        id2name[row[3]] = row[0]

    ref_read_count = {}
    for i in set(bed[0]):
        ref_read_count[i] = 0

    for index, row in matrix.iterrows():
        i1,i2 = row[0], row[1] 
        if (id2name[i1] == id2name[i2]) and (i1 != i2):
            ref_read_count[id2name[i2]] += row[2]

    max_contig = sorted(ref_read_count, key = lambda x: ref_read_count[x], reverse = True)
    for i in max_contig[:count]:
        print(i)

