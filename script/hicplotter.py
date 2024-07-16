#! /usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
import argparse
import os

VMAX = 0

def sparse2dense(sparse_matrix_df, n):
    # 创建一个全零矩阵，形状为 (max_row+1, max_col+1)
    dense_matrix = np.zeros((n+1,n+1))
    # 将稀疏矩阵中的值填充到稠密矩阵中的对称位置
    for index, row in sparse_matrix_df.iterrows():
        dense_matrix[row["row"], row["col"]] = row["value"]
        dense_matrix[row["col"], row["row"]] = row["value"]  # 填充对称位置
    return dense_matrix

def select_genome(genome_list, bedfile):
    bed = pd.read_csv(bedfile, sep = "\t", header = None, names=["sequence", "start", "end", "id"])
    genomes = pd.read_csv(genome_list, header = None,sep="\t", names = ["genome","name"])
    select_bed = bed[bed["sequence"].isin(genomes["genome"])]
    select_bed["new_id"] = [i for i in range(select_bed.shape[0])]
    print("select bed row", select_bed.shape)
    print(select_bed)
    return select_bed

def plot(dense_matrix, output, labelindex, label, ticks, vmax):
    print("start draw")
    # 绘制Hi-C热图
    font_title = fm.FontProperties(weight='bold', size=14)
    font_ticks = fm.FontProperties(weight='bold', size=12)
    cmap_red = LinearSegmentedColormap.from_list("", ['white'] + [plt.get_cmap('Reds')(i) for i in range(1, 256)])
    cmap = ListedColormap(['white'] + [plt.get_cmap('Reds')(i) for i in range(1, 256)])
    plt.figure(figsize=(16, 12))  # 设置图形大小
    plt.imshow(dense_matrix, cmap=cmap_red, vmin = 0, vmax = int(vmax),  aspect='auto', interpolation='none')
    for i in range(len(label)):
        plt.text(-3, labelindex[i], label[i], ha='right', va='center', fontproperties = font_ticks)
        plt.text(labelindex[i], dense_matrix.shape[0] , label[i], ha='right', va='top', fontproperties = font_ticks, rotation = 45)
    plt.yticks(ticks,["" for i in range(len(ticks))])
    plt.xticks(ticks,[""]*len(ticks))
    plt.grid()
    plt.colorbar(label='Value')  # 添加颜色条
    plt.title('Hi-C Matrix Heatmap', fontproperties = font_title)  # 设置标题
    plt.tight_layout()
    plt.savefig(output)

def matrix_process(matrix, bed):
    print("start read the matrix")
    # 将spase矩阵筛选出需要的genomes, 然后转换成dense矩阵
    spase_matrix_df = pd.read_csv(matrix, sep = "\t", header = None, names = ["row","col","value"])
    n_row, n_col = bed.shape
    dense_matrix_np = np.zeros((n_row, n_row))
    for index, row in spase_matrix_df.iterrows():
        if (row["row"] in bed["id"]) and (row["col"] in bed["id"]):
            x = bed.loc[bed["id"] == row[0], "new_id"]
            y = bed.loc[bed["id"] == row[1], "new_id"]
            dense_matrix_np[x, y] = row["value"]
            dense_matrix_np[y, x] = row["value"]
    print(dense_matrix_np)
    return dense_matrix_np

def label_index(bed):
    print("start obtain label")
    tmp = bed[["sequence","new_id"]]
    result = tmp.groupby("sequence").agg(['first', 'last']).reset_index()
    result.columns=["group", "start", "end"]
    return result


def partial_draw(bed, matrix, genomes, output):
    # 如果提供了需要绘制hic热图的genomes，则进入此函数进行绘制
    select_bed = select_genome(genomes, bed)
    dense_matrix = matrix_process(matrix, select_bed)
    genomes = pd.read_csv(genome_list,sep = "\t", header = None, names = ["genome","name"])
    if VMAX == 0:
        max_color_value = float(np.max(dense_matrix.diagonal())/3)
    else:
        max_color_value = VMAX
    print("vmax:", max_color_value)
    label_raw = label_index(select_bed)
    ticks = [i-1 for i in label_raw["start"]]
    label_indexs = (label_raw["end"] - label_raw["start"])/2 + label_raw["start"]-1
    tmpdf = pd.merge(label_raw, genomes, left_on="group", right_on="genome")
    print(tmpdf)
    label = tmpdf["name"]
    plot(dense_matrix, output, label_indexs, label,ticks, max_color_value)
    return 0

def all_draw(bed_file, matrix, output):
    sparse_matrix = pd.read_csv(matrix, sep="\t", header = None, names=["row", "col", "value"])
    bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=["sequence", "start", "end", "id"])
    n = bed_df.shape[0]
    dense_matrix = sparse2dense(sparse_matrix,n)
    if VMAX == 0:
        max_color_value = float(np.max(dense_matrix.diagonal())/3)
    else:
        max_color_value = VMAX
    print("start read label")
    label_df = bed_df.groupby('sequence').agg(['first', 'last']).reset_index()
    label_df.columns = ['Group', "1", '2',"3","4","5","6"]
    ticks = list(label_df["5"])
    label_index = (label_df["6"]-label_df["5"])/2 + label_df["5"] - 1
    label = list(label_df["Group"])
    plot(dense_matrix, output, label_index, label,ticks, max_color_value)
    return 0
    
    
    
    
if __name__ == "__main__":
    
    bedfile = ""
    hicpro_matrix = ""
    
    # Command parameters parse
    parser = argparse.ArgumentParser(description="Hi-C heatmap plotter")
    parser.add_argument('bed', type=str, help='The bed file corresponding matrix')
    parser.add_argument('matrix', type=str, help='HiCpro result matrix')
    parser.add_argument('-g','--genomes', type=str, help='Chromsome or contig or genomes list file need to plot. Default: [all]')
    parser.add_argument('--vmax', type = int, help = "Max value of HiC heatmap")
    parser.add_argument('-o', '--outpath', type=str, help='Output path. Default: ./')
    args = parser.parse_args()
    bedfile = args.bed
    hicpro_matrix = args.matrix
    genome_list = args.genomes
    if args.vmax:
        VMAX = args.vmax
    if args.outpath:
        output = args.outpath
    else:
        output = "./hicmap.pdf"
    if genome_list:
        partial_draw(bedfile, hicpro_matrix, genome_list, output)
    else:
        all_draw(bedfile, hicpro_matrix, output)
    print("Hi-C heatmap was drawn")
        
    
