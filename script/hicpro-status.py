#!/usr/bin/env python
import argparse
import os

def get_file_list(path, pattern):
    # 指定目录路径和文件匹配模式
    directory = path  # 将此路径替换为你的目录路径
    file_suffix = pattern  # 可以替换为你想要的文件扩展名或通配
    file_list = []
    for root, _, files in os.walk(directory):
        for filename in files:
            if filename.endswith(file_suffix) and not filename.startswith("."):
                file_list.append(os.path.join(root, filename))
    return file_list

def get_col_name(file_list):
    cols = ["Item"]
    for i in file_list:
        cols.append(i.split("/")[-2])
    return cols

def main():
    # 创建命令行解析器
    parser = argparse.ArgumentParser(description="合并多个文件并处理命令行参数")
    
    # 添加文件名参数，nargs='+' 表示接受一个或多个文件名
    parser.add_argument("path", help="hicpro 结果文件夹")
    parser.add_argument("-o","--outpath", default=".", help="output path[.]")
    
    # 解析命令行参数
    args = parser.parse_args()
    # mpairstat
    # mergestat
    # mRSstat
    patterns = ["mpairstat", "mergestat", "mRSstat"]
    for p in patterns:
        print(p)
        
        files = get_file_list(args.path+"/hic_results", p)
        for i in files:
            print("\t"+i)
        cols = get_col_name(files)
        # 在这里处理你的文件合并逻辑
        merge_files(files, cols, args.outpath+"/"+p)

def merge_files(file_list, cols, output):
    # 创建一个字典来存储第一列和后续列的映射关系
    data_dict = {}
    # 遍历文件列表
    # if output == "mRSstat" or  output=="mpairstat":
    #     for filename in file_list:
    #         with open(filename, "r") as file:
    #             for i in range(5):
    #                 file.readline()
    #             for line in file:
    #                 columns = line.strip().split()
    #                 if len(columns) == 2:
    #                     key, value = columns
    #                     if key in data_dict:
    #                         data_dict[key].append(value)
    #                     else:
    #                         data_dict[key] = [value]
    # else:
    for filename in file_list:
        with open(filename, "r") as file:
             for line in file:
                if line.startswith("#"):
                    continue
                columns = line.strip().split()
                key, value = columns[0], columns[1]
                if key in data_dict:
                    data_dict[key].append(value)
                else:
                    data_dict[key] = [value]
        
    with open(output, "w") as file:
        file.write("\t".join(cols))
        file.write("\n")
        for key, values in data_dict.items():
            file.write(key+"\t"+"\t".join(values))
            file.write("\n")

if __name__ == "__main__":
    main()
