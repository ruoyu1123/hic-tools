# hic-tools
这是一个三维基因组Hi-C文件处理的工具盒子，多个工具与Hi-Cpro对接,具体功能见本文件下文。  
## INSTALL  
```
git clone http://github.com:ruoyu1123/hic-tools.git  
cd hic-tools
export PATH=$PATH:`pwd`/script/ 
```

## REQUIRES  
python 3.6 - python 3.10  
numpy  
matplotlib  
pandas  

## hic-dishist.py 
输入hic reads的比对文件，计算valid hic reads在参考基因组中的距离分布情况  

## hicplotter.py 
输入[hicpro](https://github.com/nservant/HiC-Pro)生成的bed 和matrix，绘制hic热图，可选择绘制部分染色体（contig）。

## hic-enrich.py  
输入[hicpro](https://github.com/nservant/HiC-Pro)的结果(同上)，统计每个染色体或者contig上的Hi-C信号数量，找到Hi—C富集的参考序列。可以指定输出参考序列的个数。  

## hic-getmap.py  
输入一个Hi-C的比对文件，生成各个contig之间的Hi-C信号，该工具仅仅为粗略统计。

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


## Contact  
- Email: ruoyu1123@outlook.com 
- QQ: 729800244 

