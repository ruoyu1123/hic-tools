#!/usr/bin/env python
import networkx as nx
import pysam

class hic:
    def __init__(self, bam_path):
        self.bam_path=bam_path
        self.index = {}
        self.reference = []
        self.Enzyme_count = {}
        self.raw_map={}
        self.map_NF=nx.Graph()
    def create_map(self):
        with pysam.AlignmentFile(self.bam_path, 'rb') as bam:
            self.ref_length = bam.lengths
            self.reference = bam.references
            c = 0
            for i in bam.references:
                self.index[i] = c
                c+=1
            for read1 in bam:
                if read1.is_duplicate:
                    continue
                if read1.is_unmapped or read1.mate_is_unmapped:
                    continue
                ref1,ref2=read1.reference_name,read1.next_reference_name
                
                if ref1==ref2:
                    continue
                # if read1.is_proper_pair == False:
                #     continue
                # if read1.is_secondary:
                #     continue
                # read2 = bam.mate(read1)
                # print(read1.get_tags("MQ"))
                if read1.mapping_quality <10: #and read1.query_alignment_length<130:
                    continue
#               将两条contig的重叠部分分到前面，只取比对到第二条contig重叠部分后面的比对作为权重
                '''
                if (r1,r2) in link:
                    
                    link_info=link[(r1,r2)]
                    if link_info[1]=='-':
                        # print(link_info[1])
                        if read.mate_is_reverse and int(read.next_reference_start) > link_info[2]:
                            # print(link[(r1,r2)][3])
                            link[(r1,r2)][3]+=1
                        if read.mate_is_reverse==False and int(read.next_reference_start) < self.ref_length[r2]-link_info[2]:
                            link[(r1,r2)][3]+=1
                    if link_info[1]=='+':
                        if read.mate_is_reverse and int(read.next_reference_start) < self.ref_length[r2]-link_info[2]:
                            link[(r1,r2)][3]+=1
                        if read.mate_is_reverse==False and int(read.next_reference_start) > link_info[2]:
                            link[(r1,r2)][3]+=1
                # if (r2,r1) in link:
                #     link_info=link[(r2,r1)]
                #     if link_info[1]=='-':
                #         if read.is_reverse and read.reference_start > link_info[2]:
                #             link[(r2,r1)][3]+=1
                #         if read.is_reverse==False and read.reference_start < self.ref_length[r2]-link[(r1,r2)][2]:
                #             link[(r2,r1)][3]+=1
                #     if link_info[1]=='+':
                #         if read.is_reverse and read.reference_start < self.ref_length[r2]-link_info[2]:
                #             link[(r2,r1)][3]+=1
                #         if read.is_reverse==False and read.reference_start > link_info[2]:
                #             link[(r2,r1)][3]+=1
                '''

                # 注释掉else直接建立全局的hic信号
                # else:
                if (ref1, ref2) not in self.raw_map and (ref2, ref1) not in self.raw_map:
                    self.raw_map[(ref1, ref2)] = 1
                    continue
                elif (ref1,ref2) in self.raw_map:
                    self.raw_map[(ref1,ref2)]+=1
                    continue
                elif (ref2, ref1) in self.raw_map:
                    self.raw_map[(ref2, ref1)]+=1
                # print(read.reference_name,read.next_reference_name)
                # print(read.reference_length,read.rlen)
                # break #break zaihe
        # self.raw_map = {k:v for k, v in self.raw_map.items() if v != 1}
            # tmp.append(read)
    def normalized(self, ntype = "length"):
        if ntype == "length":
            for edge in self.raw_map:
                u, v = edge[0], edge[1]
                weight = self.raw_map[(u, v)]
                if weight != 1:
                    self.map_NF.add_edge(u, v, weight=weight/(
                        self.ref_length[u]+self.ref_length[v]))
        else:
            for edge in self.raw_map:
                u, v = edge[0], edge[1]
                self.map_NF.add_edge(u, v, weight=self.raw_map[(u,v)](
                    self.ref_length[u]+self.ref_length[v]))
                
if __name__  == "__main__":
    import sys
    import os
    if len(sys.argv) < 2:
        print("Usage:\n\thic-getmap.py bam output")
        sys.exit(0)
    shic = hic(sys.argv[1])
    shic.create_map()
    #print(shic.raw_map)
    with open(sys.argv[2], "w") as file:
        for e, v in shic.raw_map.items():
            file.write("{}\t{}\t{}\n".format(e[0], e[1], str(v)))
    print("Writing finshed!")
