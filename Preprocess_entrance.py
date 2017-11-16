#encoding:utf-8
from optparse import OptionParser
import pandas as pd 
import numpy as np
import sys
import os

#Map TF-genes to TF-bins
def gene2bin(cell):
	real_cell = cell.split("_")[0]
	TF_Gene = pd.read_table("../Data/%s/used/grn_%s_txStart.txt" %(cell,real_cell),sep = "\t")
	Gene2Bin = pd.read_table("../Data/%s/used/gene_chrom_bin_num_hg19_%s.txt" %(cell,real_cell),sep = "\t")
	Target = TF_Gene["Target"]
	Gene2Bin.index= Gene2Bin["TF"]

	selected = Gene2Bin.reindex(Target)
	selected.dropna()
	if len(selected) == len(Target):
		print ("Successfully match")
	else:
		print ("Doesn't match")

	selected.index = xrange(len(selected))


	chrom = pd.Series([item.split(':')[0] for item in selected["Bin_10kb"]])
	bin = pd.Series([item.split(':')[1] for item in selected["Bin_10kb"]])


	Bin_10kb = selected["Bin_10kb"]
	TF_Gene["Bin"] = Bin_10kb
	TF_Gene["chrom"] = chrom
	TF_Gene["bin"] = bin
	print len(TF_Gene)
	TF_Gene = TF_Gene.drop_duplicates(['Source','Bin'])
	print len(TF_Gene)
	TF_Gene.to_csv("../Data/%s/TF_Gene_mapped.txt" %(cell),sep = "\t",index = False)

def Load_dict(path):
	return np.load(path).item()

def Build_dict(cell):
	real_cell = cell.split("_")[0]
	#dict from bin to gene
	Bin2Gene = {}
	table = pd.read_table("../Data/%s/used/gene_chrom_bin_num_hg19_%s.txt" %(cell,real_cell),sep = "\t")
	for i in xrange(len(table)):
		#more than one gene in one bin
		if Bin2Gene.has_key(table["Bin_10kb"][i]):
			Bin2Gene[table["Bin_10kb"][i]] = Bin2Gene[table["Bin_10kb"][i]]+ ";" + table["TF"][i]
		else:
			Bin2Gene[table["Bin_10kb"][i]] = table["TF"][i]
	np.save("../Temp/%s/bin2gene.npy" %(cell),Bin2Gene)


	Node2Bin = {}
	Bin2Node = {} 
	count = 0
	TF_Gene = pd.read_table("../Data/%s/TF_Gene_mapped.txt" %(cell), sep = "\t")
	#TF仅使用出现在mapped后的
	TF_list = TF_Gene["Source"].unique()
	for TF in TF_list:
		Bin2Node[TF] = count
		Node2Bin[count] = TF
		count += 1
	print "There are %d TF\n" %(len(TF_list))
	
	#Bin按照出现在TF-BIN，然后INTER，然后INTRA中的合集
	Bin_list = TF_Gene["Bin"]
	inter_file = pd.read_table("../Data/%s/used/hic_%s_KR_10000_inter_top_1_genes.txt" %(cell,real_cell),sep = "\t")
	inter_file = inter_file[inter_file["bin1"] != inter_file["bin2"]]
	inter_file.to_csv("../Data/%s/remove_inter.txt" %(cell),sep = "\t",index = False)
	Bin_list = np.concatenate((Bin_list,(inter_file["bin1"])))
	Bin_list = np.concatenate((Bin_list,(inter_file["bin2"])))
	
	
	name_list = range(1,23)
	name_list.append("X")
	
	if cell == "k562" or cell == "k562_2":
		name_list.remove(9)

	for name in name_list:
		file = pd.read_table("../Data/%s/used/chr%s.txt" %(cell,name),sep = "\t")
		#drop duplicates
		file = file[file["bin1"] != file["bin2"]]
		file.to_csv("../Data/%s/remove_chr%s.txt" %(cell,name),sep = "\t",index = False)
		Bin_list = np.concatenate((Bin_list,(file["bin1"])))
		Bin_list = np.concatenate((Bin_list,(file["bin2"])))

	#drop duplicates
	Bin_list = list(set(Bin_list))
	#sort it
	Bin_list.sort()
	if Bin_list[0] == 'bin1':
		Bin_list = Bin_list[2:]
	print Bin_list[0:10]
	chrom_now = (Bin_list[0].split(":"))[0]
	print len(Bin_list)
	index_file = open("../Temp/%s/chrom_bin_index.txt" %(cell),"w")
	index_file.write("chrom\tindex_start\tindex_end\n")
	index_file.write("%s\t%d\t" %(chrom_now,590))
	for Bin in Bin_list:
		chrom = Bin.split(":")[0]
		if chrom_now!= chrom:
			chrom_now = chrom
			index_file.write("%d\n%s\t%d\t" %(count,chrom,count))
		Bin2Node[Bin] = count
		Node2Bin[count] = Bin
		count += 1
	index_file.write("%d\n" %(count))

	np.save("../Temp/%s/bin2node.npy" %(cell),Bin2Node)

	np.save("../Temp/%s/node2bin.npy" %(cell),Node2Bin)

def Change(cell):

	name_list = range(1,23)
	name_list.append("X")

	if cell == "k562" or cell == "k562_2":
		name_list.remove(9)

	for name in name_list:
		strs = cell.split("_")
		if len(strs) == 2:
			oe_num = strs[1]
		else:
			oe_num = 1
		real_cell = strs[0]
		try:
			f1 = open("../Data/%s/used/hic_%s_KR_10000_intra_subset_chr%s_oe_%s_genes.txt" %(cell,real_cell,name,str(oe_num)),"r")
		except:
			f1 = open("../Data/%s/used/hic_%s_KR_10000_intra_subset_chr%s%s_genes.txt" %(cell,real_cell,name,cell[len((cell.split("_")[0])):]),"r")

		f2 = open("../Data/%s/used/chr%s.txt" %(cell,name),"w")

		lines = f1.readlines()
		if lines[0] != "bin1\tbin2\tvalue\n":
			f2.write("bin1\tbin2\tvalue\n")
		for line in lines:
			line = line.replace(" ","\t")
			f2.write(line)
		f1.close()
		f2.close()

def GetEdge(cell):
	Bin2Node = Load_dict("../Temp/%s/bin2node.npy" %(cell))
	#Bin2Node = Load_dict("../Temp/global_bin2node2.npy")
	N =  max(Bin2Node.values())
	#HiC_Value = np.zeros((N+1,N+1),dtype = "float32")

	f = open("../Data/%s/Edge_list.txt" %(cell),"w")
	f.write("x\ty\n")
	file = open("../Data/%s/TF_Gene_mapped.txt" %(cell),"r")
	lines = file.readlines()[1:]
	for line in lines:
		a = line.split("\t")
		node1 = Bin2Node[a[0]]
		if Bin2Node.has_key(a[2]) == False:
			continue
		node2 = Bin2Node[a[2]]
		f.write("%d\t%d\n" %(node1,node2))

	name_list = range(1,23)
	name_list.append("X")

	if cell == "k562" or cell == "k562_2":
		name_list.remove(9)

	for name in name_list:
		file = open("../Data/%s/remove_chr%s.txt" %(cell,name),"r")
		lines = file.readlines()[1:]
		for line in lines:
			a = line.split("\t")
			node1 = Bin2Node[a[0]]
			node2 = Bin2Node[a[1]]
			f.write("%d\t%d\n" %(node1,node2))
			#HiC_Value[node1,node2] = a[2]


	file = open("../Data/%s/remove_inter.txt" %(cell),"r")

	lines = file.readlines()[1:]
	for line in lines:
		a = line.split("\t")
		node1 = Bin2Node[a[0]]
		node2 = Bin2Node[a[1]]
		f.write("%d\t%d\n" %(node1,node2))
		#HiC_Value[node1,node2] = a[2]
	f.close()

	a = pd.read_table("../Data/%s/Edge_list.txt" %(cell),sep = "\t")
	a = a.drop_duplicates()
	a.to_csv("../Data/%s/Edge_list.txt" %(cell),sep ="\t",index = False)
	#np.save("../Temp/%s/HiC_Value.npy" %(cell),HiC_Value)

def run(cell):
	Change(cell)
	gene2bin(cell)
	Build_dict(cell)
	GetEdge(cell)