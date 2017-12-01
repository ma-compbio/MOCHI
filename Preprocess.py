#encoding:utf-8
from optparse import OptionParser
import pandas as pd 
import numpy as np
import sys
from os import walk
import os

#Map TF-genes to TF-bins
def gene2bin(cell):
	real_cell = cell.split("_")[0]
	for (dirpath, dirnames, filenames) in walk("../Data/%s/used/" %(cell)):
		for filename in filenames:
			if "grn_%s_txStart" %(real_cell) in filename:
				TF_Gene = pd.read_table("../Data/%s/used/%s" %(cell,filename),sep="\t")

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
	TF_Gene = TF_Gene.drop_duplicates(['Source','Bin'])
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

	TF_list = TF_Gene["Source"].unique()
	for TF in TF_list:
		Bin2Node[TF] = count
		Node2Bin[count] = TF
		count += 1
	print "There are %d TF\n" %(len(TF_list))
	
	Bin_list = TF_Gene["Bin"]

	for (dirpath, dirnames, filenames) in walk("../Data/%s/used/" %(cell)):
		for filename in filenames:
			if "hic_%s_KR_10000_inter" %(real_cell) in filename:
				inter_file = pd.read_table("../Data/%s/used/%s" %(cell,filename),sep="\t")

	inter_file = inter_file[inter_file["bin1"] != inter_file["bin2"]]
	#inter_file.to_csv("../Data/%s/remove_inter.txt" %(cell),sep = "\t",index = False)

	Bin_list = np.concatenate((Bin_list,(inter_file["bin1"])))
	Bin_list = np.concatenate((Bin_list,(inter_file["bin2"])))
	
	intra = pd.read_table("../Data/%s/used/Changed_chr.txt" %(cell),sep = "\t")
	intra = intra[intra["bin1"] != intra["bin2"]]

	#intra.to_csv("../Data/%s/remove_chr.txt" %(cell,name),sep = "\t",index = False)
	Bin_list = np.concatenate((Bin_list,(intra["bin1"])))
	Bin_list = np.concatenate((Bin_list,(intra["bin2"])))

	#drop duplicates
	Bin_list = list(set(Bin_list))
	#sort it
	Bin_list.sort()

	if Bin_list[0] == 'bin1':
		print "error!"


	result_file = pd.concat([inter_file,intra])
	result_file.to_csv("../Data/%s/bin_edge.txt" %(cell),sep="\t",index=False)

	chrom_now = (Bin_list[0].split(":"))[0]
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
	real_cell = cell.split("_")[0]
	name_list = range(1,23)
	name_list.append("X")

	f2 = open("../Data/%s/used/Changed_chr.txt" %(cell),"w")
	f2.write("bin1\tbin2\tvalue\n")

	for name in name_list:
		found_flag = False
		for (dirpath, dirnames, filenames) in walk("../Data/%s/used/" %(cell)):
			for filename in filenames:
				if "hic_%s_KR_10000_intra_subset_chr%s_" %(real_cell,name) in filename:
					print filename
					f1 = open("../Data/%s/used/%s" %(cell,filename),"r")
					found_flag = True

		if not found_flag:
			print "chr%s not available for this cell" %(name)
			continue
		lines = f1.readlines()
		for line in lines:
			if line == "bin1\tbin2\tvalue\n":
				continue
			line = line.replace(" ","\t")
			f2.write(line)
		f1.close()

	f2.close()

def GetEdge(cell):
	Bin2Node = Load_dict("../Temp/%s/bin2node.npy" %(cell))
	#N =  max(Bin2Node.values())
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

	file = open("../Data/%s/bin_edge.txt" %(cell),"r")
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
run('gm12878')