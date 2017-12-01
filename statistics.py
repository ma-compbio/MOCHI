#encoding:utf-8
import pandas as pd
import numpy as np
from Build_Result import Build_Result
import os
import shutil
from scipy import stats
import random

global_TF_num = 590

def getlist(path):
	file = open(path,"r")
	lines = file.readlines()
	list = []
	for line in lines:
		list.append(line.strip())
	print list
	return list

def Load_dict(path):
	return np.load(path).item()

def Motifcount_npdd(connect,TF_num):
	#print "Start Counting"
	split_point = TF_num
	#Only focus on pro - dna interaction
	pro_connect = np.copy(connect)
	pro_connect[split_point:,split_point:] = 0
	#Only dna-dna interaction
	dna_connect = connect - pro_connect
	
	#the matrix is the num of pro that a pair of dna share,also the num of dna that a pair of pro share
	share = pro_connect.dot(pro_connect)
	tria_count_dna= share * dna_connect

	weight_count_dna =  tria_count_dna * (tria_count_dna - 1)/2 
	weight_count_dna = np.where(weight_count_dna < 0, 0, weight_count_dna)

	return (np.sum(weight_count_dna) * 0.5)

def Motifcount_pdd(connect,TF_num):
	split_point = TF_num
	connect = self.connect

	#Only focus on pro - dna interaction
	pro_connect = np.copy(connect)
	pro_connect[split_point:,split_point:] = 0
	#Only dna-dna interaction
	dna_connect = connect - pro_connect

	dna_weight = pro_connect.dot(pro_connect) * dna_connect
	
	return (np.sum(dna_weight) * 0.5)

def Add_Leftover(cell,name = ""):
	clusters = np.load("../Result/%s/clusters%s.npy" %(cell,name))
	node2bin = Load_dict("../Temp/%s/node2bin.npy" %(cell))
	gene_num_all = max(node2bin.keys())

	all_gene = range(global_TF_num,gene_num_all+1)
	gene_used = []
	for i in xrange(len(clusters)):
		cluster = clusters[i]
		bins = cluster[cluster>=global_TF_num]
		for Bin in bins:
			all_gene.remove(Bin)
			gene_used.append(Bin)

	left = np.unique(np.concatenate((np.arange(global_TF_num),all_gene),axis = -1))
	used = np.unique(np.concatenate((np.arange(global_TF_num),gene_used),axis = -1))
	clusters = list(clusters)
	clusters.append(left)
	clusters.append(used)
	clusters = np.asarray(clusters)

	np.save("../Result/%s/clusters%s_add_left.npy" %(cell,name),clusters)
	return clusters

def Split_chromosome(cell,name = ""):
	node2bin = Load_dict("../Temp/%s/node2bin.npy" %(cell))
	bin2gene = Load_dict("../Temp/%s/bin2gene.npy" %(cell))
	connect = np.load("../Temp/%s/global_connect.npy" %(cell))
	clusters = np.load("../Result/%s/clusters%s.npy" %(cell,name))
	interval_table = pd.read_table("../Temp/%s/chrom_bin_index.txt" %(cell),sep = "\t")
	new_result = []
	cluster_index = []

	for j in xrange(len(clusters)):
		cluster = clusters[j]
		cluster.sort()
		TFs = cluster[cluster<global_TF_num]
		bins = cluster[cluster>=global_TF_num]



		for i in xrange(len(interval_table)):
			start = interval_table['index_start'][i]
			end = interval_table['index_end'][i]

			index_we_care = ((bins < end) & (bins >= start))
			if np.sum(index_we_care) != 0:
				result = np.concatenate((TFs,bins[index_we_care]),axis = -1)
				local_connect = connect[result,:]
				local_connect = local_connect[:,result]
				a  = np.sum(local_connect,axis = 0)
				b = (a > 0)
				result = result[b]
				new_result.append(result)
				cluster_index.append(j)

	new_result = np.asarray(new_result)
	cluster_index = np.asarray(cluster_index)
	#print len(new_result)
	np.save("../Result/%s/clusters%s_split.npy" %(cell,name),new_result)
	np.save("../Result/%s/clusters%s_split_index.npy" %(cell,name),cluster_index)

def summary_for_cluster(cell,name = ""):

	clusters = np.load("../Result/%s/clusters%s.npy" %(cell,name))
	node2bin = Load_dict("../Temp/%s/node2bin.npy" %(cell))
	gene_num_all = max(node2bin.keys())
	bin2gene = Load_dict("../Temp/%s/bin2gene.npy" %(cell))
	clusters_self = np.load("../Result/%s/clusters_noparent%s.npy" %(cell,name))
	PPI = np.load("../Temp/PPI.npy")
	chrom_num_list = []

	complex_file = open("../Temp/PPI_complex.txt")
	lines = complex_file.readlines()
	complex_list = []
	for line in lines:
		complex_list.append(set(line.split("\t")))

	TF2index = Load_dict("../Temp/TF2index.npy")
	index2TF = Load_dict("../Temp/index2TF.npy")

	TF_num = []
	TF_num_self = []
	TF_num_inherit = []
	PPI_num_list = []
	PPI_density_list = []
	PPI_detail_list = []

	ppi_complex = []
	ppi_complex_detail = []
	ppi_complex_all = []

	random_density_list = []
	source = []




	for i in xrange(len(clusters)):
		source.append("him")
		cluster = clusters[i]
		cluster.sort()
		cluster_noparent = clusters_self[i]
		cluster_noparent.sort()

		bins = cluster[cluster>=global_TF_num]
		TF = cluster[cluster<global_TF_num]
		TF_self = cluster_noparent[cluster_noparent<global_TF_num]

		TF_for_PPI = []
		for tf in TF:
			TF_for_PPI.append(TF2index[node2bin[tf]])

		TF_set = set(index2TF[tf] for tf in TF_for_PPI)
		flag = False
		detail = ""
		all1=""
		for pro_complex in complex_list:
			intersect = list(pro_complex&TF_set)
			if len(intersect) >= 2:
				flag = True
				detail+=";".join(intersect)+"|"
				all1+=";".join(list(pro_complex))+"|"

		ppi_complex.append(flag)
		ppi_complex_detail.append(detail)
		ppi_complex_all.append(all1)
				

		local_ppi = np.copy(PPI[TF_for_PPI,:])
		local_ppi = local_ppi[:,TF_for_PPI]

		local_ppi_tril = np.tril(local_ppi,k = -1) 
		coor_x = np.where(local_ppi_tril > 0)[0] 
		coor_y = np.where(local_ppi_tril > 0)[1]

		temp_list = []
		for i in xrange(len(coor_x)):
			x = coor_x[i]
			y = coor_y[i]
			temp_list.append(index2TF[TF_for_PPI[x]]+"|"+index2TF[TF_for_PPI[y]])
		PPI_detail_list.append(";".join(temp_list))



		ppi_num = np.sum(local_ppi) * 0.5
		ppi_density = ppi_num / (len(TF) *(len(TF) - 1) * 0.5)

		PPI_num_list.append(ppi_num)
		PPI_density_list.append(ppi_density)

		random_TF = np.random.randint(0,global_TF_num,size = len(TF))
		local_ppi = np.copy(PPI[random_TF,:])
		local_ppi = local_ppi[:,random_TF]

		ppi_num = np.sum(local_ppi) * 0.5
		ppi_density = ppi_num / (len(random_TF) *(len(random_TF) - 1) * 0.5)
		random_density_list.append(ppi_density)

		chrom_now = node2bin[bins[0]].split(":")[0]
		chrom_num = 1

		for j in xrange(1,len(bins)):
			chrom = node2bin[bins[j]].split(":")[0]
			if chrom != chrom_now:
				chrom_num += 1
				chrom_now = chrom

		TF_num_self.append(len(TF_self))
		TF_num_inherit.append(len(TF) - len(TF_self))
		chrom_num_list.append(chrom_num)
		TF_num.append(len(TF))

	index = np.arange(0,len(clusters))
	ratio1= np.array(TF_num_inherit)*1.0/np.array(TF_num_self)
	ratio2 = np.array(TF_num_inherit)*1.0/np.array(TF_num)
	source = source[:-2]
	source.append("merged non-him")
	source.append("merged him")
	cluster_summary = {'index':index,'chrom_num':chrom_num_list,'TF_number':TF_num,\
						'TF_num_inherit':TF_num_inherit,'TF_self':TF_num_self,\
						'inherit/self':ratio1,'inherit/all':ratio2,\
						'ppi_num':PPI_num_list,'ppi_density':PPI_density_list,\
						'ppi_detail':PPI_detail_list,'source':source,'ppi_complex':ppi_complex,\
						'ppi_complex_detail':ppi_complex_detail,'ppi_complex_all':ppi_complex_all}
	cluster_summary = pd.DataFrame(cluster_summary,\
						columns= ['index','source','chrom_num','TF_number','TF_num_inherit',\
						'TF_self','inherit/self','inherit/all','ppi_num','ppi_density',\
						'ppi_detail','ppi_complex','ppi_complex_detail','ppi_complex_all'])
	cluster_summary.to_csv("../Statistics_Result/statistics/%s_cluster_summary.txt" %(cell),sep = "\t",index = False)
	
	#print random_density_list
	print stats.wilcoxon(PPI_density_list[0:len(clusters) - 2], random_density_list[0:len(clusters) - 2])
	print np.mean(cluster_summary['ppi_density'][0:len(cluster_summary) - 2]), cluster_summary['ppi_density'][len(cluster_summary)-1]
	
def summary_for_splitted_cluster(cell,name = ""):
	clusters_split = np.load("../Result/%s/clusters%s.npy" %(cell,name))
	node2bin = Load_dict("../Temp/%s/node2bin.npy" %(cell))
	real_cell = cell.split("_")[0]
	density = []
	subgraph_density = []
	distance = []
	TF_num_list = []
	TF_name_list = []

	Build_Result(cell,clusters_split,"../HIMS/%s/Result_split/" %(cell))

	connect = np.load("../Temp/%s/global_connect.npy" %(cell))

	if real_cell != 'nhek':
		master = set(getlist("../Statistics/master/master_tfs_subset_grn_%s.txt" %(cell.split("_")[0])))
	master_in = []
	master_list = []
	for i in xrange(len(clusters_split)):
		cluster = clusters_split[i]
		cluster.sort()

		bins = cluster[cluster>=global_TF_num]
		tfs = cluster[cluster<global_TF_num]

		bin_num = len(bins)
		TF_num = len(tfs)


		temp_list = []
		for tf in tfs:
			temp_list.append(node2bin[tf])


		if real_cell != 'nhek':
			master_in.append(len(list(set(temp_list) & master)))
			master_list.append(";".join(list(set(temp_list) & master)))
		else:
			master_in.append(np.nan)
			master_list.append(np.nan)

		TF_name_list.append(";".join(temp_list))
		TF_num_list.append(TF_num)

		edge = bin_num * (bin_num - 1) * 0.5

		local_connect = connect[bins,:]
		local_connect = local_connect[:,bins]
		if edge <= 0:
			print "error"
			density.append(np.nan)
		else:
			density.append(np.sum(local_connect) * 0.5 / edge)

		location = []
		for Bin in bins:
			location.append(int(node2bin[Bin].split(":")[1]))
		if len(location) > 0:
			left = np.min(location)
			right = np.max(location)
			distance.append((right - left)*1.0/10000)
		else:
			distance.append(0)
			print "error2"
			print cluster

		if "tri" in cell:
			subgraph = edge * TF_num
		else:
			subgraph = edge * TF_num * (TF_num-1) * 0.5
		
		local_connect = np.copy(connect[cluster,:])
		local_connect = local_connect[:,cluster]

		if "tri" in cell:
			subgraph_exist = Motifcount_pdd(local_connect,TF_num)
		else:
			subgraph_exist = Motifcount_npdd(local_connect,TF_num)

		subgraph_density.append(subgraph_exist/subgraph)

	indexes = np.load("../Result/%s/clusters%s_index.npy" %(cell,name))
	cluster_split_summary = {'index':indexes,'edge_density':density,'motif_density':subgraph_density,'gene_distance':distance,'TF_number':TF_num_list,'TF_name':TF_name_list,'master_in':master_in,'master_list':master_list}
	cluster_split_summary = pd.DataFrame(cluster_split_summary,columns = ['index','edge_density','motif_density','gene_distance','TF_number','TF_name','master_in','master_list'])
	cluster_split_summary.to_csv("../Statistics_Result/statistics/%s_cluster_split_summary.txt" %(cell),sep = "\t",index = False)

def summary_for_cell(cell,name = ""):
	result_file = open("../Statistics_Result/statistics/%s_summary.txt" %(cell),"w")
	clusters = np.load("../Result/%s/clusters%s.npy" %(cell,name))
	clusters_split = np.load("../Result/%s/clusters%s_split.npy" %(cell,name))

	node2bin = Load_dict("../Temp/%s/node2bin.npy" %(cell))
	gene_num_all = max(node2bin.keys())
	bin2gene = Load_dict("../Temp/%s/bin2gene.npy" %(cell))


	length = len(clusters_split)
	result_file.write("Number of clusters = %d(before splited)\n" %(len(clusters)))
	result_file.write("Number of clusters = %d(splited)\n" %(length))

	overlapping = 0
	bin_num_counting = 0

	for i in xrange(len(clusters)):
		cluster = clusters[i]
		flag = False
		bins = cluster[cluster>=global_TF_num]
		bin_num_counting += len(bins)
		chrom_now = node2bin[bins[0]].split(":")[0]
		chrom_num = 1
		for j in xrange(1,len(bins)):
			chrom = node2bin[bins[j]].split(":")[0]
			if chrom != chrom_now:
				chrom_num += 1
				chrom_now = chrom
				flag = True
		if flag:
			overlapping += 1

	tempfile = pd.read_table("../Statistics_Result/statistics/%s_cluster_summary.txt" %(cell),sep = "\t")
	num = len(tempfile[(tempfile["chrom_num"] == 1) & (tempfile["TF_num_inherit"] == 0)])
	print num
	result_file.write("Number of clusters with numerous chromosomes = %d\n" %(overlapping))
	result_file.write("Number of genes in a cluster = %d\nNumber of genes not in any cluster = %d\n" %(bin_num_counting,gene_num_all-bin_num_counting))
	result_file.write("Number of clusters with overlapping blah blah(splitted) = %d\n" %(length - num))

def summary_for_TF(cell,name = ""):
	clusters = np.load("../Result/%s/clusters%s.npy" %(cell,name))
	clusters_split = np.load("../Result/%s/clusters%s_split.npy" %(cell,name))

	node2bin = Load_dict("../Temp/%s/node2bin.npy" %(cell))
	gene_num_all = max(node2bin.keys())
	bin2gene = Load_dict("../Temp/%s/bin2gene.npy" %(cell))

	
	connect = np.load("../Temp/%s/global_connect.npy" %(cell))


	TF_overlap = np.zeros((global_TF_num),dtype = "int")
	TF_overlap_cluster = {}
	for i in xrange(len(clusters_split)):
		cluster = clusters_split[i]
		cluster.sort()
		tfs = cluster[cluster<global_TF_num]

		for tf in tfs:
			TF_overlap[tf] += 1
			if TF_overlap_cluster.has_key(tf):
				TF_overlap_cluster[tf] += ";%d" %(i)
			else:
				TF_overlap_cluster[tf] = "%d" %(i)

	TF_name = []
	TF_cluster = []
	for i in xrange(global_TF_num):
		TF_name.append(node2bin[i])
		if TF_overlap_cluster.has_key(i):
			TF_cluster.append(TF_overlap_cluster[i])
		else:
			TF_cluster.append(" ")



	


	TF_summary = {'Overlapping':TF_overlap,'TF_name':TF_name,'TF_overlap_cluster':TF_cluster}
	TF_summary = pd.DataFrame(TF_summary,columns = ['TF_name','Overlapping','TF_overlap_cluster'])
	TF_summary.to_csv("../Statistics_Result/statistics/%s_TF_summary.txt"%(cell),sep = "\t",index = False)

def featureremap(split_table,origin_table,name_list):
	value_list = []
	for i in xrange(len(split_table)):
		a = origin_table.loc[split_table['origin_cluster'][i],name_list]
		value_list.append(a.values)
	value_list = np.array(value_list)

	for i in xrange(len(name_list)):
		split_table[name_list[i]] = value_list[:,i]
	return split_table
	
def merge_for_cluster(cell):
	table = pd.read_table("../Statistics_Result/statistics/%s_cluster_split_summary.txt" %(cell),sep = "\t")
	repli = pd.read_table("../Statistics_Result/repli/%s.txt" %(cell),sep = "\t")
	compartment = pd.read_table("../Statistics_Result/compartment/%s.txt" %(cell),sep = "\t")
	tads = pd.read_table("../Statistics_Result/tads/%s.txt" %(cell),sep = "\t")
	loops = pd.read_table("../Statistics_Result/loops/%s.txt" %(cell),sep = "\t")

	a = loops[['chrom','gene_num','loop_num']]
	b = tads[['tads_num']]
	c = repli[['repli_mean','repli_cv']]
	d = compartment[['A_num','B_num']]
	index = pd.DataFrame(np.arange(len(table)),columns = ['index'])
	table = table.rename(columns={'index':'origin_cluster'})
	table = pd.concat([index,table,a,b,c,d],axis = 1)



	cluster_summary = pd.read_table("../Statistics_Result/statistics/%s_cluster_summary.txt" %(cell),sep = "\t")

	table = featureremap(table,cluster_summary,['TF_number','TF_num_inherit','TF_self',\
												'inherit/self','inherit/all','ppi_num',\
												'ppi_density','ppi_detail','source','ppi_complex',\
												'ppi_complex_detail','ppi_complex_all'])
	#put source to front
	cols = table.columns.tolist()
	cols.remove("source")
	cols.remove("index")
	cols = ["index","source"]+cols 
	print cols
	table = table[cols]
	table.to_csv("../Statistics_Result/new_%s_cluster_split_summary.txt" %(cell),sep = "\t",index = False)

def outputbed(cell):
	clusters = np.load("../Result/%s/clusters_add_left.npy" %(cell))
	not_used = clusters[-2]

	node2bin = Load_dict("../Temp/%s/node2bin.npy" %(cell))
	bin2gene = Load_dict("../Temp/%s/bin2gene.npy" %(cell))

	used_file = open("../Statistics_Result/statistics/%s_used.bed" %(cell),"w")
	not_used_file = open("../Statistics_Result/statistics/%s_not_used.bed" %(cell),"w")

	for i in xrange(len(clusters) - 2):
		cluster = clusters[i]
		cluster.sort()
		bins = cluster[cluster>=global_TF_num]
		for Bin in bins:
			line = node2bin[Bin]
			chrom = line.split(":")[0]
			start = int(line.split(":")[1])
			end = start + 10000
			used_file.write("%s\t%d\t%d\tcluster_%d\n" %(chrom,start,end,i))

	bins = not_used[not_used>=global_TF_num]
	for Bin in bins:
		line = node2bin[Bin]
		chrom = line.split(":")[0]
		start = int(line.split(":")[1])
		end = start + 10000
		not_used_file.write("%s\t%d\t%d\tnot_used\n" %(chrom,start,end))


def compartment(cell):
	result_file = open("../Statistics_Result/compartment/%s.txt" %(cell) ,"w")
	result_file.write("clusters\tgene_num\tA_num\tB_num\n")
	clusters = np.load("../Result/%s/clusters_add_left_split.npy" %(cell))
	node2bin = Load_dict("../Temp/%s/node2bin.npy" %(cell))
	bin2gene = Load_dict("../Temp/%s/bin2gene.npy" %(cell))
	gene_table = pd.read_table("../Data/%s/used/gene_chrom_bin_num_hg19_%s.txt" %(cell,cell.split("_")[0]),sep = "\t")
	gene_table.index = gene_table["TF"]

	compartments = pd.read_table("../Statistics/compartment/compartments_%s.bedGraph" %(cell.split("_")[0]),sep = "\t")
	count = 0
	for cluster in clusters:
		A_count = 0
		B_count = 0
		cluster.sort()
		bin_nodes = cluster[cluster >= 590]

		genes = []
		for node in bin_nodes:
			genes.append(bin2gene[node2bin[node]].split(";")[0])

		chrom_list =  gene_table["chrom"][genes].unique()
		temp1_gene = gene_table.loc[genes]

		for chrom in chrom_list:
			A_count = 0
			B_count = 0
			temp_gene = temp1_gene[temp1_gene["chrom"] == chrom]
			temp_gene.index = xrange(len(temp_gene))
			temp_compartments = compartments[compartments["chrom"] == chrom]

			for i in xrange(len(temp_gene)):
				start = temp_gene["txStart_left"][i]

				a = temp_compartments[(temp_compartments['start']<=start) & (temp_compartments['end']>=start)]
				if len(a) != 0:
					if a["value"][a.index[0]]>=0:
						A_count +=1
					else:
						B_count +=1

			result_file.write("%d\t%d\t%s\t%d\n" %(count,len(temp_gene),A_count,B_count))
		count += 1
		print "%d\r" %(count),
		sys.stdout.flush()

def loops(cell):
	result_file = open("../Statistics_Result/loops/%s.txt" %(cell) ,"w")
	result_file.write("clusters\tgene_num\tchrom\tloop_num\n")

	clusters = np.load("../Result/%s/clusters_add_left_split.npy" %(cell))
	node2bin = Load_dict("../Temp/%s/node2bin.npy" %(cell))
	bin2gene = Load_dict("../Temp/%s/bin2gene.npy" %(cell))
	gene_table = pd.read_table("../Data/%s/used/gene_chrom_bin_num_hg19_%s.txt" %(cell,cell.split("_")[0]),sep = "\t")
	gene_table.index = gene_table["TF"]

	loops = pd.read_table("../Statistics/loops/%s_loop.txt" %(cell.split("_")[0]),sep = "\t")
	count = 0
	for cluster in clusters:
		cluster.sort()
		bin_nodes = cluster[cluster >= 590]
		genes = []
		for node in bin_nodes:
			genes.append(bin2gene[node2bin[node]].split(";")[0])

		chrom_list =  gene_table["chrom"][genes].unique()
		temp1_gene = gene_table.loc[genes]

		for chrom in chrom_list:
			temp_gene = temp1_gene[temp1_gene["chrom"] == chrom]
			temp_gene.index = xrange(len(temp_gene))
			temp_loops = loops[loops["chr1"] == chrom[3:]]

			temp1 = temp_loops[(temp_loops["x1"] <= temp_gene['txStart_left'][0]) & (temp_loops["y2"] >= temp_gene['txStart_left'][0])]
			for i in xrange(1,len(temp_gene)):
				start = temp_gene['txStart_left'][i]
				temp1 = pd.concat([temp1,temp_loops[(temp_loops["x1"] <= start) & (temp_loops["y2"] >= start)]])
			temp1 = temp1.drop_duplicates(keep = 'first')
			result_file.write("%d\t%d\t%s\t%d\n" %(count,len(temp_gene),chrom,len(temp1)))
		count += 1
		#print "%d\r" %(count),
		#sys.stdout.flush()

def tads(cell):

	result_file = open("../Statistics_Result/tads/%s.txt" %(cell) ,"w")
	result_file.write("clusters\tgene_num\tchrom\ttads_num\n")

	clusters = np.load("../Result/%s/clusters_add_left_split.npy" %(cell))
	node2bin = Load_dict("../Temp/%s/node2bin.npy" %(cell))
	bin2gene = Load_dict("../Temp/%s/bin2gene.npy" %(cell))
	gene_table = pd.read_table("../Data/%s/used/gene_chrom_bin_num_hg19_%s.txt" %(cell,cell.split("_")[0]),sep = "\t")
	gene_table.index = gene_table["TF"]

	tads = pd.read_table("../Statistics/tads/%s_tads.txt" %(cell.split("_")[0]),sep = "\t")
	count = 0
	for cluster in clusters:
		cluster.sort()
		bin_nodes = cluster[cluster >= 590]
		genes = []
		for node in bin_nodes:
			genes.append(bin2gene[node2bin[node]].split(";")[0])

		chrom_list =  gene_table["chrom"][genes].unique()
		temp1_gene = gene_table.loc[genes]

		for chrom in chrom_list:
			temp_gene = temp1_gene[temp1_gene["chrom"] == chrom]
			temp_gene.index = xrange(len(temp_gene))
			temp_tads = tads[tads["chr1"] == chrom[3:]]

			
			temp1 = temp_tads[(temp_tads["x1"] <= temp_gene['txStart_left'][0]) & (temp_tads["y2"] >= temp_gene['txStart_left'][0])]
			for i in xrange(1,len(temp_gene)):
				start = temp_gene['txStart_left'][i]
				temp1 = pd.concat([temp1,temp_tads[(temp_tads["x1"] <= start) & (temp_tads["y2"] >= start)]])

			temp1 = temp1.drop_duplicates(keep = 'first')

			result_file.write("%d\t%d\t%s\t%d\n" %(count,len(temp_gene),chrom,len(temp1)))
		count += 1
		print "%d\r" %(count),
		sys.stdout.flush()

def repli(cell):

	result_file = open("../Statistics_Result/repli/%s.txt" %(cell) ,"w")
	result_file.write("cluster_origin\tclusters\tcount\trepli_mean\trepli_std\trepli_cv\trepli_min\trepli_25\trepli_50\trepli_75\trepli_max\n")

	clusters = np.load("../Result/%s/clusters_add_left_split.npy" %(cell))
	index = np.load("../Result/%s/clusters_add_left_split_index.npy" %(cell))
	node2bin = Load_dict("../Temp/%s/node2bin.npy" %(cell))
	bin2gene = Load_dict("../Temp/%s/bin2gene.npy" %(cell))
	gene_table = pd.read_table("../Data/%s/used/gene_chrom_bin_num_hg19_%s.txt" %(cell,cell.split("_")[0]),sep = "\t")
	gene_table.index = gene_table["TF"]

	repli = pd.read_table("../Statistics/repli/repli_seq_%s.txt" %(cell.split("_")[0]),sep = "\t")
	repli.index =repli["gene"]
	count = 0
	for cluster in clusters:
		cluster.sort()
		A_count = 0
		B_count = 0
		bin_nodes = cluster[cluster >= 590]
		genes = []
		for node in bin_nodes:
			genes.append(bin2gene[node2bin[node]].split(";")[0])

		value = repli["repli-seq"][genes].describe()
		result_file.write("%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" %(index[count],count,value[0],value[1],value[1]/value[0],value[2],value[3],value[4],value[5],value[6],value[7]))
		count += 1

		print "%d\r" %(count),
		sys.stdout.flush()

def sub(cell):
	result_file = open("../Statistics_Result/compartment/%ssub.txt" %(cell) ,"w")
	result_file.write("clusters\tnum\tnum2\tA1\tA2\tB1\tB2\tB3\tB4\n")
	clusters = np.load("../Result/%s/clusters_add_left_split.npy" %(cell))
	node2bin = Load_dict("../Temp/%s/node2bin.npy" %(cell))
	bin2gene = Load_dict("../Temp/%s/bin2gene.npy" %(cell))
	gene_table = pd.read_table("../Data/%s/used/gene_chrom_bin_num_hg19_%s.txt" %(cell,cell.split("_")[0]),sep = "\t")
	gene_table.index = gene_table["TF"]

	compartments = pd.read_table("../Statistics/compartment/subcompartments_%s.bed" %(cell),sep = "\t")
	print compartments["value"].unique()
	count = 0
	dict_1 = {'A1':0,'A2':1,'B1':2,'B2':3,'B3':4,'B4':5}
	for cluster in clusters:
		cluster.sort()
		compartment_count = [0,0,0,0,0,0]

		bin_nodes = cluster[cluster >= 590]
		genes = []
		for node in bin_nodes:
			genes.append(bin2gene[node2bin[node]].split(";")[0])

		chrom_list =  gene_table["chrom"][genes].unique()
		num2 = 0
		temp1_gene = gene_table.loc[genes]

		for chrom in chrom_list:
			temp_gene = temp1_gene[temp1_gene["chrom"] == chrom]
			temp_gene.index = xrange(len(temp_gene))
			temp_compartments = compartments[compartments["chrom"] == chrom]

			temp_compartments = temp_compartments.sort_values(['start'])
			temp_compartments.index = xrange(len(temp_compartments))
			length = len(temp_gene)
			temp_gene = temp_gene.sort_values(['txStart_left'])
			temp_gene.index = xrange(len(temp_gene))
			for i in xrange(len(temp_gene)):
				start = temp_gene["txStart_left"][i]

				a = temp_compartments[(temp_compartments['start']<=start) & (temp_compartments['end']>=start)]
				if len(a) != 0:
					if len(a) != 1:
						print len(a)
					if dict_1.has_key(a["value"][a.index[0]]):
						compartment_count[dict_1[a["value"][a.index[0]]]] += 1
						num2+=1
					else:
						print a
						print a["value"][a.index[0]]


		result_file.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" %(count,len(bin_nodes),num2,compartment_count[0],compartment_count[1],compartment_count[2],compartment_count[3],compartment_count[4],compartment_count[5]))
		count += 1
		'''
		print "%d\r" %(count),
		sys.stdout.flush()
		'''

def Everything(cell):
	clusters = Add_Leftover(cell)
	Build_Result(cell,clusters,"../HIMS/%s/Result_parents/" %(cell))
	clusters = Add_Leftover(cell,"_noparent")
	Build_Result(cell,clusters,"../HIMS/%s/Result_noparents/" %(cell))

	Split_chromosome(cell,"")
	Split_chromosome(cell,"_add_left")
	Split_chromosome(cell,"_noparent_add_left")
	
	compartment(cell)
	loops(cell)
	tads(cell)
	repli(cell)
	if 'gm12878' in cell:
		sub(cell)

	#新加的cluster中一个是剩下的，另一个是用过的
	summary_for_cluster(cell,"_add_left")
	summary_for_splitted_cluster(cell,"_add_left_split")
	summary_for_TF(cell,"")
	summary_for_cell(cell,"")
	#outputbed(cell)

	merge_for_cluster(cell)

#cell_list = ['gm12878','hela','huvec','nhek','k562']
#cell_list = ['gm12878_2','gm12878_4','gm12878_6','gm12878_weibull_0.5','gm12878_weibull_0.75','gm12878_weibull_0.9','gm12878_weibull_0.95']
#cell_list = ['gm12878','hela','huvec','nhek','k562','gm12878_2','gm12878_4','gm12878_6','gm12878_weibull_0.5','gm12878_weibull_0.75','gm12878_weibull_0.9','gm12878_weibull_0.95','gm12878_grn_2','hela_grn_2','huvec_grn_2','nhek_grn_2','k562_grn_2']
#cell_list = ['gm12878_grn_2','hela_grn_2','huvec_grn_2','nhek_grn_2','k562_grn_2']
#cell_list = ['gm12878_grn_3','hela_grn_3','huvec_grn_3','nhek_grn_3','k562_grn_3']
