#encoding:utf-8
import numpy as np
import pandas as pd
from numpy.linalg import eigh
#from Build_Result import Build_Result #don't need this anymore
import shutil
import sys
import os
import pickle
from optparse import OptionParser
import time

global_TF_num = 591

def read_txt(path):
	file = open(path,"r")
	lines = file.readlines()
	result_list = []
	lines = lines[1:]
	for line in lines:
		coor = tuple([int(i) for i in line.strip().split("\t")])
		result_list.append(coor)
	return result_list

class PD_Graph():
	def __init__(self, edge_list = None,connect = None, edge_weight = 1,projection = None,cell = "",motif=4,parent=None):  
		try:
			self.edge2connect(edge_list)
			#Projection means the ith row in W is the node k in original dict
			self.projection = np.array(range(self.connect.shape[0]))
			print ("Connection Initialized")
		except:
			try:
				self.connect = connect
				self.projection = projection
			except:
				print ("Error")

		if edge_weight == 0:
			self.W = np.zeros_like(self.connect)
		else:
			self.W = np.copy(self.connect) * edge_weight
		self.set_weight()

		self.parent = parent
		self.cell = cell

		self.TF = self.projection[self.projection < global_TF_num]
		self.Bin = self.projection[self.projection >= global_TF_num]
		
		try:
			self.motif = parent.motif
			self.interval_table = parent.interval_table
		except:
			print ("Root Graph")
			self.motif = motif
			self.interval_table = pd.read_table("../Temp/%s/chrom_bin_index.txt" %(cell),sep = "\t")

	def drop_zero(self):
		a  = np.sum(self.W,axis = 0)
		b = (a > 0)
		self.projection = self.projection[b]
		self.W = self.W[b,:]
		self.W = self.W[:,b]

		self.connect = self.connect[b,:]
		self.connect= self.connect[:,b] 

		self.TF = self.projection[self.projection < global_TF_num]
		self.Bin = self.projection[self.projection >= global_TF_num]
		motif = self.motif

		if motif == 4:
			self.total_num = len(self.Bin) * (len(self.Bin) - 1) * len(self.TF) * (len(self.TF)-1) * 0.25
		elif motif == 3:
			self.total_num = len(self.Bin) * (len(self.Bin) - 1) * len(self.TF) * 0.25
		a = self.W[len(self.TF):,len(self.TF):]
		self.density = np.sum(a) * 0.5 / (self.total_num)

	def auto_split(self):
		split_bins = []
		left_bins = []

		check_num = 0
		drop_list = []
		for i in range(len(self.interval_table)):
			start = self.interval_table['index_start'][i]
			end = self.interval_table['index_end'][i]

			index_we_care = ((self.projection < end) & (self.projection >= start))
			check_num += np.sum(index_we_care)
			if np.sum(index_we_care) != 0:
				a = self.W[index_we_care,:]
				connect_chr_all = np.sum(a[:,global_TF_num:])
				a = a[:,index_we_care]
				connect_chr_self = np.sum(a)

				if connect_chr_all == connect_chr_self:
					split_bins.append(np.arange(len(self.W))[index_we_care])
				else:
					left_bins.extend(list(np.arange(len(self.W))[index_we_care]))
			else:
				drop_list.append(i)

		self.interval_table = self.interval_table.drop(drop_list,axis = 0)
		self.interval_table.index = range(len(self.interval_table))

		left_bins.extend(range(len(self.TF)))
		result_index = []
		for indexs in split_bins:
			temp_list = list(indexs)
			temp_list.extend(range(len(self.TF)))
			result_index.append(np.array(temp_list))
			
		
		result_index.append(np.array(left_bins))
		if check_num != len(self.Bin):
			print ("Bin num error")
			print (check_num)
			print (len(self.Bin))
		return result_index

	def set_weight(self, pd_weight = 1,dd_weight = 1,pp_weight = 1):
		self.pp_weight = pp_weight
		self.pd_weight = pd_weight
		self.dd_weight = dd_weight

	#When W == None, the function will create a new one
	#Otherwise, it will update the existing W
	#By default it will create an undirected graph
	def edge2connect(self, edge_list):
		N = np.max(edge_list) + 1

		self.connect = np.zeros((N,N),dtype = "float")

		for coor in edge_list:
			self.connect[coor] += 1
			self.connect[coor[::-1]] += 1

	def SpectralPartition(self):
		W = np.copy(self.W)
		N = W.shape[0]
		D = np.diag(np.sum(W,axis = 0)* 1.0)
		D_2 = np.diag((np.sum(W,axis = 0)*1.0) ** -0.5)
		self.TF_num = len(self.TF)

		#get Laplacian Matrix
		L = np.dot(D_2,D-W).dot(D_2)

		#print "Get Eigen"
		#Eigen Value and Vector
		cell = self.cell
		str_N = str(N)
		if N >= 5000: 
			if os.path.isfile("../Data/%s/eigen_w%s.npy" %(cell,str_N)):
				w = np.load("../Data/%s/eigen_w%s.npy" %(cell,str_N))
				v = np.load("../Data/%s/eigen_v%s.npy" %(cell,str_N))
			else:
				w,v = eigh(L)
				np.save("../Data/%s/eigen_w%s.npy" %(cell,str_N),w)
				np.save("../Data/%s/eigen_v%s.npy" %(cell,str_N),v)
		else:
			w,v = eigh(L)

		#Eigen Vector is transposed in this function
		v= v.T

		#Sort and get the second smallest...
		index = np.argsort(w)
		w = w[index]
		v = v[index]
		z = v[1]

		#Get spectral order
		spectral_order_value = D_2.dot(z)
		spectral_order = np.argsort(spectral_order_value)
		#Sort Weight Matrix based on spectral_order
		W = W[spectral_order,:]
		W = W[:,spectral_order]
		
		W_low = np.tril(W,k = -1)
		W_sum = np.sum(W,axis = 1)
		W_lower_sum = np.sum(W_low,axis = 1)
		       
		volumes = np.cumsum(W_sum)
		a = W_sum - 2*W_lower_sum
		num_cut = np.cumsum(a) 
		total_vol = np.sum(W)
		volumes_other = total_vol - volumes
		vols = np.min([volumes,volumes_other],axis = 0)
		score = num_cut[:-1] / vols[:-1]
		
		min_index = np.argmin(score)
		self.score = score[min_index]
		self.index1 = spectral_order[0:min_index+1]
		self.index2 = spectral_order[min_index+1:]
		self.TF_num1 = np.sum(self.projection[self.index1] < global_TF_num)
		self.TF_num2 = np.sum(self.projection[self.index2] < global_TF_num)
		self.index1.sort()
		self.index2.sort()
		
		#Don't split it too much
		if (len(self.index1) - self.TF_num1 <= 10) or (len(self.index2) - self.TF_num2 <= 10) or (self.TF_num1 < 2) or (self.TF_num2) < 2:
			self.score += 1
			self.density += 1
		
	#By default, the protein index should be smaller than the bin index
	#protein-dna-dna motif
	#weight_list is for weight of protein-dna and dna-dna interaction in motif
	#split_point means:[0:split_point] = protein, [split_point:] = dna

	#Definition is the same,it makes multi-triangle more robust


	def MotifAdjacency_npdd(self):
		if self.motif == 4:
			#print "Start Counting"
			split_point = len(self.TF)
			connect = self.connect
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
			dna_weight = weight_count_dna * self.dd_weight

			tria_count_dna_left = tria_count_dna - 1
			tria_count_dna_left = np.where(tria_count_dna_left<0,0,tria_count_dna_left)
			temp = (pro_connect.dot(tria_count_dna_left))
			pro_weight = (temp+temp.T) * pro_connect * self.pd_weight

			
			pro_weight2 = np.zeros_like(self.W)
			
			for i in range(split_point):
				if split_point >= 100:    
					print ("Finishing:%d\r" %i, end="")
					sys.stdout.flush()

				share_dna = (pro_connect[i,:] + pro_connect[i+1:split_point,:])/2
				share_dna = np.where(share_dna < 1,0,share_dna)
				num = share_dna.dot(dna_connect)
				num = np.sum(num * share_dna,axis= 1) / 2
				pro_weight2[i,i+1:split_point] = num * self.pp_weight
				pro_weight2[i+1:split_point,i] = num.T * self.pp_weight
			
			self.W = self.W + dna_weight + pro_weight + pro_weight2
		elif self.motif == 3:
			split_point = len(self.TF)
			connect = self.connect

			#Only focus on pro - dna interaction
			pro_connect = np.copy(connect)
			pro_connect[split_point:,split_point:] = 0
			#Only dna-dna interaction
			dna_connect = connect - pro_connect

			dna_weight = pro_connect.dot(pro_connect) * dna_connect * self.dd_weight
			temp = pro_connect.dot(dna_connect)
			pro_weight = (temp+temp.T) * pro_connect  * self.pd_weight
			self.W = self.W + dna_weight + pro_weight

def Save_Result(clusters):
	result_list = []
	count = 0
	for graph in clusters:
		result_list.append(np.concatenate((graph.TF,graph.Bin),axis =-1))
	result = np.asarray(result_list)
	return result

def Find_Parent(clusters,connect,threshold):

	for graph in clusters:
		parent = graph
		while parent.score >= threshold:
			parent = parent.parent
		TF = parent.TF
		backup_TF = np.copy(graph.TF)
		graph.TF = TF
		graph.projection = np.unique(np.concatenate((graph.TF,graph.Bin),axis = -1))
		graph.connect = connect[graph.projection,:]
		graph.connect = graph.connect[:,graph.projection]
		graph.W = np.zeros_like(graph.connect)
		graph.MotifAdjacency_npdd()
		graph.drop_zero()
		graph.sth = graph.W[0:len(graph.TF),0:len(graph.TF)]
		a = np.sum(graph.sth,axis = 1)
		index = np.argsort(a)[::-1]
		a = a[index]
		graph.TF = graph.TF[index]
		a = a[0:-1] / a[1:]
		for i in range(len(a)):
			if a[i] >= 2:
				break
		graph.TF = np.unique(np.concatenate((graph.TF[:i],backup_TF),axis = -1))

def Multi_Spectral(clusters,threshold):
	length  = len(clusters)

	score_list = []
	for graph in clusters:
		score_list.append(graph.score)

	density_list = []
	for graph in clusters:
		density_list.append(graph.density)

	#pick the one with lowest density
	index = np.argmin(np.array(density_list))
	min_score = density_list[index]
	corresponding_score = score_list[index]
	print ("min_density = %f\tscore = %f" %(min_score,corresponding_score))

	#Stop after reaching threshold
	if min_score >= threshold:
		return

	graph = clusters[index]
	print ("Node %d to %d / %d, TF_num %d to %d / %d" %(len(graph.connect),len(graph.index1),len(graph.index2),graph.TF_num,graph.TF_num1,graph.TF_num2))
	child_list = []
	for index in [graph.index1,graph.index2]:
		index.sort()
		connect = graph.connect[index,:]
		connect = connect[:,index]
		#Projection represents global node id
		projection = graph.projection[index]
		graph_child = PD_Graph(connect = connect,edge_weight = 0,projection = projection,cell = graph.cell,parent=graph)
		
		graph_child.MotifAdjacency_npdd()
		graph_child.drop_zero()
		child_list.append(graph_child)

	if (child_list[0].density < graph.density) & (child_list[1].density < graph.density):
		graph.score += 1
		graph.density += 1
		print ("Skip this")
	else:
		for graph_child in child_list:
			indexs = graph_child.auto_split()

			if len(indexs) == 1:
				if graph_child.W.shape[0] != 0:
					graph_child.SpectralPartition()
					clusters.append(graph_child)
			else:
				print ("\nSplit into %d" %(len(indexs)))
				for index_split in indexs:
					index_split.sort()
					connect = graph_child.connect[index_split,:]
					connect = connect[:,index_split]
					projection = graph_child.projection[index_split]
					graph_child_split = PD_Graph(connect = connect,edge_weight = 0,projection = projection,cell = graph_child.cell,parent=graph)
					graph_child_split.MotifAdjacency_npdd()
					graph_child_split.drop_zero()

					if graph_child_split.W.shape[0] != 0:
						graph_child_split.parent = graph
						graph_child_split.SpectralPartition()
						clusters.append(graph_child_split)
			

	if len(clusters) > length:
		clusters.remove(graph)

	Multi_Spectral(clusters,threshold)

def run(cell,stop_threshold,parent_threshold,total,motif):
	edge_list = read_txt("../Data/%s/Edge_list.txt" %(cell))
	graph = PD_Graph(edge_list = edge_list,edge_weight = 0,cell = cell,motif = motif)
	np.save("../Temp/%s/global_connect.npy" %(cell),graph.connect)
	global_connect = np.copy(graph.connect)
	print (graph.connect.shape)

	if (os.path.isfile("../Data/%s/W.npy" %(cell))) & (total == False):
		graph.W = np.load("../Data/%s/W.npy" %(cell))
	else:
		graph.MotifAdjacency_npdd()
		np.save("../Data/%s/W.npy" %(cell),graph.W)

	graph.drop_zero()
	graph.SpectralPartition()
	print (graph.W.shape)

	
	clusters = [graph]
	Multi_Spectral(clusters,stop_threshold)

	

	result = Save_Result(clusters)
	np.save("../Result/%s/clusters_noparent.npy" %(cell),result,allow_pickle=True)
	Find_Parent(clusters,global_connect,parent_threshold)
	result = Save_Result(clusters)
	np.save("../Result/%s/clusters.npy" %(cell),result,allow_pickle=True)

	from statistics import Everything
	Everything(cell)