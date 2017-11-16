#encoding:utf-8
import numpy as np
import pandas as pd
from numpy.linalg import eigh
from Build_Result import Build_Result
import shutil
import sys
import os
import pickle
from optparse import OptionParser
import time
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
	def __init__(self, edge_list = None,connect = None, edge_weight = 1,projection = None,cell = ""):  
		try:
			self.edge2connect(edge_list)
			#Projection means the ith row in W is the node k
			self.projection = np.array(range(self.connect.shape[0]))
			print "Connection Initialized"
		except:
			try:
				self.connect = connect
				self.projection = projection
			except:
				print "Error"

		if edge_weight == 0:
			self.W = np.zeros_like(self.connect)
		else:
			self.W = np.copy(self.connect) * edge_weight
		self.set_weight()

		self.parent = None
		self.cell = cell

		self.TF = self.projection[self.projection < 590]
		self.Bin = self.projection[self.projection >= 590]


		self.interval_table = pd.read_table("../Temp/%s/chrom_bin_index.txt" %(cell),sep = "\t")

	def drop_zero(self):
		a  = np.sum(self.W,axis = 0)
		b = (a > 0)
		self.projection = self.projection[b]
		self.W = self.W[b,:]
		self.W = self.W[:,b]

		self.connect = self.connect[b,:]
		self.connect= self.connect[:,b] 

		self.TF = self.projection[self.projection < 590]
		self.Bin = self.projection[self.projection >= 590]


		self.density = []
		for i in xrange(len(self.interval_table)):
			start = self.interval_table['index_start'][i]
			end = self.interval_table['index_end'][i]

			Bin_num = np.sum((self.projection < end) & (self.projection >= start))
			if Bin_num >= 10:
				index_we_care = ((self.projection < end) & (self.projection >= start))
				self.total_num = Bin_num * (Bin_num - 1) * len(self.TF)*(len(self.TF)-1) * 0.25
				
				a = self.W[index_we_care,:]
				a = a[:,index_we_care]
				self.density.append(np.sum(a) * 1.0 / (2* self.total_num))
		if len(self.density) == 0:
			self.density = 1 
		else:
			self.density = np.max(self.density)

		#print "density!"
		#print self.density
		#self.total_num = len(self.TF)*(len(self.TF)-1)*len(self.Bin)*(len(self.Bin)-1)*0.25
		#self.density = np.sum(self.W) * 1.0 / 12 / self.total_num

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
		if N >= 3000: 
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
		#相当于Weight的不重复计数,i >= j
		W_low = np.tril(W,k = -1)
		#i 有多少motif
		W_sum = np.sum(W,axis = 1)
		#i 和比自己小的有多少motif
		W_lower_sum = np.sum(W_low,axis = 1)
		       
		#到目前为止，能形成的motif个数
		volumes = np.cumsum(W_sum)
		#减去两份是因为，之前i-1的时候已经算了一次，要补回来 
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
		self.TF_num1 = np.sum(self.projection[self.index1] < 590)
		self.TF_num2 = np.sum(self.projection[self.index2] < 590)
		
		
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

		tra_count_dna_left = tria_count_dna - 1
		tra_count_dna_left = np.where(tra_count_dna_left<0,0,tra_count_dna_left)
		temp = (pro_connect.dot(tra_count_dna_left))
		pro_weight = (temp+temp.T) * pro_connect * self.pd_weight

		
		pro_weight2 = np.zeros_like(self.W)
		
		#加速版：
		for i in xrange(split_point):
			#太多的情况下还是给个进度条
			if split_point >= 50:    
				print "Finishing:%d\r" %(i),
				sys.stdout.flush()

			share_dna = (pro_connect[i,:] + pro_connect[i+1:split_point,:])/2
			share_dna = np.where(share_dna < 1,0,share_dna)
			num = share_dna.dot(dna_connect)
			#除以2是因为，k与l的组合在k和l分别计算了一次
			num = np.sum(num * share_dna,axis= 1) / 2
			pro_weight2[i,i+1:split_point] = num * self.pp_weight
			pro_weight2[i+1:split_point,i] = num.T * self.pp_weight
		
		self.W = self.W + dna_weight + pro_weight + pro_weight2


def Save_Result(cell,clusters,path):
	result_list = []
	count = 0
	for graph in clusters:
		result_list.append(np.concatenate((graph.Bin, graph.TF),axis =-1))
		#np.save("./Matrix/"+str(count)+".npy",graph.sth)
		count += 1
	result = np.asarray(result_list)
	np.save("../Result/%s/clusters.npy" %(cell),result)

	Build_Result(cell,result,path)
	return result

def Dump_Graph(clusters):
	count = 0
	for graph in clusters:
		f = open("./"+str(count)+".pkl","wb")
		pickle.dump(graph,f)
		count += 1

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
		for i in xrange(len(a)):
			if a[i] >= 2:
				break
		graph.TF = np.unique(np.concatenate((graph.TF[:i],backup_TF),axis = -1))

def Multi_Spectral(clusters,threshold):

	score_list = []
	length  = len(clusters)
	for graph in clusters:
		score_list.append(graph.density)

	#选出如果被切了，score最低的
	index = np.argmin(np.array(score_list))
	min_score = score_list[index]
	print "min_score = %f" %(min_score)
	#都不怎么能切了就不切
	if min_score >= threshold:
		return

	graph = clusters[index]
	print "Node %d to %d / %d, TF_num %d to %d / %d" %(len(graph.connect),len(graph.index1),len(graph.index2),graph.TF_num,graph.TF_num1,graph.TF_num2)

	for index in [graph.index1,graph.index2]:
		index.sort()

		connect = graph.connect[index,:]
		connect = connect[:,index]
		#Projection是global的，都表示的最终映射到哪个节点
		projection = graph.projection[index]
		graph_child = PD_Graph(connect = connect,edge_weight = 0,projection = projection,cell = graph.cell)
		
		graph_child.MotifAdjacency_npdd()
		graph_child.drop_zero()
		'''
		if graph_child.density <= graph.density:
			print graph_child.density
			print graph.density
			graph.score += 1
			break
		'''
		#如果这个graph中没node了就不再放入list
		if graph_child.W.shape[0] != 0:
			graph_child.parent = graph
			graph_child.SpectralPartition()
			clusters.append(graph_child)

	if len(clusters) > length:
		clusters.remove(graph)
	Multi_Spectral(clusters,threshold)

def run(cell,stop_threshold,parent_threshold,total):
	edge_list = read_txt("../Data/%s/Edge_list.txt" %(cell))
	graph = PD_Graph(edge_list = edge_list,edge_weight = 0,cell = cell)
	np.save("../Temp/%s/global_connect.npy" %(cell),graph.connect)
	global_connect = np.copy(graph.connect)
	print graph.connect.shape

	if (os.path.isfile("../Data/%s/W.npy" %(cell))) & (total == False):
		graph.W = np.load("../Data/%s/W.npy" %(cell))
	else:
		graph.MotifAdjacency_npdd()
		np.save("../Data/%s/W.npy" %(cell),graph.W)
	'''
	bin_num = len(graph.W) - 590
	whole_graph_density = np.sum(graph.W[590:,590:]) * 0.5 / (590*589*0.5*bin_num*(bin_num-1)*0.5)
	print "whole_graph_density:%f\n" %(whole_graph_density)
	'''
	graph.drop_zero()
	graph.SpectralPartition()

	print graph.W.shape

	
	clusters = [graph]
	Multi_Spectral(clusters,stop_threshold)

	try:
		shutil.rmtree("../Result/%s/Result_noparents/" %(cell))
		shutil.rmtree("../Result/%s/Result_parents/" %(cell))
	except:
		print "First time"

	os.makedirs("../Result/%s/Result_noparents/" %(cell))
	os.makedirs("../Result/%s/Result_parents/" %(cell))

	result = Save_Result(cell,clusters,"../Result/%s/Result_noparents/" %(cell))
	np.save("../Result/%s/clusters_noparent.npy" %(cell),result)
	Find_Parent(clusters,global_connect,parent_threshold)
	Save_Result(cell,clusters,"../Result/%s/Result_parents/" %(cell))
