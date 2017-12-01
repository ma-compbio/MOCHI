import pandas
import numpy as np
import os
def Load_dict(path):
	return np.load(path).item()

def Build_Result(cell,cluster_result,path = "./Result/"):
	try:
		shutil.rmtree(path)
	except:
		print "First Time"
	os.makedirs(path)

	#HiC_Value = np.load("../Temp/%s/HiC_Value.npy" %(cell))
	Node2Bin = Load_dict("../Temp/%s/node2bin.npy" %(cell))
	Bin2Gene = Load_dict("../Temp/%s/bin2gene.npy" %(cell))
	connect = np.load("../Temp/%s/global_connect.npy" %(cell))
	connect = np.triu(connect,k = 1)
	print "Number of clusters = %d" %(len(cluster_result))
	size_list = []
	print len(cluster_result)
	for i in xrange(len(cluster_result)):
		cluster = cluster_result[i]
		size_list.append(len(cluster))
		f = open(path+str(i)+".txt","w")
		f.write("Source\tTarget\tWeight\tEdgeType\tSourceType\tTargetType\n")
		cluster = np.asarray(cluster)
		cluster = np.sort(cluster)
		cluster = cluster.astype(int)

		local_connect = connect[cluster,:]
		local_connect = local_connect[:,cluster]

		targetindex = np.nonzero(local_connect)
		for k in xrange(len(targetindex[0])):
			SourceType = Node2Bin[cluster[targetindex[0][k]]]
			TargetType = Node2Bin[cluster[targetindex[1][k]]]
			
			Target = Bin2Gene[TargetType]
			if cluster[targetindex[0][k]] < 590:
				Source = SourceType
				EdgeType = 'D'
			else:
				Source = Bin2Gene[SourceType]
				EdgeType = 'U'
			#Weight = max(HiC_Value[cluster[targetindex[0][k]],cluster[targetindex[1][k]]],HiC_Value[cluster[targetindex[1][k]],cluster[targetindex[0][k]]])
			Weight = 0
			f.write("%s\t%s\t%f\t%s\t%s\t%s\n" %(Source,Target,Weight,EdgeType,SourceType,TargetType))
		f.close()
