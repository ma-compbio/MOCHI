#-*- coding: utf-8 -*-  
from optparse import OptionParser
import os
import Preprocess
import Graph

def parse_args():
	parser = OptionParser(usage="Graph Partition", add_help_option=False)
	parser.add_option("-c","--cell",default = 'gm12878')
	parser.add_option("-s", "--stop_threshold", default = 0.7)
	parser.add_option("-p", "--parent_threshold",default = 0.45)
	parser.add_option("-t", "--total",default = False)
	parser.add_option("-m", "--motif",default = 4)

	(opts, args) = parser.parse_args()
	return opts

def makepath(cell):
	if os.path.exists("../Temp/%s" %(cell)) == False:
		os.makedirs("../Temp/%s" %(cell))

	if os.path.exists("../Result/%s" %(cell)) == False:
		os.makedirs("../Result/%s" %(cell))


	for result in ['compartment','tads','loops','repli','statistics']:
		if os.path.exists("../Statistics_Result/%s" %(result)) == False:
			os.makedirs("../Statistics_Result/%s" %(result))



def run(cell,stop_threshold,parent_threshold,total,motif):

	#Build the Path
	makepath(cell)

	if (not os.path.isfile("../Data/%s/Edge_list.txt" %(cell))) or total:
		Preprocess.run(cell)

	Graph.run(cell,stop_threshold,parent_threshold,total,motif)








def main():
	opts = parse_args()
	run(opts.cell,opts.stop_threshold,opts.parent_threshold,opts.total,opts.motif)
if __name__ == '__main__':
	main()