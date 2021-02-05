#-*- coding: utf-8 -*-  
from optparse import OptionParser
import os
import Preprocess
import Graph

def parse_args():
	parser = OptionParser(usage="Graph Partition", add_help_option=False)
	parser.add_option("-c","--cell", type="string", default='gm12878')
	parser.add_option("-s", "--stop_threshold", type="float", default=0.7)
	parser.add_option("-p", "--parent_threshold", type="float", default=0.45)
	parser.add_option("-t", "--total", type="string", default="False")
	parser.add_option("-m", "--motif", type="int", default=4)

	(opts, args) = parser.parse_args()
	return opts

def makepath(cell):
	print ("making path")
	if not os.path.exists("../Temp/%s/" %(cell)):
		os.makedirs("../Temp/%s/" %(cell))

	if not os.path.exists("../Result/%s/" %(cell)):
		os.makedirs("../Result/%s/" %(cell))


	for result in ['compartment','tads','loops','repli','statistics']:
		if not os.path.exists("../Statistics_Result/%s" %(result)):
			os.makedirs("../Statistics_Result/%s" %(result))



def run(cell,stop_threshold,parent_threshold,total,motif):
	#Build the Path
	makepath(cell)
	if (not os.path.isfile("../Data/%s/Edge_list.txt" %(cell))) or total:
		Preprocess.run(cell)

	Graph.run(cell,stop_threshold,parent_threshold,total,motif)








def main():
	opts = parse_args()

	#Quick workaround because optparse doesn't natively support boolean types.
	#The cmdline interface should eventually be changed so that putting "-t" makes opts.total True, otherwise its False.
	#For now though, want to maintain consistency with pre-existing documentation and usage.
	if (opts.total == "True"):
		opts.total = True
	else:
		opts.total = False

	run(opts.cell,opts.stop_threshold,opts.parent_threshold,opts.total,opts.motif)

if __name__ == '__main__':
	main()