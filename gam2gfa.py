import sys
import stream
import logging
import vg_pb2
from collections import Counter
from collections import defaultdict
import collections
from collections import OrderedDict, namedtuple
from collections import defaultdict

file_input = sys.argv[1]
#file_out = argv[2]
out = open(file_input + '.gfa', 'w')

nodes_list = set()
with stream.open(str(file_input), "rb") as istream:
	for data in istream:
		g = vg_pb2.Alignment()
		g.ParseFromString(data) 
		for i in range(0,len(g.path.mapping)-1):
			node1 = g.path.mapping[i].position.node_id
			node2 = g.path.mapping[i+1].position.node_id
			nodes_list.add(node1)
			nodes_list.add(node2)
			if g.path.mapping[i].position.is_reverse == True and g.path.mapping[i+1].position.is_reverse == True:
				out.write("L" + "\t" + str(node1)+"\t" + '-' + "\t" + str(node2) + "\t" + '-' +"\t" + "0M" + "\n")
			if g.path.mapping[i].position.is_reverse == False and g.path.mapping[i+1].position.is_reverse == True:
				out.write("L" + "\t" + str(node1)+"\t" + '+' + "\t" + str(node2) + "\t" + '-' +"\t" + "0M" + "\n")
			if g.path.mapping[i].position.is_reverse == True and g.path.mapping[i+1].position.is_reverse == False:
				out.write("L" + "\t" + str(node1)+"\t" + '-' + "\t" + str(node2) + "\t" + '+' +"\t" + "0M" + "\n")
			if g.path.mapping[i].position.is_reverse == False and g.path.mapping[i+1].position.is_reverse == False:
				out.write("L" + "\t" + str(node1)+"\t" + '+' + "\t" + str(node2) + "\t" + '+' +"\t" + "0M" + "\n")


for i in nodes_list:
	out.write("S" + "\t" + str(i) + "\t" + "A" + "\n")
	

			

