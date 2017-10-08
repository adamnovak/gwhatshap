import sys
import stream
import logging
import vg_pb2
from collections import Counter
from collections import defaultdict
import collections
from collections import OrderedDict, namedtuple
from collections import defaultdict
# assumption ... all S' and before L's
#filename = sys.argv[1]
#out = sys.argv[2]

d={}
count=1



def vg_graph_reader(vg_file):
	node_seq_list= defaultdict()
	edge_connections_from = defaultdict(set)
	edge_connections_to = defaultdict(set)
	with stream.open(str(vg_file), "rb") as istream:
		for data in istream:
			l = vg_pb2.Graph()
			l.ParseFromString(data)
			for i in range(len(l.node)):
				index = l.node[i].id
				seq = l.node[i].sequence
				node_seq_list[index]=seq
			for j in range(len(l.edge)):
				from_edge = getattr(l.edge[j], "from")
				edge_connections_from[from_edge].add(l.edge[j].to)
				edge_connections_to[l.edge[j].to].add(from_edge)
	return node_seq_list, edge_connections_from, edge_connections_to


coverage_monitor = defaultdict(int)
with stream.open('out.new.gam', "rb") as istream:
	for data in istream:
		g = vg_pb2.Alignment()
		g.ParseFromString(data)
		for i in range(0,len(g.path.mapping)):
			node = g.path.mapping[i].position.node_id
			coverage_monitor[node]+=1

													
count1 =0
count =0
node_seq_list, edge_connections_from, edge_connections_to = vg_graph_reader('component0.vg')
for k,v in node_seq_list.items():
	#print(v)
	if len(edge_connections_from[k])==0 and coverage_monitor[k] > 30:
		print(coverage_monitor[k])
		count=count+1
	if len(edge_connections_to[k]) ==0 and coverage_monitor[k] > 30:
		print('hello')
		print(coverage_monitor[k])
		count1=count1+1
	#if len(edge_connections_to[k]) ==0 and len(edge_connections_from[k]) ==0:
	#	count1=count1+1
		
print(count)	
print(count1)

