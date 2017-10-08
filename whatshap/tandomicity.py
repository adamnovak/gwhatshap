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

#bubbles_start = set()
#with stream.open('assembly_graph.P.int.remn2n.X_100.trans' ,"rb") as istream:
#	for data in istream:
#		l = vg_pb2.SnarlTraversal()
#		l.ParseFromString(data)
#		if l.snarl.start.backward == True:
#			start_node = l.snarl.end.node_id
#		else:
#			start_node = l.snarl.start.node_id
#		bubbles_start.add(start_node)
#print(len(bubbles_start))

nodes = set()
edge_connections = defaultdict(set)

gfafile = open('assembly_graph.P.int.remn2n.X_100.view.gfa', "rb")
for line in gfafile:
	var = line.split('\t')
	if var[0] == 'S':
		nodes.add(int(var[1]))
	if var[0] == 'L':
		edge_connections[int(var[1])].add(int(var[3]))


multiplicity_bubbles = defaultdict(list)
read_details = defaultdict(list)
with stream.open('out.new.gam', "rb") as istream:
	for data in istream:
		g = vg_pb2.Alignment()
		g.ParseFromString(data)
		for i in range(0,len(g.path.mapping)):
			node = g.path.mapping[i].position.node_id
			if node in nodes:
				multiplicity_bubbles[node].append(g.name)
			tmp = '+'
			if g.path.mapping[i].position.is_reverse == 'True':
				tmp = '-'
			node_tmp = str(node)+ "_" + str(tmp)
			read_details[g.name].append(node_tmp)

count=0
repeaticity = defaultdict()
repeaticity_read_support = defaultdict()
repeaticity_i = defaultdict()
repeaticity_read_support_i = defaultdict()
for k,v in multiplicity_bubbles.items():
	tmp = defaultdict(set)
	for item,count in collections.Counter(v).items():
		tmp[count].add(item)
	max_val = 0 
	item_val = set()
	l1_max = 0
	# how many reads support multiplicity l1
	for l1,l2 in tmp.items():
		if max_val < len(l2):
			max_val = len(l2)
			item_val = l2
			l1_max = l1
	if max_val > 5 and l1_max > 1: # minimum read support and minimum multiplicty
		repeaticity[k] = l1_max
		repeaticity_read_support[k] = item_val
		print(item_val, max_val, l1_max, k, v)
	if len(set(v)) > 38:
		cnv = len(set(v))/38.0
		repeaticity_i[k] = cnv
		repeaticity_read_support_i[k] = set(v)
		print('interpersed repeat')
		true_chrms = set()
		for i in set(v):
			var = var = i.split("_")[2]
			true_chrms.add(var)
		print(set(v), k, cnv, len(true_chrms))
		#if int(round(cnv)) == len(true_chrms):
		#	count=count+1

gfaout = open('repeastresolved.gfa', 'w')
gfafile = open('assembly_graph.P.int.remn2n.X_100.view.gfa', "rb")
all_ready_done = set()

for 
for line in gfafile:
	var = line.split('\t')
	if var[0] == 'L':
		if int(var[1]) in repeaticity: # if node in repated list
			all_ready_done.add(int(var[1]))
			for neighbors in edge_connections(int(var[1])):
			for read in repeaticity_read_support[int(var[1])]:				
				str_tmp1 = var[1]+'+'
				str_tmp2 = var[1]+'-'
				for i,j in zip(read_details[read], read_details[read][1:]):
					gfaout.write('L' + "\t" + i.split("_")[0] + "\t" + i.split("_")[1] + j.split("_")[0] + "\t" + j.split("_")[1] + "\n")
				indices1 = [i for i, x in enumerate(read_details[read]) if x == str_tmp1] # which reads support the repeat and find indices in the list
				indices2 = [i for i, x in enumerate(read_details[read]) if x == str_tmp2] # which reads support the repeat and find indices in the list
			
					


# Now find the shortest common superstring between strings with maximum multiplicity.
		

	

