import sys
import networkx as nx
from collections import defaultdict

G=nx.DiGraph()

filename = sys.argv[1]
out = sys.argv[3]
variants_file = sys.argv[2]

output= open(out, 'w')
with open(filename) as fp:
	for line in fp:
		var=line.split(	)
		if(var[0] == 'S'):
			G.add_node(var[1])
		else:
			G.add_edge(var[1],var[3])


def longest_path(G):
	dist = {} # stores [node, distance] pair
	for node in nx.topological_sort(G):
		# pairs of dist,node for all incoming edges
		pairs = [(dist[v][0]+1,v) for v in G.pred[node]] 
		if pairs:
			dist[node] = max(pairs)
		else:
			dist[node] = (0, node)
	node,(length,_)  = max(dist.items(), key=lambda x:x[1])
	path = []
	while length > 0:
		path.append(node)
		length,node = dist[node]
	return list(reversed(path))

locus_branch_mapping_tmp=defaultdict(list)

with open(variants_file) as fp:
	for line in fp:
		if line[0]==">":
			var=line.split("/")
			bubbleID = int(var[0][9:])
			tmpvector = line.split(" ")
			branchID = int(var[1].split(" ")[0])
			tmp=list()
			for n in range(len(tmpvector)-2,0,-1):
				if tmpvector[n]!="VertexIDs:":
					tmp.append(tmpvector[n])
				else:
					break
			locus_branch_mapping_tmp[bubbleID].append(tmp)

#print(locus_branch_mapping_tmp)


bubble_mapping=defaultdict()
for key,value in locus_branch_mapping_tmp.items():
	if len(value) ==2:
		if len(list(value[0]))==1 and len(list(value[1]))==1:
			str1=" ".join(str(x) for x in value[0])
			str2 = " ".join(str(x) for x in value[1])
			bubble_mapping[str1] = str2
print(nx.weakly_connected_components(G))
tmp=[]
#print(nx.weakly_connected_component_subgraphs(G))
for c in nx.weakly_connected_component_subgraphs(G):
	if len(longest_path(c)) >0:
		for i in longest_path(c):
			if i in bubble_mapping:
				tmp.append(bubble_mapping[i])
			else:
				tmp.append(i)
#print(longest_path(G))
#print(tmp)
#paths = list(nx.shortest_simple_paths(G, 'DB9GZKS1:345:H5LNFBCXX:2:1104:10355:46030/2', 'DB9GZKS1:346:H5LCMBCXX:2:1114:10126:98756/2'))
#print(max(paths, key=len))
#print(paths)

with open(filename) as fp:
	for line in fp:
		var=line.split(	)
		if(var[0] == 'S'):
			if var[1] in tmp:
				output.write(line)
		else:
			if var[1] in tmp and var[3] in tmp:		
				if (tmp.index(var[1])-tmp.index(var[3]))==1 or (tmp.index(var[1])-tmp.index(var[3]))==-1:
					output.write(line)



