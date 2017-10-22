import sys
import networkx as nx

G=nx.DiGraph()

filename = sys.argv[1]
out = sys.argv[2]


output= open(out, 'w')
with open(filename) as fp:
	for line in fp:
		var=line.split(	)
		if(var[0] == 'S'):
			G.add_node(var[1])
		else:
			G.add_edge(var[1],var[3])


#print(G.edges())
for x in list(nx.weakly_connected_component_subgraphs(G)):
	output.write(">")
	output.write(' '.join(nx.topological_sort(x)))

#print(x[0].nodes())

