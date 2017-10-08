import sys
import stream
import logging
import vg_pb2
from collections import Counter
from collections import defaultdict
import collections
from collections import OrderedDict, namedtuple
from collections import defaultdict
import networkx as nx


bubble_nodes = set()
bubble_start = defaultdict()
bubble_end = set()
with stream.open('assembly_graph.P.int.remn2n.X_100.chrXIII.trans' ,"rb") as istream:
	for data in istream:
		l = vg_pb2.SnarlTraversal()
		l.ParseFromString(data)
		start = l.snarl.start.node_id
		end = l.snarl.end.node_id
		if l.snarl.start.backward == True:
			bubble_start[end] = start
		else:
			bubble_start[start] = end

		for i in range(0,len(l.visits)):
			bubble_nodes.add(l.visits[i].node_id)

print(len(bubble_nodes))


G = nx.DiGraph()
with stream.open('assembly_graph.P.int.remn2n.X_100.chrXIII.vg', "rb") as istream:
	for data in istream:
		l = vg_pb2.Graph()
		l.ParseFromString(data)
		for j in range(len(l.edge)):
			from_edge = getattr(l.edge[j], "from")
			to_edge = l.edge[j].to
			if from_edge in bubble_start: # collapse bubbles
				G.add_node(from_edge)
				G.add_node(bubble_start[from_edge])
				G.add_edge(from_edge, bubble_start[from_edge])
			elif from_edge not in bubble_nodes and to_edge not in bubble_nodes:
				G.add_edge(from_edge, to_edge) 
				G.add_node(from_edge)
				G.add_node(to_edge)

file_out = open('collapsed_bubbles.gfa', 'w')

for i in G.nodes():
	file_out.write("S" + '\t' + str(i) + "\t" + "A" + "\n")

for j in G.edges():
	file_out.write("L" + '\t' + str(j[0]) + "\t" + '+' + "\t" + str(j[1]) + "\t" + '+' + "\t" + '0M' + "\n")


def is_subset(chains,tmp):
	for i in range(len(chains)):
		if set(tmp).issubset(chains[i]):
				return True
	return False

tmp = [] # temporary storage for chains
chains = [] # stores all found chains
nodes = G.nodes() # save the labels of nodes
prev_node =  None # Holds the node which other nodes will be checked if they are a subset of.
for i in range(len(nodes)): # each node is considered as a superset and others will be compared to it
	tmp.append(nodes[i]) # append this node as the head of the chain
	prev_node = nodes[i] # store it to compare other nodes to it
	for j in range (i+1,len(nodes)): # looping the rest of the nodes
		if set(G.neighbors(nodes[j])).issubset(set(((G.neighbors(prev_node))))) : # if it is a subset
			tmp.append(nodes[j]) # append to the chain
			prev_node = nodes[j] # store to check if the rest of the nodes are subset of this one

	if not (is_subset(chains,tmp)): # After finishing the chain check it is not a subset of an already created chain
		#print(chains)
		chains.append(list(reversed(tmp))) # add it to the chains
	tmp =[] # reset
	prev_node = None
print(len(chains))

count = 0
for i in chains:
	if len(i) > 200:
		count+=1

print("I am greater than 200")
print(count)
#print(chains)







			

