import collections, sys
from itertools import groupby
import string
import networkx as nx


def parseGFA(f):
    rev = dict(('+-','-+','MM'))
    nodes = {}
    G = collections.defaultdict(dict)
    for line in f:
        if len(line) == 0:
            pass
        c = line[0]
        if c == 'H':
            pass
        elif c == 'S':
            l = line.split()
            nodes[l[2]]=l[1]
        elif c == 'L':
            l = line.split()
            i,j = l[1],l[3]
            o = (l[2],l[4],l[5])
            if j not in G[i]:
                G[i][j] = o
            if i not in G[j]:
                G[j][i] = rev[o[1]],rev[o[0]]

    for i in nodes:
        if i not in G:
            d = G[i]
    #print(nodes)
    return G,nodes

def create_graph(f):
	G=nx.DiGraph()
	for line in f:
		if len(line) == 0:
			pass
		c = line[0]
		if c == "H":
			pass
		else c == "S":
			line.split("\t")
			G.add_node()


if __name__ == "__main__":
    gfafile = open(sys.argv[1])
    
    G,nodes = parseGFA(gfafile)
    print(nx.connected_components(G))
    #fasta_iter(sys.argv[2],nodes,G)
