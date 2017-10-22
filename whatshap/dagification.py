import sys
import networkx as nx

G=nx.DiGraph()

filename = sys.argv[1]
out = sys.argv[2]

#[['DB9GZKS1:346:H5LCMBCXX:1:1108:9448:53577/1'], ['DB9GZKS1:345:H5LNFBCXX:2:2204:18574:80957/2'], ['DB9GZKS1:346:H5LCMBCXX:2:1204:15368:17140/2'], ['DB9GZKS1:345:H5LNFBCXX:1:2112:16983:19144/1'], ['DB9GZKS1:346:H5LCMBCXX:1:2208:17435:47435/1'], ['DB9GZKS1:346:H5LCMBCXX:1:2207:14553:15515/2'], ['DB9GZKS1:346:H5LCMBCXX:2:1112:19538:84361/2'], ['DB9GZKS1:345:H5LNFBCXX:2:1212:15375:83118/2'], ['DB9GZKS1:345:H5LNFBCXX:2:1212:10434:40185/2'], ['DB9GZKS1:346:H5LCMBCXX:2:1101:9874:98930/1'], ['DB9GZKS1:345:H5LNFBCXX:1:1103:8254:68288/1'], ['DB9GZKS1:346:H5LCMBCXX:2:2103:3765:71967/1']]
output= open(out, 'w')
with open(filename) as fp:
	for line in fp:
		var=line.split(	)
		if(var[0] == 'S'):
			G.add_node(var[1])
		else:
			G.add_edge(var[1],var[3])

for i in G.nodes_with_selfloops():
	G.remove_edge(i,i)

for i in nx.isolates(G):
	G.remove_node(i)
#print G.edges()
cycles = list(nx.simple_cycles(G))
print(cycles)

# greedy approach to remove cycles in graph
while(len(cycles) != 0):
	#print(cycles)
	tmp=[]
	for i in cycles:
		i.append(i[0])
		tmp1=set()
		for first, second in zip(i, i[1:]):
			tmp1.add((first,second))

		tmp.append(tmp1)
	#print(tmp)
	#print(set.intersection(*tmp))
	if len(list(set.intersection(*tmp)))>0:
		edges_to_remove = list(set.intersection(*tmp))
		ele = list(edges_to_remove[0])
		#print(ele)
		G.remove_edge(ele[0], ele[1])
	cycles = list(nx.simple_cycles(G))
	print(cycles)

with open(filename) as fp:
	for line in fp:
		var=line.split(	)
		if(var[0] == 'S'):
			#if var[1] in G.nodes():
			output.write(line)
		if(var[0] == 'L'):
			#print(tuple([var[1],var[3]]))
			if tuple([var[1],var[3]]) in G.edges():
				output.write(line)

