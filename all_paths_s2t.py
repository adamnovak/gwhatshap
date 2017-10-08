import networkx as nx
def main():
	nodeToNodes = dict()
	nodeToNodes['A'] = ['B']
	nodeToNodes['B'] = ['E', 'C']
	nodeToNodes['C'] = ['E', 'D']
	nodeToNodes['D'] = ['A']

	print(str(getAllSimplePaths('A', 'E', nodeToNodes)))

	G = nx.DiGraph()
	for k,v in nodeToNodes.items():
		G.add_node(k)
		for v1 in v:
			G.add_node(v1)
			G.add_edge(k,v1)
	paths = nx.all_simple_paths(G, source='A', target='E', cutoff=1000)
	print(list(paths))
	for k,v in G.out_degree().items():
		if v ==1:
			print('hello')
			print(k)
	print(G.in_degree())


 
#
# Return all distinct simple paths from "originNode" to "targetNode".
# We are given the graph in the form of a adjacency list "nodeToNodes".
#
def getAllSimplePaths(originNode, targetNode, nodeToNodes):
	return helpGetAllSimplePaths(targetNode, [originNode], list(originNode), nodeToNodes, list())
 
#
# Return all distinct simple paths ending at "targetNode", continuing
# from "currentPath". "usedNodes" is useful so we can quickly skip
# nodes we have already added to "currentPath". When a new solution path
# is found, append it to "answerPaths" and return it.
#	
def helpGetAllSimplePaths(targetNode, currentPath, usedNodes, nodeToNodes, answerPaths):
	lastNode = currentPath[-1]
	if lastNode == targetNode:
		answerPaths.append(list(currentPath))
	else:
		for neighbor in nodeToNodes[lastNode]:
			if usedNodes.count(neighbor) <=1:
				currentPath.append(neighbor)
				usedNodes.append(neighbor)
				#if len(currentPath) > 5:
				#	return None
				helpGetAllSimplePaths(targetNode, currentPath, usedNodes, nodeToNodes, answerPaths)
				usedNodes.remove(neighbor)
				currentPath.pop()
	return answerPaths
 
if __name__ == '__main__':
	main()

