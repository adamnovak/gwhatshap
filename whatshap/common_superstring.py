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

bubbles_start = set()
with stream.open('assembly_graph.P.int.remn2n.X_100.chrXIII.trans' ,"rb") as istream:
#with stream.open('assembly_graph.P.int.remn2n.X_100.trans' ,"rb") as istream:
	for data in istream:
		l = vg_pb2.SnarlTraversal()
		l.ParseFromString(data)
		if l.snarl.start.backward == True:
			start_node = l.snarl.end.node_id
		else:
			start_node = l.snarl.start.node_id
		bubbles_start.add(start_node)


multiplicity_bubbles = defaultdict(list)
read_details = defaultdict(list)
with stream.open('../out.new.chrXIII.gam', "rb") as istream:
	for data in istream:
		g = vg_pb2.Alignment()
		g.ParseFromString(data)
		tmp = []
		for i in range(0,len(g.path.mapping)):
			node = g.path.mapping[i].position.node_id
			if node in bubbles_start:
				multiplicity_bubbles[node].append(g.name)
			if node in tmp:
				if node is bubbles_start:
					bubbles_start.remove(node)
			tmp.append(node)
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
	#if max_val > 5 and l1_max > 1: # minimum read support and minimum multiplicty
	if l1_max > 1:
		#print(k)
		bubbles_start.remove(k)
		#repeaticity[k] = l1_max
		#repeaticity_read_support[k] = item_val
		#print(item_val, max_val, l1_max, k, v)
	if len(set(v)) > 38:
		#print(k)
		bubbles_start.remove(k)

#print(len(bubbles_start))

aln_paths = []
with stream.open('../out.new.chrXIII.gam', "rb") as istream:
	for data in istream:
		g = vg_pb2.Alignment()
		g.ParseFromString(data)
		tmp = []
		for i in range(0,len(g.path.mapping)):
			node = g.path.mapping[i].position.node_id
			if node in bubbles_start:
				tmp.append(node)
		aln_paths.append(tmp)

print(len(aln_paths))


def findOverlappingPair(str1, str2, str):
	max = -1
	len1 = len(str1)
	len2 = len(str2)

	for i in range(1, min(len1, len2)+1):
		if str1[len1-i:] == str2[:i]:
			if max < i:
				max = i
				str.clear()
				str.extend(str1 + str2[i:])

	for i in range(1, min(len1, len2)+1):
		if str1[:i] == str2[len2-i:]:
			if max < i:
				max = i
				str.clear()
				str.extend(str2 + str1[i:])
	print(max)
	return max

def findShortestSuperstring(arr):
	n_seqs = len(arr)
	while n_seqs != 1:
		max = -1
		l = 0
		r = 0
		result = []

		for i in range(0, n_seqs):
			for j in range(i+1, n_seqs):
				str = []
				res = findOverlappingPair(arr[i], arr[j], str)
				if max < res:
					max = res
					result = str
					l = i
					r = j
		n_seqs -= 1

		if max == -1:
			arr[0] += arr[n_seqs]
			print(arr[0])
		else:
			arr[l] = result
			arr[r] = arr[n_seqs]


	return arr[0]

fileout = open('common_superstring.txt', 'w')
fileout.write(findShortestSuperstring(aln_paths))
print(findShortestSuperstring(aln_paths))
		

	

