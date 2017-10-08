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
filename = sys.argv[1]
#out = sys.argv[2]

d={}
count=1
with open(filename) as fp:
	for line in fp:
		var=line.rstrip()
		edge = var+"_"+next(fp).rstrip()
		#print(edge)
		d[edge]=defaultdict() # periodicity in a read, read support


with stream.open('ouralns.SK1.gam', "rb") as istream:
	for data in istream:
		g = vg_pb2.Alignment()
		g.ParseFromString(data)
		#if g.name != "S1_SK1_110":
		#	continue
		tmp=[]
		for i in range(0,len(g.path.mapping)-1):
			edge1 = str(g.path.mapping[i].position.node_id)+"_"+str(g.path.mapping[i+1].position.node_id) # go over nodes in a mapping
			edge2 = str(g.path.mapping[i+1].position.node_id)+"_"+str(g.path.mapping[i].position.node_id) # go over nodes in a mapping
			#if edge1 in d or edge2 in d:
			#print(edge1)
			if edge1 in d:
				#print('in1')
				tmp.append(edge1)
			else:
				#print('in2')
				#print(g.name)
				tmp.append(edge2)
		#print(tmp)
		for item, count in collections.Counter(tmp).items():
			if count >=1:
				#print(str(item))
				if item in d:
					#print('yes')
					if count not in d[item]:
						d[item][count] = 0
					else:
						x = d[item][count]+1
						d[item][count] = x


													
count =0

for k,v in d.items():
	tmp = set()
	for k1,v1 in v.items():
		tmp.add(v1)

	if len(list(tmp))>1:
		print('hello')	
		if max(tmp) < 3:
			print(count)
			count = count +1
		
print(count)	

