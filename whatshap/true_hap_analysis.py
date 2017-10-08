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


d=defaultdict(list)
count=1


with stream.open('true_hap1.gam', "rb") as istream:
	for data in istream:
		g = vg_pb2.Alignment()
		g.ParseFromString(data)
		#if g.name != "S1_SK1_110":
		#	continue
		tmp=[]
		for i in range(0,len(g.path.mapping)):
			node = g.path.mapping[i].position.node_id
			d[node].append(g.name)
		#print(tmp)

count=0
max_val = 0
for k,v in d.items():
	for item, count in collections.Counter(v).items():
		if len(v) > 1:
			print(str(k) + "," +str(item) + "," + str(count))
		if len(set(v)) > 1:
			count+=1
			x = len(set(v))
			if x > max_val:
				max_val = x


print('max nodes with interchromosomal')
print(count)

print('max chromosomes that share any node')
print(max_val)





													


