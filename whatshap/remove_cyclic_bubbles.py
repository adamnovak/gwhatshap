import sys
import stream
import vg_pb2
from collections import Counter
from collections import defaultdict
import collections
from collections import OrderedDict, namedtuple

# assumption ... all S' and before L's
filename = sys.argv[1]
f = open(filename, 'w')
#out = sys.argv[2]

duplicated_nodes =set()
count=1
with stream.open('out.new.gam', "rb") as istream:
	for data in istream:
		g = vg_pb2.Alignment()
		g.ParseFromString(data)
		read_count =0
		tmp =[]
		for i in range(0,len(g.path.mapping)):
			node = g.path.mapping[i].position.node_id
			if node in tmp:
				duplicated_nodes.add(node)
			tmp.append(node)

for i in list(duplicated_nodes):
	f.write(str(i)+ ",")
			

print(len(list(duplicated_nodes)))
													
	

