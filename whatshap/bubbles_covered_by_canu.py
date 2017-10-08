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






nodes = set()
with stream.open('canu_new.contigs.gam', "rb") as istream:
	for data in istream:
		g = vg_pb2.Alignment()
		g.ParseFromString(data)
		for i in range(0,len(g.path.mapping)):
			node = g.path.mapping[i].position.node_id
			nodes.add(node)


bubbles_start = set()
covered_by_canu =set()
with stream.open('assembly_graph.P.int.remn2n.X_100.trans' ,"rb") as istream:
	for data in istream:
		l = vg_pb2.SnarlTraversal()
		l.ParseFromString(data)
		if l.snarl.start.backward == True:
			start_node = l.snarl.end.node_id
		else:
			start_node = l.snarl.start.node_id
		bubbles_start.add(start_node)
		if l.snarl.end.node_id in nodes or l.snarl.start.node_id in nodes:
			covered_by_canu.add(l.snarl.end.node_id)


print(len(bubbles_start))
print(len(covered_by_canu))


