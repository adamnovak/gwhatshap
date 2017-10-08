import sys
import stream
import vg_pb2
from collections import Counter
from collections import defaultdict
import collections
from collections import OrderedDict, namedtuple

# assumption ... all S' and before L's

out = sys.argv[1]
f = open(out,'w')

bubble_to_remove = set()
bubbles_dict_trans = defaultdict(int)
with stream.open('component54.trans', "rb") as istream:
	for data in istream:
		l = vg_pb2.SnarlTraversal()
		l.ParseFromString(data)
		tmp = str(l.snarl.start.node_id)
		bubbles_dict_trans[tmp]=0




with stream.open('out.new.chrI.gam', "rb") as istream:
	for data in istream:
		g = vg_pb2.Alignment()
		g.ParseFromString(data)
		for i in range(0,len(g.path.mapping)):
			node = str(g.path.mapping[i].position.node_id)
			if node in bubbles_dict_trans:
				bubbles_dict_trans[node]+=1

for k,v in bubbles_dict_trans.items():
	f.write(str(k)+" "+ str(v)+"\n")




													
	

