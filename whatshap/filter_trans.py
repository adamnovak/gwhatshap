import sys
import stream
import logging
import vg_pb2
from collections import Counter
from collections import defaultdict
import collections
from collections import OrderedDict, namedtuple
from collections import defaultdict

nodes_list = set()
for line in open('alignment_8_tig00000547_pilon.gam.nodes'):
	nodes_list.add(int(line.rstrip()))


ostream = stream.open('alignment_8_tig00000547_pilon.gam.gfa.vg.ordered.filtered.trans', 'wb')
with stream.open('alignment_8_tig00000547_pilon.gam.gfa.vg.ordered.trans' ,"rb") as istream:
	for data in istream:
		l = vg_pb2.SnarlTraversal()
		l.ParseFromString(data)
		if l.snarl.start.node_id in nodes_list:
			ostream.write(l)

ostream.close()
			

