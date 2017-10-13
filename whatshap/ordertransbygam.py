rep_bubbles_filename = sys.argv[1]
gam_filename = sys.argv[2]
trans_filename = sys.argv[3]

#out = sys.argv[2]

d={}
count=1

rep_nodes = []
with open(rep_bubbles_filename) as fp:
	for line in fp:
		var=line.rstrip()
		rep_nodes.append(var)

canu_gam_order = []
with stream.open(str(gam_filename), "rb") as istream:
	for data in istream:
		g = vg_pb2.Alignment()
		g.ParseFromString(data)
		for i in range(0,len(g.path.mapping)):
			node = g.path.mapping[i].position.node_id
			if node not in rep_nodes:
				canu_gam_order.append(node)

bubbles_start = defaultdict()
with stream.open('assembly_graph.P.int.remn2n.X_100.chrXIII.trans' ,"rb") as istream:
	for data in istream:
		l = vg_pb2.SnarlTraversal()
		l.ParseFromString(data)
		if l.snarl.start.backward == True:
			start_node = l.snarl.end.node_id
		else:
			start_node = l.snarl.start.node_id
		if start_node not in rep_nodes:
			bubbles_start[start_node] = l


sorted_dict = OrderedDict()
index_map = {v: i for i, v in enumerate(canu_gam_order)}
for i in sorted(a.items(), key=lambda pair: index_map[pair[0]]):
	sorted_dict[i[0]] = i[1]


ostream = stream.open(str(sys.argv[4]), 'wb')

for k,v in sorted_dict.items():
	ostream.write(v)
