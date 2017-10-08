import pyfaidx
from xopen import xopen
import stream
import logging
from . import vg_pb2
from collections import Counter
import networkx as nx
from collections import defaultdict
from .core import ReadSet, Read
from functools import reduce
import operator as op


from contextlib import ExitStack
from .vcf import VcfReader, PhasedVcfWriter
from . import __version__
from .core import ReadSet, readselection, Pedigree, PedigreeDPTable, NumericSampleIds, PhredGenotypeLikelihoods
from .graph import ComponentFinder
from .pedigree import (PedReader, mendelian_conflict, recombination_cost_map,
                       load_genetic_map, uniform_recombination_map, find_recombination)
from .bam import BamIndexingError, SampleNotFoundError
from .timer import StageTimer
from .variants import ReadSetReader, ReadSetError

__author__ = "Shilpa Garg, Tobias Marschall"

logger = logging.getLogger(__name__)

from heapq import heappush, heappop
from itertools import count

import networkx as nx

"""
consider only top-k paths from complex bubbles.
"""

def k_shortest_paths(G, source, target, k=1, weight='weight'):
	"""Returns the k-shortest paths from source to target in a weighted graph G.
	"""
	if source == target:
		return ([0], [[source]]) 
	   
	length, path = nx.single_source_dijkstra(G, source, target, weight=weight)
	if target not in length:
		raise nx.NetworkXNoPath("node %s not reachable from %s" % (source, target))
		
	lengths = [length[target]]
	paths = [path[target]]
	c = count()		
	B = []						
	G_original = G.copy()	
	
	for i in range(1, k):
		for j in range(len(paths[-1]) - 1):			
			spur_node = paths[-1][j]
			root_path = paths[-1][:j + 1]
			
			edges_removed = []
			for c_path in paths:
				if len(c_path) > j and root_path == c_path[:j + 1]:
					u = c_path[j]
					v = c_path[j + 1]
					if G.has_edge(u, v):
						edge_attr = G.edge[u][v]
						G.remove_edge(u, v)
						edges_removed.append((u, v, edge_attr))
			
			for n in range(len(root_path) - 1):
				node = root_path[n]
				# out-edges
				for u, v, edge_attr in G.copy().edges_iter(node, data=True):
					G.remove_edge(u, v)
					edges_removed.append((u, v, edge_attr))
				
				if G.is_directed():
					# in-edges
					for u, v, edge_attr in G.in_edges_iter(node, data=True):
						G.remove_edge(u, v)
						edges_removed.append((u, v, edge_attr))
			
			spur_path_length, spur_path = nx.single_source_dijkstra(G, spur_node, target, weight=weight)			
			if target in spur_path and spur_path[target]:
				total_path = root_path[:-1] + spur_path[target]
				total_path_length = get_path_length(G_original, root_path, weight) + spur_path_length[target]				
				heappush(B, (total_path_length, next(c), total_path))
				
			for e in edges_removed:
				u, v, edge_attr = e
				G.add_edge(u, v, edge_attr)
					   
		if B:
			(l, _, p) = heappop(B)		
			lengths.append(l)
			paths.append(p)
		else:
			break
	
	return (lengths, paths)

def get_path_length(G, path, weight='weight'):
	length = 0
	if len(path) > 1:
		for i in range(len(path) - 1):
			u = path[i]
			v = path[i + 1]
			
			length += G.edge[u][v].get(weight, 1)
	
	return length 

"""
output the possible allele-pairs for a bubble.
"""
def ncr(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, range(n, n-r, -1))
    denom = reduce(op.mul, range(1, r+1))
    return numer//denom

"""
Input: Phase variants from Locus file and aligned reads from GAM file.

It creates an association between het variants and read alignments. 

Output: The optimal partitioning is written to standard output.
"""
def vg_reader(variants_file, gam_file, ordering_file):
	"""
	input: reads locus and GAM file output from vg.
	output: sorted readset for core DP.
	assumptions: 
	1. locus file consists of linear ordering of simple bubbles only and hence sorted.
	2. paths in the locus should be covered by atleast one pacbio read.
	2. GAM file is sorted and restricted to locus file.
	3. files consists of all connected components.
	4. add variant only when it identifies the branch uniquely.
	"""
	#create an ordering list.
	ordering_list=[]
	with open(ordering_file) as fp:
		tmp_len=len(ordering_list)
		for line in fp:
			if line[0]==">":
				var=line.split(" ")
				for i in var:
					ordering_list.append(tmp_len+int(i))
	# create a dictionary of branches for each locus based on variantsfile.
	locus_branch_mapping_tmp=defaultdict(list)
	alleles_per_pos= defaultdict()

	with open(variants_file) as fp:
		for line in fp:
			if line[0]==">":
				var=line.split("/")
				bubbleID = int(var[0][9:])
				tmpvector = line.split(" ")
				branchID = int(var[1].split(" ")[0])
				tmp=set()
				for n in range(len(tmpvector)-2,0,-1):
					if tmpvector[n]!="VertexIDs:":
						tmp.add(tmpvector[n])
					else:
						break
				locus_branch_mapping_tmp[bubbleID].append(tmp)

	# create a list nodes from bubbles, one node from each bubble.	
	tmp_bubbles=[]
	for key, value in locus_branch_mapping_tmp.items():
		tmp_bubbles.append(list(value[0])[0])
	# restrict odering list to only bubbles nodes.
	odering_bubbles = sorted(set(ordering_list) & set(tmp_bubbles), key = ordering_list.index)
	# sorted bubbles
	locus_branch_mapping1=defaultdict(list)
	for key, value in locus_branch_mapping_tmp.items():
		new_key = odering_bubbles.index(list(value[0])[0])
		locus_branch_mapping1[new_key] = value


	# top-k paths in bubbles
	G=nx.Graph()
	count=0
	for key, value in locus_branch_mapping1.items():
		if len(value)>=5: # hardcoded k- value here.
			for paths in value:
				for i in range(0,len(paths)-1):
					if i==0:
						counter=counter-1
						G.add_node(str(counter))
						G.add_edge(list(paths)[counter], list(paths)[i], weight=0)
					if i==len(paths)-2:
						counter=counter-1
						G.add_node(str(counter))
						G.add_node(list(paths)[len(paths)-1])
						G.add_edge(list(paths)[len(paths)-1], list(paths)[counter], weight=0)					
					G.add_node(list(paths)[i])
					G.add_edge(list(paths)[i], list(paths)[i+1], weight=0)
	
	# read gam file to count read support for each egde.
	with stream.open(str(gam_file), "rb") as istream:
		for data in istream:
			aln = vg_pb2.Alignment()
			aln.ParseFromString(data)
			node_list=[]
			for i in aln.path.mapping: # go over the mapping in a read
				#if i.position.is_reverse==True:
				#	break
				node = i.position.node_id # go over nodes in a mapping
				node_list.append(node)
			for i in range(0, len(node_list)-1):
				G[u][v]['weight']+=1

	counter =0
	for key, value in locus_branch_mapping1.items():
		if len(value)< 5:
			locus_branch_mapping[key]=value
		else:
			counter1=counter-1
			counter= counter-2
			result = k_shortest_paths(G, counter1, counter2, 5) # harcoded here.
			counter=counter-2
			for i in rasult[1]:
				del i[0]
				del i[-1]
			locus_branch_mapping[key] = result[1] # check once of list of lists is ok, or list of sets is required.

	for i in locus_branch_mapping:
		alleles_per_pos[i] = len(locus_branch_mapping[i])
				
	# create a dictionary of branches for each locus based on locus file.
	#locus_branch_mapping=defaultdict()
	#locus_count=0
	#alleles_per_pos= defaultdict()
	#for l in locus_list:
		#per_locus_position=[]
		#print(l.allele)
		#if len(l.allele) > 1:
			#for paths in l.allele:
				#tmp=set()
				#if len(paths.mapping) > 1: # donot consider homozygous locus, only with multiple paths.
					#for mappings in paths.mapping:
						#tmp.add(mappings.position.node_id)
					#per_locus_position.append(tmp)
			#alleles_per_pos[locus_count] = len(per_locus_position)
			#locus_branch_mapping[locus_count]=per_locus_position
			#locus_count=locus_count+1
			
	alleles_per_pos= defaultdict()
	for k,v in locus_branch_mapping.items():
		alleles_per_pos[k]=len(v)		


	# key is the values in locus_branch_mapping and value is triplet(locus, branch, alleles to be flipped)
	# if a simple bubbles consists of complicated path, then there would atleast one node that uniquely determine the branch.
	reverse_mapping= defaultdict(list)
	for k,v in locus_branch_mapping.items():
		if len(v) > 1: # more than one branch
			for i,b in enumerate(v):
				if len(b) > 0:
					for p,j in enumerate(b): # reverse_mapping for every node.
						reverse_mapping[j].append([k,i, len(v)])


	# extract reads from GAM file associated with the locus and create a sorted readset.
	readset=ReadSet()  
	for g in gam_list:
		# hard-coded source id, mapping quality and other values.
		read=Read(g.name, 0, 0, 0) # create read for each read alignment
		prev_tmp=[]
		prev_locus= -1
		for i in g.path.mapping: # go over the mapping in a read
			node = i.position.node_id # go over nodes in a mapping
			if node in reverse_mapping and reverse_mapping[node][0][0]!=-1: # handle start and sink node.
				qualities = [10]* reverse_mapping[node][0][2]
				tmp=[]
				node_inf= [tuple(i[0:2]) for i in reverse_mapping[node]]
				for i in node_inf:
					tmp.append(i)
				interset_tmp= list(set(tmp).intersection(set(prev_tmp)))
				if len(tmp)==1 and len(prev_tmp)==0: # simple bubble with only two paths.
					qualities[tmp[0][1]] = 0
					read.add_variant(tmp[0][0], tmp[0][1], qualities) # if any new locus of branch is encountered, enter a variant in read.
					locus = tmp[0][0]
				elif len(prev_tmp) > 0 and len(set(tmp).intersection(set(prev_tmp)))==1: # for complicated bubbles, but with Top-N paths. combination of some nodes uniquely determine branch.
					qualities[interset_tmp[0][1]] = 0
					read.add_variant(interset_tmp[0][0], interset_tmp[0][1], qualities)
					locus= interset_tmp[0][0]
				elif len(prev_tmp) > 0 and len(tmp) ==1: # complicated bubbles, but one node can uniquely identify the branch.
					qualities[tmp[0][1]] = 0
					read.add_variant(tmp[0][0], tmp[0][1], qualities)
					locus = tmp[0][0]

				if prev_locus!=locus:
					prev_tmp = []
				else:
					for i in tmp:
						prev_tmp.append(i)
		if len(read) >= 2:
			readset.add(read)
	for read in readset:
			read.sort()
	readset.sort()
	return readset, alleles_per_pos

def run_phaseg(gam_file):
	"""
	Run WhatsHap.

	gam_file -- path to GAM file
	locus_file -- path to input variants
	"""
	timers = StageTimer()
	recombrate=1.26
	all_heterozygous = False
	distrust_genotypes = True
	timers.start('overall')
	#logger.info("This is WhatsHap %s running under Python %s", __version__, platform.python_version())
	with ExitStack() as stack:
		logger.info('Using uniform recombination rate of %g cM/Mb.', recombrate)
		all_reads, alleles_per_pos = vg_reader(locus_file, gam_file)
		print(all_reads)
		selected_indices = readselection(all_reads, 15)
		selected_reads = all_reads.subset(selected_indices)
		accessible_positions = sorted(selected_reads.get_positions())
		pedigree = Pedigree(NumericSampleIds())
		# compute the number of alleles at each position.
		alleles_per_accessible_pos =[]
		genotype_likelihoods = []
		for pos in accessible_positions:
			if pos in alleles_per_pos:
				n_alleles = alleles_per_pos[pos]  
				possible_genotypes = n_alleles +  ncr(n_alleles, 2)
				genotype_likelihoods.append(None if all_heterozygous else PhredGenotypeLikelihoods([0]* possible_genotypes))
		# random input of genotypes, since distrust_genotypes is always ON.
		pedigree.add_individual('individual0', [0]* len(accessible_positions), genotype_likelihoods)
		recombination_costs = uniform_recombination_map(recombrate, accessible_positions)
		# Finally, run phasing algorithm
		print(selected_reads)
		dp_table = PedigreeDPTable(selected_reads, recombination_costs, pedigree, distrust_genotypes, accessible_positions)
		superreads_list, transmission_vector = dp_table.get_super_reads()
		cost = dp_table.get_optimal_cost()
		print(superreads_list[0])
		print(cost)
		read_partitions = dp_table.get_optimal_partitioning()
		print(read_partitions)


def add_arguments(parser):
	arg = parser.add_argument
	# Positional arguments
	#arg('locus_file', metavar='LOCUS', help='variants in LOCUS file to phase')
	arg('gam_file', metavar='PHASEINPUT',
	    help='read alignments in GAM file ')
	#TODO: add odering file.
	#TODO: add k parameter

def main(args):
	run_phaseg(**vars(args))

