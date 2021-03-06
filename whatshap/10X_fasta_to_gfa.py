#! /usr/bin/python

import Bio
import Bio.SeqIO
import sys
import argparse

parser = argparse.ArgumentParser(description='Convert Supernova assembler "raw" style FASTA output to GFA format')
parser.add_argument('file',help="The FASTA file to convert")
args = parser.parse_args()

handle = Bio.SeqIO.parse(args.file,"fasta")

vertex = dict() ##All the vertexes 

print("H	VN:Z:1.0") ##GFA header

i=0
for record in handle:
		s = record.description.find("edges=")
		d = record.description[s:].find(" ") + s
		EDGE_ID = record.description[s+6:d] ##Plus 6 removes edges=
		print("\t".join(["S",EDGE_ID,str(record.seq)]))
		
		s = record.description.find("left=")
		d = record.description[s:].find(" ") + s
		LEFT_VERTEX = record.description[s+5:d] ##Plus 5 removes left=
		
		s = record.description.find("right=")
		d = record.description[s:].find(" ") + s
		RIGHT_VERTEX = record.description[s+6:d] ##Plus 6 removes right=		

		try:
			vertex[LEFT_VERTEX]["out"].append(EDGE_ID)
		except KeyError:
			vertex[LEFT_VERTEX] = dict()
			vertex[LEFT_VERTEX]["in"] = []
			vertex[LEFT_VERTEX]["out"] = []
			vertex[LEFT_VERTEX]["out"].append(EDGE_ID)
		try:
			vertex[RIGHT_VERTEX]["in"].append(EDGE_ID)
		except KeyError:
			vertex[RIGHT_VERTEX] = dict()
			vertex[RIGHT_VERTEX]["in"] = []
			vertex[RIGHT_VERTEX]["out"] = []
			vertex[RIGHT_VERTEX]["in"].append(EDGE_ID)
handle.close()
sys.stderr.write('Done with FASTA file\n')
sys.stderr.write('Now writing links...\n')
for key in vertex:
	for i in vertex[key]["in"]:
		for j in vertex[key]["out"]:
			print("\t".join(["L",i,"+",j,"+","0M"]))
