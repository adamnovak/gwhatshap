import re
import sys

#should nicely separate CIGAR entries
cigar_pat = re.compile(r"\d+[MIDNSHP=X]{1}")

def cigar2end( left,cigar ):
  """Return right-most position of aligned read."""
  #store info about each CIGAR category
  counts={ "M":0, #M 0 alignment match (can be a sequence match or mismatch)
           "I":0, #I 1 insertion to the reference
           "D":0, #D 2 deletion from the reference
           "N":0, #N 3 skipped region from the reference
           "S":0, #S 4 soft clipping (clipped sequences present in SEQ)
           "H":0, #H 5 hard clipping (clipped sequences NOT present in SEQ)
           "P":0, #P 6 padding (silent deletion from padded reference)
           "=":0, #= 7 sequence match
           "X":0, #X 8 sequence mismatch
        }
  #split cigar entries
  m=[]
  for centry in cigar_pat.findall(cigar):
    ccount  = int(centry[:-1])
    csymbol = centry[-1]
    counts[csymbol] = ccount
    m.append(counts["M"])
    #print(counts["M"])
  #get number of aligned 'reference' bases
  aligned = counts["M"] + counts["D"] + counts["N"] #+ counts["="] + counts["X"]
  right   = left + aligned
  return max(m)


filename = sys.argv[1]
minlist=[]

with open(filename) as fp:
  for line in fp:
    var=line.rstrip()
    if var!='*':
      x = cigar2end(0,var)
      minlist.append(x)

print(min(minlist))


