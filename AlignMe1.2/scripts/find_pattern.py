#!/usr/bin/python

import sys, re

if len( sys.argv) < 4:
    print "\nUSAGE:", sys.argv[0], " FASTA  PATTERN  OUTFILE \n\n"
    exit()


infile = sys.argv[1]
pattern = sys.argv[2]
outfile = sys.argv[3]

out = open( outfile, "w")

seq = ""
with open( infile) as f:
    for l in f:
        if l[0] != ">":
            seq = seq + l.rstrip()

for m in re.finditer( pattern, seq):
         out.write( str( m.start() + 1) + " " + str( m.end()) + "  " + seq[m.start():m.end()] + "\n")

out.close()
