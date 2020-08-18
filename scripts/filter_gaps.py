#!/usr/bin/python

import sys, re

if len( sys.argv) < 3:
    print "\nUSAGE:", sys.argv[0], " FASTA  OUTFILE \n\n"
    exit()


infile = sys.argv[1]
outfile = sys.argv[2]

out = open( outfile, "w")

with open( infile) as f:
    for l in f:
        if l[0] != ">":
            l = l.replace( "-","")
            l = l.replace( ".","")
        out.write( l)


out.close()
