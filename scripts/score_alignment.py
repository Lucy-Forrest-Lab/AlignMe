#!/usr/bin/python

import sys, os
from Bio import AlignIO

if len( sys.argv) < 3:
    print "USAGE:\n", sys.argv[0]," ALIGNMENT_FILE SIMILARITY_MATRIX\n"
    exit()

def ScoreAlignment( align, matrix, gap_opening, gap_extension):
    length = len( align[0])
    sum = 0
    prev_was_gap = False
    for i in range( 0, length):
        if align[0][i] == "-" or align[0][i] == "." or align[1][i] == "-" or align[1][i] == ".":
            if prev_was_gap:
                sum -= gap_extension
            else:
                sum -= gap_opening
                prev_was_gap = True
        else:
            sum += matrix[ align[0][i] + align[1][i]]
            prev_was_gap = False
    return sum

### READ PARAMETERS ###

matrix = {}
similarity_file = sys.argv[2]
with open( similarity_file) as f:
    header = f.readline().split()
    lines = f.readlines()
    for i in range( 0, len(lines)):
        cols = lines[i].split()
        for j in range( 1, len(cols)):
            matrix[ header[i] + header[j-1]] = float( cols[j])

### DEFINITIONS ###
gap_opening = 10
gap_extension = 1

# read alignment using Biopython
align = AlignIO.read( sys.argv[1],"clustal")
            

print "score:", ScoreAlignment( align, matrix, gap_opening, gap_extension)


