#!/usr/bin/python

import sys, os
from Bio import AlignIO

if len( sys.argv) < 4:
    print "USAGE:\n", sys.argv[0]," CMD_FILE  DELTA_STRENGTH  ANCHOR_TEMPLATE_FILE SIMILARITY_MATRIX\n"
    print "\tCMD_FILE: executable file containing the command calling AlignMe"
    print "\tNOTE:  the anchor file inside the command needs to be named 'anchors.txt'!\n"
    print "\tANCHOR_TEMPLATE_FILE: anchor definitions, SHOULD NOT be called 'anchors.txt' and in the same directory\n"
    print "\tSIMILARITY_MATRIX: e.g. BLOSUM matrix\n\n"
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

def ExtractScore( file):
    with open( file) as f:
        l = f.readline()
        c = l.split()
        if len(c) == 0:
            print "did not find last value!"
        return float( c[-1])


### READ PARAMETERS ###

cmd = sys.argv[1]
delta = float( sys.argv[2])
anchor_file = sys.argv[3]

matrix = {}
if len(sys.argv) > 4:
    similarity_file = sys.argv[4]
    with open( similarity_file) as f:
        header = f.readline().split()
        lines = f.readlines()
        for i in range( 0, len(lines)):
            cols = lines[i].split()
            for j in range( 1, len(cols)):
                matrix[ header[i] + header[j-1]] = float( cols[j])

anchors = []
### READ ANCHORS (only first two numbers in each line) ###
with open( anchor_file) as f:
    for l in f:
        cols = l.split()
        if len(cols) >= 2:
            anchors.append( [ cols[0], cols[1] ] )

### DEFINITIONS ###
gap_opening = 10
gap_extension = 1
strength = 0.0
out = "anchors.txt"
max_dead = 50
penalty = -10

nr_anchors = str( len( anchors))
shift_sum = -1         # define shift sum as negative for first while condition
prev = shift_sum
deadlock = 0

if out == anchor_file:
    print "anchor file should not be called: ", out
    exit()

w = open( "opt.txt", "w")
w.write( "# anchor-strength   nr-shifts   anchored_seq-similarity   whole-seq_seq-similarity   AlignMe_score\n")

### MAIN LOOP ###
# continue loop until anchors are all matched or nothing happens for max_dead
while shift_sum != 0 and deadlock < max_dead:
    
    # write anchor file for alignment with current strength
    f = open( out, "w")
    f.write( nr_anchors + "\n")
    for a in anchors:
        f.write( a[0] + " " + a[1] + " " + str( strength) + "\n")
    f.close()

    # perform alignment
    os.system( "./" + cmd)
    
    # read alignment using Biopython
    align = AlignIO.read("aligned_sequences.aln","clustal")
            
    # calc number of shifts and similarity score of anchored positions
    shift_sum = 0
    sim_score = 0

    # anchor loop
    for a in anchors:
        # helpers
        i = int(a[0]) - 1
        j = int(a[1]) - 1
        gapi = 0
        gapj = 0
        ci = 0
        cj = 0

        # counting gaps
        for c in align[0]:
            if c == "-":
                gapi += 1
            else:
                if ci == i:
                    break
                else:
                    ci += 1

        for c in align[1]:
            if c == "-":
                gapj += 1
            else:
                if cj == j:
                    break
                else:
                    cj += 1

        # sum of position plus number of gaps
        toti = i + gapi
        totj = j + gapj

        print i, ci, gapi, toti
        print j, cj, gapj, totj 

        # the shift score
        shift_sum += abs( toti - totj )

        # the similarity score
        if align[0][toti] != "-" and align[1][toti] != "-":
            sim_score += matrix[ align[0][toti] + align[1][toti]]
        else:
            sim_score += penalty
        if align[0][totj] != "-" and align[1][totj] != "-":
            sim_score += matrix[ align[0][totj] + align[1][totj]]
        else:
            sim_score += penalty
        # end of anchor loop

    # avoid infinite loop
    if shift_sum == prev:
        deadlock += 1
    else:
        prev = shift_sum

    # main output
    w.write( str(strength) + " " + str(shift_sum) + "  " + str( 0.5 * sim_score) + "  " + str( ScoreAlignment( align, matrix, gap_opening, gap_extension)) + "  " + str( ExtractScore( "aligned_sequences.aln")) + "\n")

    strength += delta
    # end of loop
w.close()


