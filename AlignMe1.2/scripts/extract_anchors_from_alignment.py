#!/usr/bin/python

from Bio import AlignIO
from optparse import OptionParser
from os import path

parser = OptionParser()

parser.add_option( "-f", "--alignement_file", dest="alignment_file", metavar="FILE",
                   help="file containing the alignment,                             \nNOTE: this is the ONLY required flag!                     others are optional")

parser.add_option( "-o", "--output_file", metavar="FILE", default="extracted.txt",
                   help="file for output of selected anchors")

parser.add_option( "-a", "--nr_allowed_gaps", default=0, type="int", dest="allowed_gaps", 
                   help="number of gaps that are allowed in a given window for it to be considered")

parser.add_option( "-p", "--penalty", default=0, type="float", dest="penalty",
                   help="penalty for a gap within window")

parser.add_option( "-w", "--window_size", default=11, type="int", 
                   help="size of the window for which averaged scales values are summed, NOTE: even numbered values will lead to odd numbered windows anyway")

parser.add_option( "-t", "--threshold", default=0, type="float", 
                   help="")

parser.add_option( "-s", "--scale", default="../scales/HWvH.txt", dest="scale_file", metavar="FILE",
                   help="file containing (hydrophobicity) scale that is used for scoring windows")

parser.set_defaults( filter="--below_threshold")

parser.add_option( "--below_threshold", action="store_const",const="below", dest="filter",
                   help="select all windows with total scores below threshold for anchors (default behaviour)")

parser.add_option( "--above_threshold", action="store_const",const="above", dest="filter",
                   help="select all windows with total scores above threshold for anchors")

parser.set_defaults( merge="--average")

parser.add_option( "--average",action="store_const", const="average", dest="merge",
                   help="averages both values in aligned positions for score calculation (defaul behaviour)")

parser.add_option( "--min", action="store_const", const="min", dest="merge",
                   help="uses minimum value of aligned positions for score calculation")

parser.add_option( "--max", action="store_const", const="max", dest="merge",
                   help="uses maximum value of aligned positions for score calculation")

(options, args) = parser.parse_args()

print parser.parse_args()

if not options.alignment_file:
    print "\n\nERROR: no alignment file given! bailing out. \n\n"
    parser.print_help()
    exit()

if not path.isfile( options.scale_file):
    print "\n\nERROR: scale file not found, consider using '-s'\n\n"
    parser.print_help()
    exit()

half_size = int( options.window_size / 2)

# read alignment
align = AlignIO.read( options.alignment_file ,"clustal")

# read scale
scale = {}
with open( options.scale_file ) as f:
    for l in f:
        cols = l.split()
        scale[ cols[0]] = float( cols[1])

# open output
w = open( options.output_file ,"w")
pro = open( "profile.txt", "w")

# array that contains anchor information as '1', '0' else
array = ['0'] * len( align[0])

# main loop over alignment
for i in range( 0, len( align[0])):

    # sum over window
    sum = 0
    for j in range( max( 0, i - half_size), min( len(align[0]), i + half_size + 1)):
        a = align[0][j]
        b = align[1][j]
        if a != "-" and b != "-" and a != "." and b != ".":
            if options.merge == "--average":
                sum += 0.5 * ( scale[a] + scale[b])
            elif options.merge == "--min":
                sum += min( scale[a], scale[b])
            elif options.merge == "--max":
                sum += max( scale[a], scale[b])
        else:
            sum += options.penalty   # gap penalty
    pro.write( str( sum) + "\n")

    # set array values according to sum, collecting anchored alignment positions
    if options.filter == "--below_threshold":
        if sum < options.threshold:
            for j in range( max( 0, i - half_size), min( len(align[0]), i + half_size + 1)):
                if align[0][j] != "-" and align[1][j] != "-":
                    array[j] = 1   # ANCHOR
    elif options.filter == "--above_threshold":
        if sum > options.threshold:
            for j in range( max( 0, i - half_size), min( len(align[0]), i + half_size + 1)):
                if align[0][j] != "-" and align[1][j] != "-":
                    array[j] = 1   # ANCHOR

pro.close()

# translate alignment positions into sequence positions and output 
counter = 0
for val in array:
    if val == 1:
        seqA = align[0][:counter+1].format( "fasta")
        seqA = seqA[seqA.find("\n")+1:]
#        print seqA
        posA = counter - seqA.count("-") + 1

        seqB = align[1][:counter+1].format( "fasta")
        seqB = seqB[ seqB.find("\n")+1:]
#        print seqB
        posB = counter - seqB.count("-") + 1
 
#        print counter, posA, posB
        w.write( str( posA) + " " + str( posB) + "  0 \n")
    counter += 1

w.close()
