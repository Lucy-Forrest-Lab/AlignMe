#!/usr/bin/python

import sys

w = open( sys.argv[2], "w")


with open( sys.argv[1]) as f:
    for l in f:
        l = l.replace( ".","-")
        w.write(l)

