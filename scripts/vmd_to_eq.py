# Copyright (c) 2012 Thomas Spura <thomas.spura@gmail.com>
import sys
toA = 0.5291772108

if len(sys.argv) != 2:
    print "First argument must be vmd file."
    sys.exit(1)

for i,line in enumerate(open(sys.argv[1])):
    if i == 0:
        na = int(line)
        print na/3, na, 1
    elif i == 1:
        box = [float(item) for item in line.split()[1:]]
        print box[0]/toA, box[1]/toA, box[2]/toA
    else:
        box = [float(item) for item in line.split()[1:]]
        print box[0]/toA, box[1]/toA, box[2]/toA,
