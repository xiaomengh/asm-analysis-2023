#!/usr/bin/env python

import sys
first = 0

def gt2af(gt):
    fields = [float(val) for val in gt.split(':')]
    if fields[0] == 0:
        return str(-1);
    else:
        return str(fields[2] / fields[0])

for line in sys.stdin:
    if first == 0:
        print(line.rstrip())
        first = 1
        continue

    tokens = line.split('\t')
    var = [tokens[0]]
    afs = [gt2af(gt) for gt in tokens[1:]]
    print("\t".join(var + afs))



