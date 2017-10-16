#!/usr/bin/env python

fname='/Users/hudaiber/Projects/NewSystems/data/CRISPR/Cas1/genes_and_flanks/cnt_files.txt'
K=5


cnts = [int(l.split()[0]) for l in open(fname).readlines()]
from math import factorial as f

def count_f(n,r):
	return f(n)/(f(n-r)*f(r))

total_original = 0
total_splitted = 0

for K in [5,4,3,2]:
    print "K=",K
    tsum = 0
    for n in cnts:
        tsum += count_f(n,K) if n >= K else 0
    total_original += tsum
    print 'Original size:', tsum

    tsum = 0
    for n in cnts:
        if n < K:
            continue
        tsum += count_f(n,K) if n < 20 else 2*count_f(n/2, K)
    total_splitted += tsum
    print 'Splitted size:', tsum

print "Total original,", total_original
print "Total splitted,", total_splitted
print total_original/float(total_splitted)