__author__ = 'hudaiber'

import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/lib/BioPy/')
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/lib/BioPy/')
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')

from lib.utils import tools as t

cdd_file = '/Users/hudaiber/data/CDD/all_Prok1402.ccp.csv'

gnm2weight = t.map_genome2weight()

profile2count = {}
profile2weight = {}

missing = []

for l in open(cdd_file):
    terms = l.split(',')
    org = terms[1]
    profile = terms[6]
    if org in gnm2weight:
        t.update_dictionary(profile2count, profile, 1)
        t.update_dictionary(profile2weight, profile, gnm2weight[org])
    else:
        missing.append(org)

print "Missing from weights:"
for gnm in set(missing):
    print gnm

print "Finished scanning"
with open('/Users/hudaiber/data/CDD/profile2weight.tab','w') as outf:
    outf.write("#Profile\tweight\tcount\n")
    for k,v in profile2weight.items():
        outf.write("%s\t%f\t%s\n"%(k,v, profile2count[k]))
print "Done"