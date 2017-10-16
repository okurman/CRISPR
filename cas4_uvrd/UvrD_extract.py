#!/home/hudaiber/env2.7/bin/python
__author__ = 'hudaiber'

import os
import sys
if sys.platform == 'darwin':
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))
elif sys.platform == 'linux2':
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))

import global_variables as gv
sys.path.append(gv.project_code_path)

from lib.utils import hhpred as hp
from lib.utils import tools as t

cas4_profiles = {"cd09637", "cd09659", "cls000170", "COG1468", "COG4343", "pfam01930", "pfam06023"}

# if __name__=="__main__":

# 	hhrs_dir = '/home/hudaiber/Projects/NewSystems/data/UvrD/except_prok1603/COG0210/hhpred/'

# 	for root, dirs, files in os.walk(hhrs_dir):
# 		for f in files:
# 			if not f.endswith('.seq.hhr'):
# 				continue

# 			hits = hp.hhsearch_parse(os.path.join(hhrs_dir, f))

# 			profiles = [h.profile for h in hits]
# 			# print len(profiles), f, profiles
# 			if 'COG0210' in profiles:
# 				print len(profiles), f, profiles


if __name__=="__main__":

    work_dir            = '/home/hudaiber/Projects/NewSystems/data/UvrD/paralogs/'
    # prok1603_ccp        = '/home/hudaiber/Projects/NewSystems/data/UvrD/paralogs/Prok1603.ccp.csv'
    prok1603_ccp        = '/dev/shm/Prok1603.ccp.csv'
    prok1603_except_dir = '/panfs/pan1/patternquest/Projects/NewSystems/data/UvrD/except_prok1603/COG0210/hhpred/'

    pk_gi2pr = {}
    for l in open(prok1603_ccp):
        terms = l.strip().split(',')
        gi=terms[0]
        profile=terms[6]
        t.update_dictionary_set(pk_gi2pr, gi, profile)
        
    outfmt = "%s\t%s\t%s"
    print outfmt % ("#GI", "Genome", "Cas4 fusion")
    print "#--------------------------------------------"
    
    all_lines = open(os.path.join(work_dir, 'from_not_prok1603/uvrd_filtered.pty')).readlines() + \
                open(os.path.join(work_dir, 'from_prok1603/uvrd.pty')).readlines()
                
    for pty_line in open(os.path.join(work_dir, 'from_prok1603/uvrd.pty')).readlines():
        
        parts = pty_line.strip().split()
        genome = parts[3]
        gi = parts[-1]
        acc = parts[-2]
        
        fusion=0
        
        if pk_gi2pr[gi].intersection(cas4_profiles):
            fusion = 1
            
        print outfmt % (gi, genome, str(fusion))


    for pty_line in open(os.path.join(work_dir, 'from_not_prok1603/uvrd_filtered.pty')):

        found = False
        fusion = 0
        profiles = None

        parts = pty_line.strip().split()
        genome = parts[3]
        gi = parts[-1]
        acc = parts[-2]

        _f = os.path.join(prok1603_except_dir, "%s.seq.hhr" % acc.replace(".", "_"))

        if os.path.exists(_f):
            hits = hp.hhsearch_parse(_f)
            found = True
        else:
            _f = os.path.join(prok1603_except_dir, "gi_%s.seq.hhr" % gi)
            if os.path.exists(_f):
                hits = hp.hhsearch_parse(_f)
                found = True

        if not found:
            continue

        profiles = {hit.profile for hit in hits}

        if 'COG0210' not in profiles:
            continue

        if profiles.intersection(cas4_profiles):
            fusion=1

        print outfmt % (gi, genome, str(fusion))

