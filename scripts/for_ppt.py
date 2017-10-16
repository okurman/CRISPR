#!/usr/bin/env python
__author__ = 'Sanjarbek Hudaiberdiev'

import os
import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')

import global_variables as gv
sys.path.append(gv.project_code_path)

import lib.utils.tools as t


# print "Bacteria"
# files = ['duplets_merged_within.p.bz2', 'triplets_merged_within.p.bz2',
#          'quadruplets_merged_within.p.bz2', 'pentaplets_merged_within.p.bz2']
# for f in files:
#     fpath = os.path.join(gv.project_data_path,'Bacteria/pickle/100000/',f)
#     # print fpath
#     kplets = t.load_compressed_pickle(fpath)
#     print f, len(kplets)


print "CRISPR"
files = ['duplets_seed.p.bz2', 'triplets_seed.p.bz2', 'quadruplets_seed.p.bz2', 'pentaplets_seed.p.bz2']
for f in files:
    fpath = os.path.join(gv.project_data_path,'CRISPR/pickle/cas1_2',f)
    kplets = t.load_compressed_pickle(fpath)
    print f, len(kplets)



