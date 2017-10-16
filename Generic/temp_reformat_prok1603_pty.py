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
from lib.utils import tools as t
import glob

root_path = '/home/hudaiber/data/Prok1603/pty/'

print "Loading profile2genes"
profile2gene = t.map_cdd_defense2gene_name()

print "Loading ccp"
map_file = '/panfs/pan1.be-md.ncbi.nlm.nih.gov/prokdata/db_tmp/Prok1603/Annotation/Prok1603.ccp.csv'
gi2cdd = t.map_gid2cdd(map_file)

print "Loading CDD"
cdd2gene = t.map_cdd_defense2gene_name()

cnt=0

cas4=0
cas4_uvrd=0
uvrd=0

for dir in os.listdir(root_path):

    cnt += 1

    org_genes = ""

    for file in glob.glob(os.path.join(root_path, dir, "*.pty")):

        file = os.path.basename(file)
        new_file = os.path.join(root_path, dir, ".".join(file.split('.')[:-1]) + ".pty2")

        with open(os.path.join(root_path, dir, file)) as inf:
            with open(new_file, 'w') as outf:

                for l in inf:

                    if l.startswith("#"):
                        outf.write(l)
                        continue

                    l = l.strip()
                    gi = l.split()[-1]

                    profiles = gi2cdd[gi] if gi in gi2cdd else ""
                    gene_names = " ".join(cdd2gene[profile] if profile in cdd2gene else profile
                                          for profile in gi2cdd[gi].split(" ")
                                          ) if gi in gi2cdd else ""

                    if not profiles:
                        outf.write(l + "\n")
                    else:
                        outf.write(l + "\t" + profiles  + "\t" + gene_names + "\n")

    # for file in os.listdir(os.path.join(root_path, dir)):
    #
    #     with open(os.path.join(root_path, dir, file)) as inf:
    #
    #         for l in inf:
    #             l = l.strip()
    #
    #             gi = l.split()[-1]
    #
    #             gene_names = " ".join(cdd2gene[profile] if profile in cdd2gene else profile
    #                                   for profile in gi2cdd[gi].split()
    #                                   ) if gi in gi2cdd else ""
    #
    #             org_genes += " "+gene_names
    #
    # org_genes = org_genes.lower()
    #
    # if 'cas4' in org_genes and 'uvrd' in org_genes:
    #
    #     cas4_uvrd += 1
    #     continue
    #
    # if 'cas4' in org_genes and 'uvrd' not in org_genes:
    #
    #     cas4 += 1
    #     continue
    #
    # if 'cas4' not in org_genes and 'uvrd' in org_genes:
    #
    #     uvrd += 1


    # if cnt % 100 == 0:
    #     print cnt, cas4, cas4_uvrd, uvrd

# print cnt, cas4, cas4_uvrd, uvrd