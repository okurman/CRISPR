#!/usr/bin/env python
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

import numpy as np
import matplotlib.pyplot as plt

desc_file_fname = os.path.join(gv.project_data_path, 'cas4/baiticity/Clusters.descr')

id2gene = {}

for l in open(desc_file_fname):
    parts = l.split('\t')
    id = parts[0]
    genes = parts[2]
    id2gene[id] = genes

baits, cas_hits, others = [], [], []

cas4_hits_fname = os.path.join(gv.project_data_path, 'cas4/baiticity/cas4_ratio_sorted.tsv')

all_hits = open(cas4_hits_fname).readlines()

for l in all_hits:

    parts = l.split()
    id = parts[0]

    genes = id2gene[id]

    if 'cas4' in genes:
        baits.append(l)
    elif 'cas' in genes:
        cas_hits.append(l)
    else:
        others.append(l)

save_file = os.path.join(gv.project_data_path, 'cas4/baiticity/cas_hits.tsv')
with open(save_file, 'w') as fout:
    for l in cas_hits:
        id, in_loci, everywhere = l.strip().split('\t')
        try:
            icity = float(in_loci)/float(everywhere)

        except:
            continue
        fout.write("%s\t%s\t%s\t%f\n" % (id, in_loci, everywhere, icity))

save_file = os.path.join(gv.project_data_path, 'cas4/baiticity/bait_hits.tsv')
with open(save_file, 'w') as fout:
    for l in baits:
        id, in_loci, everywhere = l.strip().split('\t')
        try:
            icity = float(in_loci)/float(everywhere)
        except:
            continue
        fout.write("%s\t%s\t%s\t%f\n" % (id, in_loci, everywhere, icity))

save_file = os.path.join(gv.project_data_path, 'cas4/baiticity/other_hits.tsv')
with open(save_file, 'w') as fout:
    for l in others:
        id, in_loci, everywhere = l.strip().split('\t')
        try:
            icity = float(in_loci)/float(everywhere)
        except:
            continue
        fout.write("%s\t%s\t%s\t%f\n" % (id, in_loci, everywhere, icity))

plt.figure(figsize=(20,10))

ax = plt.gca()

fnames = ['other_hits.tsv', 'cas_hits.tsv', 'bait_hits.tsv']
labels = ['others', 'cas (~cas4)', 'cas4']
colors = ['blue', 'green', 'red']

for i in range(3):

    icity, in_loci = [], []

    fname = os.path.join(gv.project_data_path, 'cas4/baiticity/', fnames[i])

    for l in open(fname):
        parts = l.split()
        _icity = float(parts[3])

        if _icity > 1:
            continue
        
        in_loci.append(int(parts[1]))
        icity.append(_icity)

    icity = np.asarray(icity)
    in_loci = np.asarray(in_loci)

    plt.scatter(np.log10(in_loci), icity, color=colors[i], label=labels[i])



# icity, in_loci = [], []

# for l in open('hits_cas.tsv'):
#     parts = l.split()
#     in_loci.append(int(parts[1]))
#     icity.append(float(parts[3]))

# icity = np.asarray(icity)
# in_loci = np.asarray(in_loci)

# plt.scatter(np.log10(in_loci), icity, color='green', label='cas (~cas4)')

# icity, in_loci = [], []

# for l in open('hits_bait.tsv'):
#     parts = l.split()
#     in_loci.append(int(parts[1]))
#     icity.append(float(parts[3]))

# icity = np.asarray(icity)
# in_loci = np.asarray(in_loci)

# plt.scatter(np.log10(in_loci), icity, color='red', label='cas4')

plt.plot([0,3], [1.0, 1.0], color='black')
plt.plot([0,3], [0, 0], color='black')
plt.plot([0,3], [0.5, 0.5], color='black', linestyle='--')
plt.legend()
plt.xlim([0,3])

plt.xlabel("Occurrence in loci", fontsize=16, fontweight='bold')
plt.ylabel("Cas4icity", fontsize=16, fontweight='bold')
plt.title("Cas4icity plot", fontsize=16, fontweight='bold')

locs, labels = plt.xticks()
plt.xticks(locs, map(lambda x: "%d" % np.power(10,x), locs))

plt.savefig(os.path.join(gv.project_data_path, 'cas4/baiticity/baiticity.png'))