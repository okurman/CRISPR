__author__ = 'hudaiber'

import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')
import global_variables as gv
import os


class Cluster(object):

    def __init__(self, id, size, files):
        self.id = id
        self.size = size
        if isinstance(files, tuple):
            self.files = files
        else:
            self.files = set(files)


def get_singleton_loci(dataset):

    if dataset == 1:
        fname = os.path.join(gv.project_data_path, 'CRISPR/redundancy_elimination/cas1_1_singleton_loci.txt')
    elif dataset == 2:
        fname = os.path.join(gv.project_data_path, 'CRISPR/redundancy_elimination/cas1_2_singleton_loci.txt')
    elif dataset == 3:
        fname = os.path.join(gv.project_data_path, 'CRISPR/redundancy_elimination/crispr_singleton_loci.txt')
    return set([l.strip() for l in open(fname).readlines()])


def get_clustered_loci(dataset):

    if dataset == 1:
        fname = os.path.join(gv.project_data_path, 'CRISPR/redundancy_elimination/cas1_1_collapsed_loci.txt')
    elif dataset == 2:
        fname = os.path.join(gv.project_data_path, 'CRISPR/redundancy_elimination/cas1_2_collapsed_loci.txt')
    elif dataset == 3:
        fname = os.path.join(gv.project_data_path, 'CRISPR/redundancy_elimination/crispr_collapsed_loci.txt')

    clusters = []
    with open(fname) as f:
        for l in f:
            _id, _size, _files = l.split('\t')
            _files = set(_files.strip().split(','))
            clusters.append(Cluster(_id, _size, _files))

    return clusters