__author__ = 'hudaiber'

import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/lib/BioPy/')
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/lib/BioPy/')
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')

import dm_tools as dt
import global_variables as gv
sys.path.append(gv.project_code_path)
import os
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree, fcluster, leaders, leaves_list
import matplotlib.pyplot as plt
from lib.utils import tools as t
import time


class Locus(object):

    def __init__(self, file_name):

        self.file_name = file_name

        _genes = dt.get_wgs_file(file_name)
        self.genes = _genes

        _forward = set()
        _reverse = set()

        for i in range(len(_genes)):
            _gene = _genes[i]
            for _cogid in _gene.cogid.split(','):
                _forward.update((_cogid,))
                if i == len(_genes)-1:
                    continue
                _next_gene = _genes[i+1]
                for _next_cogid in _next_gene.cogid.split(','):
                    _forward.update(("%s-%s"%(_cogid, _next_cogid),))

        _genes.sort(reverse=True)

        for i in range(len(_genes)):
            _gene = _genes[i]
            for _cogid in _gene.cogid.split(','):
                _reverse.update((_cogid,))
                if i == len(_genes)-1:
                    continue
                _next_gene = _genes[i+1]
                for _next_cogid in _next_gene.cogid.split(','):
                    _reverse.update(("%s-%s"%(_cogid, _next_cogid),))

        self.forward_set = _forward
        self.reverse_set = _reverse

    @staticmethod
    def calculate(first, second):

        score_intersection = sum([0.5 if '-' in term else 1 for term in first.intersection(second)])
        score_union        = sum([0.5 if '-' in term else 1 for term in first.union(second)])

        return score_intersection / score_union

    def score(self, other):

        ff = self.calculate(self.forward_set, other.forward_set)
        fr = self.calculate(self.forward_set, other.reverse_set)
        rf = self.calculate(self.reverse_set, other.forward_set)
        rr = self.calculate(self.reverse_set, other.reverse_set)

        return max(ff, fr, rf, rr)


def jackard_weighted_scores(loci, save_file=None):

    M = np.zeros((len(loci), len(loci)))
    no_of_loci = len(loci)

    print "Starting score calculations"
    tic = time.time()
    for i in range(no_of_loci):
        if i % 1000 == 0:
            print i, "time for last block: %f" % (time.time() - tic)
            tic = time.time()
        for j in range(i, no_of_loci):
            M[i, j] = loci[i].score(loci[j])
    print "Score calculations finished"

    if save_file:
        np.savez_compressed(os.path.join(save_file, 'jw_scores.npz'), data=M)

    return M

class Node(object):

    def __init__(self, other, parent=None):
        self.count = other.count
        self.dist = other.dist
        self.id = other.id
        self.is_leaf = other.is_leaf
        self.parent = parent
        self.left = None
        self.right = None


def clone_graph(root_node, parent=None):

    new_node = Node(root_node, parent)

    if new_node.is_leaf():
        return new_node
    else:
        new_node.left  = clone_graph(root_node.left,  new_node)
        new_node.right = clone_graph(root_node.right, new_node)
        return new_node


def get_leaves(root_node):

    if root_node.is_leaf():
        return [root_node.id]
    else:
        return get_leaves(root_node.left) + get_leaves(root_node.right)


def get_nodes(root_node):

    if root_node.is_leaf():
        return [root_node]
    else:
        return get_nodes(root_node.left) + [root_node] + get_nodes(root_node.right)


def break_down(root_node, size_limit=1000):

    if root_node.count<size_limit:
        return [root_node]
    else:
        return break_down(root_node.left, size_limit) + break_down(root_node.right, size_limit)


if __name__=='__main__':

    files_path = os.path.join(gv.project_data_path,'CRISPR/Cas1/genes_and_flanks/wgs/')
    pickle_path = os.path.join(gv.project_data_path, 'CRISPR/pickle/')

    print "Loading data"
    Z = np.load('/Users/hudaiber/Projects/NewSystems/data/CRISPR/pickle/upgma_distance_array.npz' ).items()[0][1]
    jwd = np.load('/Users/hudaiber/Projects/NewSystems/data/CRISPR/pickle/jw_distances.npz').items()[0][1]
    fname = os.path.join(pickle_path, 'loci.p.bz2')
    loci = t.load_compressed_pickle(fname)

    print "Starting "
    root = to_tree(Z)
    root = clone_graph(root)
    nodes = get_nodes(root)
    id2node = {node.id: node for node in nodes}

    leaf_ids = leaves_list(Z)

    cnt = 0
    i = 0
    threshold = 0.5
    clusters = []
    total_count = 0

    pool = []
    print "Starting merging"

    while True:
        cur_node = id2node[leaf_ids[i]]
        parent_dist = cur_node.parent.dist

        while parent_dist < threshold:
            cur_node = cur_node.parent
            parent_dist = cur_node.parent.dist

        cur_leaf_ids = get_leaves(cur_node)

        if len(cur_leaf_ids) > 1000:
            descendants = break_down(cur_node, size_limit=1000)
            for _n in descendants:
                pool.append([id for id in get_leaves(_n)])
        else:
            pool.append([id for id in cur_leaf_ids])
        total_count += cur_node.count
        # idx = np.array(cur_leaf_ids)
        # print len(idx), cur_node.id

        # M_cluster = M[np.ix_(idx, idx)]
        # max_i, max_j = np.unravel_index(M_cluster.argmax(), M_cluster.shape)

        i += len(cur_leaf_ids)

        if i >= len(leaf_ids)-1:
            print "Finished"
            print cnt, total_count
            break
        cnt += 1

    to_collapse = [l for l in pool if len(l)>1]
    singles = [l for l in pool if len(l)==1]
    print "To collapse:", len(to_collapse)
    print "Singles", len(singles)

    to_collapse = sorted(to_collapse, key=lambda x: len(x), reverse=True)

    out_fmt = "%s\t%d\t%s\n"

    with open('collapsed_loci.txt', 'w') as fout:
        for cur_list in to_collapse:

            cur_loci = sorted( [loci[id] for id in cur_list], key= lambda x: len(x.genes), reverse=True)
            file_names = [l.file_name.split('/')[-1] for l in cur_loci]
            fout.write(out_fmt%(file_names[0], len(file_names), " ".join(file_names)))

    with open('singleton_loci.txt', 'w') as fout:

        [fout.write(loci[id[0]].file_name.split('/')[-1]+"\t1\n") for id in singles]

    # for i, l in enumerate(to_collapse[:10]):
    #     print i, len(l)