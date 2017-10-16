#!/usr/bin/env python
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
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree, leaves_list
import matplotlib.pyplot as plt
import scipy.spatial.distance as ssd
import time
from lib.utils import tools as t

sys.setrecursionlimit(100000)

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

    no_of_loci = len(loci)
    M = np.zeros((no_of_loci, no_of_loci))

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
        np.savez_compressed(save_file, data=M)

    return M


def generate_data():

    print "Loading loci"
    loci = [Locus(os.path.join(files_path, f)) for f in os.listdir(files_path)]
    loci = [locus for locus in loci if len(locus.genes) > 2 ]

    fname = os.path.join(pickle_path, 'loci.p.bz2')
    t.dump_compressed_pickle(fname, loci)
    # loci = t.load_compressed_pickle(fname)
    out_file = os.path.join(pickle_path, 'jw_scores.npz')
    jackard_weighted_scores(loci, out_file)


def build_linkage(M, distance_file=None, linkage_file=None):
    """
    Input:
    M: nxn score matrix
    """

    M += np.transpose(M)
    M = np.negative(np.log(M))
    np.fill_diagonal(M, 0)
    inf_idx = np.isinf(M)
    M[inf_idx] = 100  #the maximum non-inf value was 4.7.

    if distance_file:
        np.savez_compressed(distance_file, data=M)

    print("Clustering with distance array")
    M_array = ssd.squareform(M)
    Z = linkage(M_array, method='average')

    if linkage_file:
        print("Dumping result")
        np.savez_compressed(linkage_file, data=Z)

    return Z


def fancy_dendrogram(*args, **kwargs):

    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata


def plot_dendrogram(Z, report_path):

    print("Plotting")
    plt.figure(figsize=(100, 100))
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('loci')
    plt.ylabel('distance')

    fancy_dendrogram(
    Z,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=4.,  # font size for the x axis labels
    annotate_above=10,
    max_d=2
    )

    plt.savefig(os.path.join(report_path, 'dendrogram_distance_array.eps'), format='eps', dpi=900)
    plt.savefig(os.path.join(report_path, 'dendrogram_distance_array.pdf'), format='pdf')


def get_leaves(root_node):

    if root_node.is_leaf():
        return [root_node.id]
    else:
        return get_leaves(root_node.left) + get_leaves(root_node.right)


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

    if root_node.count < size_limit:
        return [root_node]
    else:
        return break_down(root_node.left, size_limit) + break_down(root_node.right, size_limit)


if __name__=='__main__':

    files_path = os.path.join( gv.project_data_path, 'CRISPR/datasets/crispr/wgs/')
    pickle_path = os.path.join(gv.project_data_path, 'CRISPR/pickle/crispr/')
    report_path = os.path.join(gv.project_data_path, 'CRISPR/reports/crispr/dendrograms/')
    if not os.path.exists(report_path):
        os.mkdir(report_path)

    # generate_data()

    score_file = os.path.join(pickle_path, 'jw_scores.npz')
    distance_file = os.path.join(pickle_path, 'jw_distances.npz')
    linkage_file = os.path.join(pickle_path, 'upgma_linkage.npz')
    locus_file = os.path.join(pickle_path, 'loci.p.bz2')

    loci = t.load_compressed_pickle(locus_file)
    # M = np.load(score_file).items()[0][1]

    print("Generating/Loading linkage file")
    Z = build_linkage(M, distance_file, linkage_file)
    Z = np.load(linkage_file).items()[0][1]
    print("Starting to form the clusters")

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

        i += len(cur_leaf_ids)

        if i >= len(leaf_ids)-1:
            print "Finished"
            print cnt, total_count
            break
        cnt += 1

    to_collapse = [l for l in pool if len(l) > 1]
    singles = [l for l in pool if len(l) == 1]
    print "To collapse:", len(to_collapse)
    print "Singles", len(singles)

    to_collapse = sorted(to_collapse, key=lambda x: len(x), reverse=True)

    out_fmt = "%d\t%d\t%s\n"

    with open(os.path.join(gv.project_data_path,'CRISPR/redundancy_elimination/','crispr_collapsed_loci.txt'), 'w') as fout:
        for cnt, cur_list in enumerate(to_collapse):
            cur_loci = sorted([loci[id] for id in cur_list], key= lambda x: len(x.genes), reverse=True)
            file_names = [l.file_name.split('/')[-1] for l in cur_loci]
            fout.write(out_fmt%(cnt+1, len(file_names), ",".join(file_names)))

    with open(os.path.join(gv.project_data_path,'CRISPR/redundancy_elimination/','crispr_collapsed_loci_rep.txt'), 'w') as fout:
        for cnt, cur_list in enumerate(to_collapse):
            cur_loci = sorted([loci[id] for id in cur_list], key= lambda x: len(x.genes), reverse=True)
            file_names = [l.file_name.split('/')[-1] for l in cur_loci]
            fout.write(out_fmt%(cnt+1, len(file_names), cur_loci[0].file_name.split('/')[-1]))

    with open(os.path.join(gv.project_data_path,'CRISPR/redundancy_elimination/','crispr_singleton_loci.txt'), 'w') as fout:

        [fout.write(loci[id[0]].file_name.split('/')[-1]+"\n") for id in singles]