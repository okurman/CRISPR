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

import dm_tools as dt
import global_variables as gv
sys.path.append(gv.project_code_path)

import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree, leaves_list
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.spatial.distance as ssd
import time
from lib.utils import tools as t
sys.setrecursionlimit(100000)
import pandas as pd
import xlsxwriter as x
import lib.utils.reporting as r
import shutil as sh
from operator import itemgetter

from lib.db.generic import map_profiles_id2code_code2def
(_, profile_code2def) = map_profiles_id2code_code2def('cas')




def generate_data_feature(feature_definition_file, method, loci=None, save_file=None):

    if not method:
        print "Provide a method for score calculations"
        sys.exit()

    if not loci:
        print "Loading loci"
        loci = [Locus(os.path.join(files_path, f)) for f in os.listdir(files_path)]
        loci = [locus for locus in loci if len(locus.genes) > 2 ]

    feature_labels, feature_weights = [], []

    with open(feature_definition_file) as inf:
        for l in inf:
            terms = l.strip().split()
            feature_labels.append(terms[0])
            feature_weights.append(terms[2])

    # set feature weights to presence/absebce values instead
    feature_weights = [1 if w>0 else 0 for w in feature_weights]

    for locus in loci:
        locus.set_features(feature_labels, feature_weights)

    no_of_loci = len(loci)
    M = np.zeros((no_of_loci, no_of_loci))

    print "Starting score calculations"
    for i in range(no_of_loci):
        for j in range(i+1, no_of_loci):
            dot_score = 1.0*np.dot(loci[i].feature_weights, loci[j].feature_weights)

            if dot_score == 0:
                M[i, j] = 0
                continue

            if method=='dot':
                M[i, j] = dot_score
            elif method=='union':
                union_length = np.sum(np.logical_or(loci[i].feature_weights, loci[j].feature_weights))
                M[i, j] = dot_score/union_length
            elif method == 'min':
                min_length = min(sum(loci[i].feature_weights), sum(loci[j].feature_weights))
                M[i, j] = dot_score/min_length
            else:
                print "Invalid method for dot score calculation!"
                sys.exit()

    print "Score calculations finished"

    if save_file:
        np.savez_compressed(save_file, data=M)

    if method=='dot':
        M = M / np.max(M)

    return M, feature_labels, feature_weights





def filter_feature_profiles(thr_occurrence, thr_crispricity):

    M = pd.read_csv(os.path.join( gv.project_data_path, 'cas1402/crispricity.tab'),
                    sep="\t")

    M = M[(M['Occurrence in CRISPR loci'] >= thr_occurrence) & (M['Crispricity'] >= thr_crispricity)]
    out_file = os.path.join( gv.project_data_path, 'cas1402/profiles_%d_%.2f.tab'%(thr_occurrence, thr_crispricity))
    M.to_csv(out_file, sep="\t", index=False, header=None)


def classify_by_scores(M, threshold, loci):

    M_array = ssd.squareform(M)
    # main linkage structure for upgma
    # print "Building linkage"
    Z = linkage(M_array, method='average')

    # Z = np.load(linkage_file).items()[0][1]
    # print "plotting dendogram"
    # plot_dendrogram(Z, report_path)

    root = to_tree(Z)
    root = clone_graph(root)

    nodes = get_nodes(root)
    id2node = {node.id: node for node in nodes}

    leaf_ids = leaves_list(Z)

    cnt = 0
    i = 0
    total_count = 1

    pool = []
    # print "Starting merging"

    while True:
        cur_node = id2node[leaf_ids[i]]
        parent_dist = cur_node.parent.dist

        while parent_dist < threshold:
            cur_node = cur_node.parent
            parent_dist = cur_node.parent.dist

        cur_leaf_ids = get_leaves(cur_node)

        pool.append([id for id in cur_leaf_ids])

        total_count += cur_node.count

        i += len(cur_leaf_ids)

        if i >= len(leaf_ids)-1:
            break
        cnt += 1

    to_collapse = [l for l in pool if len(l) > 1]
    singles = [l[0] for l in pool if len(l) == 1]

    to_collapse = sorted(to_collapse, key=lambda x: sum(gnm2weight[loci[i].organism] for i in x), reverse=True)

    sum_errors = []
    entropies = []
    weights = []
    to_collapse_retval = []

    cluster_ind = 0

    for cluster in to_collapse:

        cluster_ind += 1
        type2cnt = {}
        type2wgt = {}

        cluster_files = [loci[id].file_name.split('/')[-1] for id in cluster]

        cluster_weight = 0

        for _f in cluster_files:

            file_weight = gnm2weight[file2org[_f]]
            cluster_weight += file_weight

            if _f not in file2crispr_type:
                t.update_dictionary(type2cnt, "NA", 1)
                t.update_dictionary(type2wgt, "NA", file_weight)
                continue
            for _type in file2crispr_type[_f]:
                t.update_dictionary(type2cnt, _type, 1)
                t.update_dictionary(type2wgt, _type, file_weight)

        _weights = np.array(type2wgt.values(), dtype=np.float)

        sum_errors.append(np.sum(np.square(_weights - np.mean(_weights))))

        _weights /= np.sum(_weights)
        entropy = -1 * np.sum(_weights * np.log(_weights))
        entropies.append(entropy)
        weights.append(cluster_weight)

        to_collapse_retval.append((cluster, type2cnt, type2wgt, entropy))

    sum_errors = np.average(sum_errors)

    entropies = np.array(entropies)
    weights = np.array(weights)

    average_entropy = np.sum(entropies * weights) / np.sum(weights)
    sum_errors = np.sum(sum_errors * weights) / np.sum(weights)

    return singles, to_collapse_retval, sum_errors, average_entropy


def plot_block(block):

    _fname = block[0].strip()

    thresholds = []
    singles = []
    clusters = []
    mean_errors = []
    entropies = []

    for l in block[1:]:

        terms = l.strip().split('\t')

        thresholds.append(terms[0])
        singles.append(terms[1])
        clusters.append(terms[2])
        mean_errors.append(terms[3])
        entropies.append(terms[4])

    thresholds = np.asarray(thresholds,dtype=np.float)
    singles = np.asarray(singles, dtype=np.int)
    clusters = np.asarray(clusters, dtype=np.int)
    mean_errors = np.asarray(mean_errors, dtype=np.float)
    entropies = 1000*np.asarray(entropies, dtype=np.float)

    plt.plot(thresholds, singles)
    plt.plot(thresholds, clusters)
    plt.plot(thresholds, mean_errors)
    plt.plot(thresholds, entropies)

    _cnt = _fname.split('_')[1]
    _crispricity = _fname.split()[0].split('_')[2][:-4]
    _profiles = _fname.split()[1]

    plt.title("Occurrence:%s, CRISPRicity:%s, profiles: %s"%(_cnt, _crispricity, _profiles))
    plt.grid(True)
    plt.legend(['Singletons', 'Clusters', r'$ \left< \sum_{i=1}^N Error^2 \right> $',r'$ \left<\textit{I}\right> x(10^3) $'],loc='upper left')
    plt.xticks(thresholds, [str(t) for t in thresholds], rotation='vertical')
    plt.xlabel('Clustering thresholds')


def plot_results(data_file_name = 'results.txt', image_file_name='results.png'):

    plt.figure(figsize=(50,20))
    plt.rc('text', usetex=True)

    font = {'family': 'serif',
            'weight': 'bold',
            'size': 22}

    plt.rc('font', **font)


    block = []
    i = 1

    for l in open(data_file_name).readlines():
        if not block:
            block.append(l)
            continue

        if l.startswith('#'):
            plt.subplot(2,3,i)
            plt.tight_layout()
            plot_block(block)
            i += 1
            block = [l]
            continue

        block.append(l)

    plt.subplot(2,3,i)
    plot_block(block)

    plt.savefig(image_file_name)


def process_dot_product(thresholds, method, loci=[]):

    results_data_file_name = os.path.join(gv.project_data_path,'cas1402/results_%s.txt'%method)
    results_image_file_name = os.path.join(gv.project_data_path,'cas1402/clustering_%s.png'%method)

    # outf = open(results_data_file_name, 'w')

    for thr_occurrence, thr_crispricity in thresholds:

        _fname = 'profiles_%d_%.2f.tab'%(thr_occurrence, thr_crispricity)

        feature_definition_file = os.path.join(gv.project_data_path, 'cas1402/profiles/', _fname)
        print feature_definition_file

        if not os.path.exists(feature_definition_file):
            raise FileNotFoundError("File doesn't exist: %s" % feature_definition_file)

        M, feature_labels, _ = generate_data_feature(feature_definition_file, method=method, loci=loci)

        M += np.transpose(M)

        M = -1 * np.log(M)
        M[np.diag_indices_from(M)] = 0
        M[np.where(M==np.inf)] = 50

        # For dendrograms
        # M_array = ssd.squareform(M)
        # Z = linkage(M_array, method='average')
        # dendrogram_file = os.path.join(gv.project_data_path, 'cas1402/dendrograms/dendrogram_%s_%d_%.2f.pdf'%(method, thr_occurrence, thr_crispricity))
        # print "Printing dendrogram:", dendrogram_file
        # plot_dendrogram(Z, dendrogram_file)
        # End dendrograms

        # save_file = os.path.join(pickle_path, 'cosine_distances.npz')
        # np.savez_compressed(save_file, data=M)
        # M = np.load(save_file).items()[0][1]

        outf.write("#%s\t%d\n"%(_fname, len(feature_labels)))

        for threshold in [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3,4,5,6,7,8,9,10]:
            # print threshold
            singles, clusters, sum_errors, entropies = classify_by_scores(M, threshold, loci)
            outf.write("%f\t%d\t%d\t%f\t%f\n" % (threshold, len(singles), len(clusters), sum_errors, entropies))
            outf.flush()

    outf.close()
    plot_results(results_data_file_name, results_image_file_name)

    print "Results saved:", results_data_file_name
    print "Results saved:", results_image_file_name


def process_jackard_distances(loci):

    # print "Calculating Jackard scores"
    # M = jackard_weighted_scores(loci)
    #
    # save_file = os.path.join(gv.project_data_path, 'cas1402/pickle/jw_scores.npz')
    #
    # print "Saving JW scores"
    # np.savez_compressed(save_file, data=M)
    #
    # # M = np.load(save_file).items()[0][1]
    #
    # M += np.transpose(M)
    # M = np.negative(np.log(M))
    # np.fill_diagonal(M, 0)
    # inf_idx = np.isinf(M)
    # M[inf_idx] = 100
    #
    result_file = os.path.join(gv.project_data_path,'cas1402/results_jw.txt')
    #
    # outf = open(result_file,'w')
    #
    # for threshold in [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3,4,5,6,7,8,9,10]:
    # # for threshold in [1.5, 1.75]:
    #     singles, clusters, sum_errors, entropies = classify_by_scores(M, threshold, loci)
    #     outf.write("%f\t%d\t%d\t%f\t%f\n" % (threshold, len(singles), len(clusters), sum_errors, entropies))
    #     print threshold
    #     print "Clusters",len(clusters)
    #     print "Singles",len(singles)
    #     outf.flush()
    #
    # outf.close()

    lines = open(result_file).readlines()
    plt.figure(figsize=(30,20))
    thresholds = []
    singles = []
    clusters = []
    mean_errors = []
    entropies = []


    for l in lines:
        terms = l.strip().split('\t')
        thresholds.append(terms[0])
        singles.append(terms[1])
        clusters.append(terms[2])
        mean_errors.append(terms[3])
        entropies.append(terms[4])

    thresholds = np.asarray(thresholds,dtype=np.float)
    singles = np.asarray(singles, dtype=np.int)
    clusters = np.asarray(clusters, dtype=np.int)
    mean_errors = np.asarray(mean_errors, dtype=np.float)/100
    entropies= np.asarray(entropies, dtype=np.float) * 1000

    plt.plot(thresholds, singles)
    plt.plot(thresholds, clusters)
    plt.plot(thresholds, mean_errors)
    plt.plot(thresholds, entropies)

    plt.title("Jackard scores")
    plt.grid(True)
    plt.legend(['Singletons','Clusters','SSE', 'I(x1K)'], loc='upper left')
    plt.xticks(thresholds,[str(t) for t in thresholds], rotation='vertical')
    plt.xlabel('Clustering thresholds')

    font = {'family': 'serif',
            'weight': 'bold',
            'size': 22}

    plt.rc('font', **font)

    # plt.ylabel('Respective values')

    plt.savefig(os.path.join(gv.project_data_path, 'cas1402/clustering_jw.png'))


def report_feature_clusters_jw(loci):

    M = generate_data_jackard()
    M += np.transpose(M)
    M = -1 * np.log(M/np.max(M))
    M[np.diag_indices_from(M)] = 0
    M[np.where(M==np.inf)] = 50

    reports_dir_base = os.path.join(gv.project_data_path, 'cas1402/reports/')

    for threshold in [3, 4, 6]:

        repors_dir = reports_dir_base + 'jw_%.2f'%threshold
        print repors_dir
        if os.path.exists(repors_dir):
            sh.rmtree(repors_dir)
        os.mkdir(repors_dir)
        
        singles, cluster_packs, sum_errors = classify_by_scores(M, threshold, loci)
        generate_cluster_reports(singles, cluster_packs, loci, repors_dir, [], [])


def report_feature_clusters_dot(loci, thresholds_pack, method):

    thr_occ, thr_crisp, cluster_thresholds = thresholds_pack
    feature_definition_file = os.path.join(gv.project_data_path, 'cas1402/profiles/profiles_%d_%.2f.tab' % (thr_occ, thr_crisp))

    M, feature_labels, feature_weights = generate_data_feature(feature_definition_file, method, loci=loci)
    M += np.transpose(M)
    print "Max score:",np.max(M)
    M = -1 * np.log(M)
    M[np.diag_indices_from(M)] = 0
    M[np.where(M==np.inf)] = 50

    reports_dir_base = os.path.join(gv.project_data_path, 'cas1402/reports/')

    for threshold in cluster_thresholds:

        repors_dir = reports_dir_base + 'dot_%s_%d_%.2f_%.2f'%(method, thr_occ, thr_crisp, threshold)
        print "Thresholds:", thr_occ, thr_crisp, threshold
        print repors_dir
        if os.path.exists(repors_dir):
            sh.rmtree(repors_dir)
        os.mkdir(repors_dir)

        singles, cluster_packs, sum_errors, entropies = classify_by_scores(M, threshold, loci)

        _local_thresholds_pack = (thr_occ, thr_crisp, threshold)

        generate_cluster_reports(singles, cluster_packs, loci, repors_dir, feature_labels, feature_weights, method, _local_thresholds_pack)


def singleton_reports(loci, thresholds_pack, method):

    thr_occ, thr_crisp, cluster_threshold = thresholds_pack
    feature_definition_file = os.path.join(gv.project_data_path, 'cas1402/profiles/profiles_%d_%.2f.tab' % (thr_occ, thr_crisp))

    M, feature_labels, feature_weights = generate_data_feature(feature_definition_file, method, loci=loci)
    M += np.transpose(M)
    M = -1 * np.log(M)
    M[np.diag_indices_from(M)] = 0
    M[np.where(M==np.inf)] = 50

    reports_dir_base = os.path.join(gv.project_data_path, 'cas1402/reports/')

    repors_dir = reports_dir_base + 'dot_singletons_%s_%d_%.2f_%.2f'%(method, thr_occ, thr_crisp, cluster_threshold)
    print "Thresholds:", thr_occ, thr_crisp, cluster_threshold
    print repors_dir
    if os.path.exists(repors_dir):
        sh.rmtree(repors_dir)
    os.mkdir(repors_dir)

    singles, cluster_packs, sum_errors, entropies = classify_by_scores(M, cluster_threshold, loci)

    weights = np.array([gnm2weight[locus.organism] for locus in loci])
    M = M*weights

    avg_dist = np.sum(M,1)/(2*np.sum(weights))
    avg_dist = avg_dist[singles]
    loci = [loci[i] for i in singles]
    sort_idx = np.argsort(avg_dist)

    single_groups = [sort_idx[:20], sort_idx[-20:]]

    generate_singleton_reports(single_groups, loci, repors_dir, feature_labels)


def profile_statistics_after_clustering(loci, thresholds_pack, method):

    thr_occ, thr_crisp, cluster_threshold = thresholds_pack
    feature_definition_file = os.path.join(gv.project_data_path, 'cas1402/profiles/profiles_%d_%.2f.tab' % (thr_occ, thr_crisp))

    out_dir = os.path.join(gv.project_data_path, 'cas1402/profile_statistics/')

    M, feature_labels, feature_weights = generate_data_feature(feature_definition_file, method, loci=loci)
    M += np.transpose(M)
    M = -1 * np.log(M)
    M[np.diag_indices_from(M)] = 0
    M[np.where(M==np.inf)] = 50

    singles, cluster_packs, sum_errors, entropies = classify_by_scores(M, cluster_threshold, loci)

    singleton_counts = np.zeros(len(feature_labels))
    cluster_counts = np.zeros(len(feature_labels))

    singletons_weight = 0
    clusters_weight = 0

    for _i in singles:
        locus = loci[_i]
        _weight = gnm2weight[locus.organism]
        singleton_counts += locus.feature_weights * _weight
        singletons_weight += _weight

    for cluster_pack in cluster_packs:
        (cluster, _, _, _) = cluster_pack
        for _i in cluster:
            locus = loci[_i]
            _weight = gnm2weight[locus.organism]
            cluster_counts += locus.feature_weights * _weight
            clusters_weight += _weight

    plt.figure()

    cluster_counts /= (clusters_weight + singletons_weight)
    singleton_counts /= (clusters_weight + singletons_weight)

    sort_idx = np.argsort(cluster_counts)[::-1]
    xs = np.array(range(len(cluster_counts)))
    plt.bar(xs, cluster_counts[sort_idx], color='b')
    plt.bar(xs, singleton_counts[sort_idx], color='r')

    plt.savefig(os.path.join(out_dir, "profiles_count.png"))


def profile_statistics_after_clustering_2(loci, thresholds_pack, method):

    thr_occ, thr_crisp, cluster_threshold = thresholds_pack
    feature_definition_file = os.path.join(gv.project_data_path, 'cas1402/profiles/profiles_%d_%.2f.tab' % (thr_occ, thr_crisp))

    out_dir = os.path.join(gv.project_data_path, 'cas1402/profile_statistics/')

    M, feature_labels, feature_weights = generate_data_feature(feature_definition_file, method, loci=loci)
    M += np.transpose(M)
    M = -1 * np.log(M)
    M[np.diag_indices_from(M)] = 0
    M[np.where(M==np.inf)] = 50

    singles, cluster_packs, sum_errors, entropies = classify_by_scores(M, cluster_threshold, loci)

    weights = np.array([gnm2weight[locus.organism] for locus in loci])
    M = M * weights

    avg_dist = np.sum(M, 1) / (2 * np.sum(weights))
    avg_dist = avg_dist[singles]

    singleton_counts = np.zeros(len(feature_labels))
    cluster_counts = np.zeros(len(feature_labels))

    singletons_weight = 0
    clusters_weight = 0

    for _i in singles:
        locus = loci[_i]
        _weight = gnm2weight[locus.organism]
        singleton_counts += locus.feature_weights * _weight
        singletons_weight += _weight

    for cluster_pack in cluster_packs:
        (cluster, _, _, _) = cluster_pack
        for _i in cluster:
            locus = loci[_i]
            _weight = gnm2weight[locus.organism]
            cluster_counts += locus.feature_weights * _weight
            clusters_weight += _weight

    plt.figure()

    cluster_counts /= (clusters_weight + singletons_weight)
    singleton_counts /= (clusters_weight + singletons_weight)

    sort_idx = np.argsort(cluster_counts)[::-1]
    xs = np.array(range(len(cluster_counts)))
    plt.bar(xs, cluster_counts[sort_idx], color='b')
    plt.bar(xs, singleton_counts[sort_idx], color='r')

    plt.savefig(os.path.join(out_dir, "profiles_count.png"))


def generate_singleton_reports(single_groups, loci, reports_dir, feature_labels):


    print "Generating report files"
    ind = 0

    for single_group in single_groups:

        ind += 1
        xls_file_name = os.path.join(reports_dir, '%d.xls' % ind)

        group_loci = sorted([loci[_i] for _i in single_group])

        params = {}

        params['xls_file_name']         = xls_file_name
        params['loci']                  = group_loci
        params['profile_code2def']      = profile_code2def
        params['gnm2weight']            = gnm2weight
        params['feature_labels']        = feature_labels
        params['file2crispr_type'] = file2crispr_type

        r.write_to_xls_generic_singleton_loci(params)


def generate_cluster_reports(cluster_packs, loci, reports_dir, feature_labels, method, thresholds_pack):

    if not feature_labels:
        local_features = True
    else:
        local_features = False

    thr_occ, thr_crisp, cluster_threshold = thresholds_pack

    summary_file = os.path.join(reports_dir,
                                'summary_%s_%d_%.2f_%.2f.xls' % (method, thr_occ, thr_crisp, cluster_threshold))

    workbook = x.Workbook(summary_file)
    worksheet = workbook.add_worksheet()

    header_format = workbook.add_format()
    header_format.set_font_size(12)
    header_format.set_bold()
    header_format.set_align('center')
    worksheet.set_column(4,5,50)
    worksheet.write_row(0, 0, ["File name", "Weight", "Loci", "Entropy", "systems weight", "systems count"], header_format)

    print "Generating report files"
    ind = 0

    weights = np.zeros(len(cluster_packs))
    entropies = np.zeros(len(cluster_packs))

    for outer_i in range(len(cluster_packs)):

        (cluster, type2count, type2weight, entropy) = cluster_packs[outer_i]

        ind += 1
        cl_files = [os.path.basename(loci[i].file_name) for i in cluster]

        weight = sum([gnm2weight[file2org[file]] for file in cl_files])

        weights[outer_i] = weight
        entropies[outer_i] = entropy

        crispr_cas_types_count = " ; ".join([k+":"+str(v) for (k,v) in sorted(type2count.items(), key=itemgetter(1), reverse=True)])
        crispr_cas_types_weight = " ; ".join([k+":"+str(v) for (k,v) in sorted(type2weight.items(), key=itemgetter(1), reverse=True)])

        xls_file_name = os.path.join(reports_dir, '%d.xls' % ind)

        worksheet.write_row(ind+1, 0, ['%d.xls'%ind,
                                       weight,
                                       len(cl_files),
                                       entropy,
                                       crispr_cas_types_weight,
                                       crispr_cas_types_count,
                                       " "])

        cl_loci = sorted([loci[_i] for _i in cluster], key = lambda x: gnm2weight[x.organism], reverse=True)

        local_profile2weight = {}
        for locus in cl_loci:
            for gene in locus.genes:
                for profile in gene.cogid.split(','):
                    t.update_dictionary(local_profile2weight, profile, gnm2weight[locus.organism])

        global_profile2weight = t.map_global_cdd_profile_count()

        if local_features:
            feature_labels = [ k for k,v in local_profile2weight.items() if v/weight >= 0.5 ]

        params = {}

        params['xls_file_name']         = xls_file_name
        params['loci']                  = cl_loci
        params['weight']                = weight
        params['profile_code2def']      = profile_code2def
        params['gnm2weight']            = gnm2weight
        params['feature_labels']        = feature_labels
        params['file2crispr_type']      = file2crispr_type
        params['local_profile2weight']  = local_profile2weight
        params['global_profile2weight'] = global_profile2weight

        r.write_to_xls_generic_loci(params)

    worksheet.write_row(ind+3, 0, ['Average entropy'], header_format)
    worksheet.write_row(ind+3, 1, [np.sum(weights*entropies)/np.sum(weights)])

    worksheet.write_row(ind + 4, 0, ['Exp(Average entropy)'], header_format)
    worksheet.write_row(ind + 4, 1, [np.exp(np.sum(weights * entropies) / np.sum(weights))])


def plot_dot_score_dendrograms(loci, thresholds_list):

    for thr_occ, thr_crsp in thresholds_list:

        print "Starting for:", thr_occ, thr_crsp
        feature_definition_file = os.path.join(gv.project_data_path, 'cas1402/profiles_%d_%.2f.tab'%(thr_occ, thr_crsp))

        M, feature_labels, feature_weights = generate_data_feature(feature_definition_file, loci=loci)
        M += np.transpose(M)
        M = -1 * np.log(M/np.max(M))
        M[np.diag_indices_from(M)] = 0
        M[np.where(M==np.inf)] = 50

        M_array = ssd.squareform(M)
        Z = linkage(M_array, method='average')

        dendogram_file = os.path.join(gv.project_data_path, 'cas1402/reports/dnd_dot_%d_%.2f.pdf'%(thr_occ, thr_crsp))
        print dendogram_file
        plot_dendrogram(Z, dendogram_file)


if __name__=='__main__':

    # thresholds = [(0, 0.01), (0, 0.05), (0, 0.2) ,(0, 1), (1, 0), (1, 0.6), (1, 0.4), (5, 0.5), (10, 0.6), (10, 0.4)]
    thresholds = [(0, 0.01), (1, 0.4),  (1, 0.6), (5, 0.5), (10, 0.6), (10, 0.4)]
    # thresholds = [(5, 0.5), (10, 0.6), (10, 0.4)]

    # for thr_occ, thr_crisp in thresholds:
    #     filter_feature_profiles(thr_occ, thr_crisp)

    files_path  = os.path.join(gv.project_data_path, 'cas1402/files/')
    pickle_path = os.path.join(gv.project_data_path, 'cas1402/pickle/')
    report_path = os.path.join(gv.project_data_path, 'cas1402/reports/')

    # print "Loading CDD"
    #
    # gi2annotation = t.map_gid2cdd()
    #
    # print "Loading loci"
    # loci = [Locus(os.path.join(files_path, f), gi2annotation) for f in os.listdir(files_path)]
    # loci = [locus for locus in loci if len(locus.genes) > 2]
    #
    # t.dump_compressed_pickle('loci.p.bz2', loci)

    # print "Loading loci"
    # loci = t.load_compressed_pickle('loci.p.bz2')

    # process_jackard_distances(loci)
    process_jackard_distances([])
    sys.exit()
    # loci = t.load_compressed_pickle('/Volumes/hudaiber/Projects/NewSystems/code/Generic/loci.p.bz2')

    # method = "union"
    # results_data_file_name = os.path.join(gv.project_data_path,'cas1402/results_%s.txt'%method)
    # results_image_file_name = os.path.join(gv.project_data_path,'cas1402/clustering_%s.png'%method)
    # plot_results(results_data_file_name, results_image_file_name)

    thresholds = [(5, 0.5)]
    process_dot_product(thresholds, 'union', loci)
    # process_dot_product(thresholds, 'min', loci)
    sys.exit()

    # plot_dot_score_dendrograms(loci)
    # sys.exit()

    # report_feature_clusters_jw(loci)
    # thresholds_pack = (5, 0.5, (0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.75, 4))
    # thresholds_pack = (5, 0.5, (0.75, ))

    # method = 'min'
    # print "thresholds_pack: ",thresholds_pack
    # report_feature_clusters_dot(loci, thresholds_pack, method)

    # method = 'union'
    # print "thresholds_pack: ",thresholds_pack
    # report_feature_clusters_dot(loci, thresholds_pack, method)

    # thresholds_pack = (5, 0.5, 0.75)
    # method = 'union'
    # print "thresholds_pack: ", thresholds_pack
    # singleton_reports(loci, thresholds_pack, method)

    thresholds_pack = (5, 0.5, 0.75)
    method = 'union'
    print "thresholds_pack: ", thresholds_pack
    profile_statistics_after_clustering(loci, thresholds_pack, method)
    
    sys.exit()
    # process_jackard_distances(loci)

    # filter_feature_profiles(20,0.40)

    # out_path = '/Volumes/hudaiber/Projects/NewSystems/data/cas1402/17x13/'
    # files_13xls = ["393.pty", "1223.pty", "2938.pty", "2791.pty", "2774.pty", "67.pty", "169.pty", "2066.pty", "1482.pty", "810.pty", "2020.pty", "70.pty", "2000.pty", "2003.pty", "2004.pty", "2934.pty", "808.pty", "2559.pty", "2561.pty"]
    # files_17xls = ["1219.pty", "1507.pty", "740.pty", "2753.pty", "2938.pty", "1122.pty", "2771.pty", "2066.pty", "2173.pty", "486.pty", "2018.pty", "774.pty", "721.pty", "2010.pty", "2007.pty"]
    #
    # feature_definition_file = os.path.join(gv.project_data_path,
    #                                        'cas1402/profiles/profiles_%d_%.2f.tab' % (5, 0.5))
    #
    # all_files = set(files_13xls + files_17xls)
    #
    # case_loci = [l for l in loci if os.path.basename(l.file_name) in all_files]
    #
    # label_file = open(out_path+'labels.txt','w')
    #
    # [label_file.write("%s\n"%os.path.basename(l.file_name)) for l in case_loci]
    # label_file.close()
    #
    # method = 'union'
    # M, feature_labels, feature_weights = generate_data_feature(feature_definition_file, method, loci=case_loci)
    # M += np.transpose(M)
    #
    # score_file = open(out_path+'union.txt','w')
    #
    # for i in range(M.shape[0]):
    #     score_file.write("%s\n"% "\t".join([str(j) for j in M[i,:]]))
    #
    # score_file.close()