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

from lib.utils import tools as t

sys.setrecursionlimit(100000)

import lib.utils.clustering.scores as scores
import lib.utils.clustering.reporting as reporting
import lib.utils.clustering.dendrogram as dendrogram


from lib.utils import Locus


def process_jackard_distances(loci):

    score_file = os.path.join(gv.project_data_path, 'cas1402/pickle/jw_scores.npz')
    M = scores.jackard_weighted_scores(loci, score_file)

    M += np.transpose(M)
    M = np.negative(np.log(M))
    np.fill_diagonal(M, 0)
    inf_idx = np.isinf(M)
    M[inf_idx] = 100

    result_file = os.path.join(gv.project_data_path, 'cas1402/results_jw.txt')

    if not os.path.exists(result_file) or os.path.getsize(result_file) == 0:

        clustering_results = []

        for threshold in [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 4, 5, 6, 7, 8, 9, 10]:
            print "Threshold:", threshold
            singles, clusters, sum_errors, entropies = dendrogram.classify_by_scores(M, threshold, loci)

            clustering_results.append([threshold, len(singles), len(clusters), sum_errors, entropies])

        outf = open(result_file, 'w')
        for i in range(len(result_file)-1):
            print i, clustering_results[i]
            outf.write("%f\t%d\t%d\t%f\t%f\n" % (clustering_results[i][0],
                                                 clustering_results[i][1],
                                                 clustering_results[i][2],
                                                 clustering_results[i][3],
                                                 clustering_results[i][4]))

        print "Wrote the results to file:", result_file
        outf.close()

    else:

        clustering_results = [l.strip().split() for l in open(result_file).readlines()]

    clustering_results = np.asarray(clustering_results, dtype=np.float)

    thresholds = clustering_results[:,0]
    singles = clustering_results[:,1]
    clusters = clustering_results[:,2]
    mean_errors = clustering_results[:,3] / 100
    entropies = clustering_results[:,4] * 1000

    plt.figure(figsize=(30,20))

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

    plt.savefig(os.path.join(gv.project_data_path, 'cas1402/clustering_jw.png'))


def get_feature_labels(feature_definition_file):

    feature_labels = []

    for id in open(feature_definition_file).readlines():
        id = id.strip()
        feature_labels.append("CLUSTER_%s_%s" % (id, cluster_id2effective_size[id]))

    return feature_labels


def process_dot_product(thresholds, method, loci=[]):

    for thr_occurrence, thr_baiticity in thresholds:

        results_data_file_name = os.path.join(gv.project_data_path,
                                              'cas4/clustering_experiments_results_%s_%d_%.2f.txt' %
                                              (method, thr_occurrence, thr_baiticity))

        outf = open(results_data_file_name, 'w')

        _fname = "profiles_%d_%.2f.tab" % (thr_occurrence, thr_baiticity)
        feature_definition_file = os.path.join(gv.project_data_path, 'cas4/profiles/', _fname)

        feature_labels = get_feature_labels(feature_definition_file)

        print "Running for: thr_occurrence: %f, thr_baiticity: %f" % (thr_occurrence, thr_baiticity)

        M, _ = scores.generate_dot_product_score_matrix(feature_labels, method=method, loci=loci)

        M += np.transpose(M)
        M = -1 * np.log(M)
        M[np.diag_indices_from(M)] = 0
        M[np.where(M==np.inf)] = 100

        outf.write("#%s\t%d\n"%(_fname, len(feature_labels)))

        print "Clustering with threshold:"
        for threshold in [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3,4,5,6,7,8,9,10]:
            print "         ", threshold
            try:
                singles, clusters, entropies = dendrogram.classify_by_scores_cas4(M, threshold, loci)
            except:
                break

            outf.write("%f\t%d\t%d\t%f\n" % (threshold, len(singles), len(clusters), entropies))
            outf.flush()

    outf.close()

    print "Results saved:", results_data_file_name


def reporting_dot_product_run(loci, thresholds_pack):

    # thresholds_pack = (100, 0.5, (1.25, 1.5))
    method = 'union'
    print "thresholds_pack: ", thresholds_pack

    thr_occurrence, thr_baiticity, _ = thresholds_pack
    _fname = "profiles_%d_%.2f.tab" % (thr_occurrence, thr_baiticity)
    feature_definition_file = os.path.join(gv.project_data_path, 'cas4/profiles/', _fname)

    feature_labels = get_feature_labels(feature_definition_file)

    reporting.report_clustering_dot_product(loci, thresholds_pack, method, feature_labels)


def reporting_jackard_distance(loci):

    thresholds = (4.0, 3.0, 2.0, 1.5, 1.0)

    reporting.report_clustering_jw(loci, thresholds)


if __name__ == '__main__':

    files_path  = os.path.join(gv.project_data_path, 'cas4/files/')
    pickle_path = os.path.join(gv.project_data_path, 'cas4/pickle/')
    report_path = os.path.join(gv.project_data_path, 'cas4/reports/')

        # print "Loading loci"
        #
        # def_file = os.path.join(gv.project_data_path, 'cas4/profiles/defenseProfiles.tab')
        # profile2gene={}
        #
        # for l in open(def_file):
        #     terms = l.split('\t')
        #     profile = terms[0]
        #     gene_names = terms[3].split(',')
        #     if len(gene_names)>1:
        #         profile2gene[profile] = gene_names[1]
        #     else:
        #         profile2gene[profile] = gene_names[0]
        #
        # cdd_profile2gene = t.map_cdd_profile2gene_name()
        # cdd_profile2gene.update(profile2gene)
        #
        # loci = [Locus(os.path.join(files_path, f), file_format='generic', profile2gene=cdd_profile2gene) for f in os.listdir(files_path)]

    cluster_description_file = os.path.join(gv.project_data_path, 'cas4/baiticity/cluster_id_size.tsv')
    cluster_id2effective_size = {l.split()[0]: l.strip().split()[1] for l in open(cluster_description_file)}

    # method = 'union'
    # thresholds = [(0, 0.5), (10, 0.5)]
    # thresholds = [(0, 0.5), (0, 0.5)]
    # process_dot_product(thresholds, method, loci)
    #
    # results_data_file_name = os.path.join(gv.project_data_path, 'cas4/clustering_results.tsv')
    # results_image_file_name = os.path.join(gv.project_data_path, 'cas4/clustering_results.png')

    # reporting.plot_results(results_data_file_name, results_image_file_name)
    # threshold_pack = (10, 0.5, (1, 1.25, 1.5, 1.75, 2))
    threshold_pack = (10, 0.5, (2,))
    reporting_dot_product_run(loci, threshold_pack)