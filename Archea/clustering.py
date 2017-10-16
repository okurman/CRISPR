 #!/usr/bin/python
 # -*- coding: utf-8 -*-

__author__ = 'hudaiber'


import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')


import global_variables as gv
import os
import numpy as np
from sklearn.cluster import KMeans
import pickle
 from lib.db.archea import archea_db as adb
 from lib.utils import tools as t
 import lib.utils.classes as cl

from operator import itemgetter

def cluster_neighborhoods(n_clusters, feature_profiles, conserved_profiles, n_jobs=None):
    print "Generating data matrix for clustering"

    feature_size = len(feature_profiles)
    data_matrix = np.zeros((len(conserved_profiles), feature_size))

    for i in range(len(conserved_profiles)):
        for j in range(len(feature_profiles)):
            if feature_profiles[j] in conserved_profiles[i][1:6]:
                data_matrix[i, j] = 1
    data_matrix *= 100
    if not n_jobs:
        n_jobs = n_clusters/10
    estimator = KMeans(n_clusters=n_clusters, n_jobs=n_jobs)
    print "Starting clustering job(s):"
    print "Input data:", data_matrix.shape
    print "No of clusters: %d, no of jobs: %d" % (n_clusters, n_jobs)
    predictions = estimator.fit_predict(data_matrix)

    estimator_file = os.path.join(gv.project_data_path, 'Archea/clustering/models_predictions/', 'clustering_estimator_%d.p' % n_clusters)
    predictions_file = os.path.join(gv.project_data_path, 'Archea/clustering/models_predictions/', 'clustering_predictions_%d.p' % n_clusters)
    # data_matrix_file = 'files/data_matrix.p'

    print "Clustering finished. Writing the results to files:"
    print "     ", estimator_file
    print "     ", predictions_file
    pickle.dump(estimator, open(estimator_file, 'w'))
    pickle.dump(predictions, open(predictions_file, 'w'))
    # pickle.dump(data_matrix, open(data_matrix_file, 'w'))


def clustering_postprocess(n_clusters, conserved_profiles):

    predictions_file = os.path.join(gv.project_data_path, 'Archea/clustering/models_predictions/clustering_predictions_%d.p' % n_clusters)
    predictions = pickle.load(open(predictions_file))

    all_cluster_profiles = []
    all_cluster_neighborhoods = []

    for i in range(n_clusters):
        indcs = np.where(predictions == i)
        cluster_mapped = [conserved_profiles[i] for i in indcs[0]]
        cluster_profiles = set()
        cluster_neighborhood_ids = []

        for nbr in cluster_mapped:
            cluster_profiles = cluster_profiles | set(nbr[1:6])
            cluster_neighborhood_ids.append(nbr[0])
        cluster_profiles = sorted(list(cluster_profiles))
        all_cluster_profiles.append(cluster_profiles)
        all_cluster_neighborhoods.append(cluster_neighborhood_ids)

    cluster_path = os.path.join(gv.project_data_path, 'Archea/clustering', str(n_clusters))
    if not os.path.exists(cluster_path):
        os.mkdir(cluster_path)

    cluster_profiles_file = os.path.join(cluster_path, "clustered_profiles.txt")
    cluster_neighborhood_file = os.path.join(cluster_path, "clustered_neighborhoods.txt")

    with open(cluster_profiles_file, "w") as f:
        f.write("Cluster no\tCommunity\n")
        cnt = 1
        for cl in all_cluster_profiles:
            f.write(str(cnt) + "\t" + " ".join(cl)+"\n")
            cnt += 1

    with open(cluster_neighborhood_file, "w") as f:
        f.write("Cluster no\tNeighborhood IDs\n")
        cnt = 1
        for cl in all_cluster_neighborhoods:
            f.write(str(cnt) + "\t" + " ".join([str(i) for i in cl])+"\n")
            cnt += 1

    return all_cluster_profiles


def get_popular_profiles(kplet_rows, limit_to):
    profile_stats = {}
    for r in kplet_rows:
        wgt = r[7]
        cnt = r[6]
        for profile in r[1:6]:
            if profile in profile_stats:
                profile_stats[profile].count += cnt
                profile_stats[profile].weight += wgt
            else:
                profile_stats[profile] = cl.ProfileCount(cnt, wgt)

    profile_weights = [(k, v.weight) for k, v in profile_stats.items()]
    profile_weights = sorted(profile_weights, key=itemgetter(1), reverse=True)
    for p in profile_weights[:limit_to]:
        print p


if __name__ == '__main__':

    # profile_weights = open('files/profile_weights.tab').readlines()
    # top_500_profiles = np.asarray([l.strip().split()[1] for l in profile_weights[:500]])
    #
    # conserved_profiles_file = os.path.join(gv.project_data_path, 'Archea/genes_and_flanks/win_10/kplets/pentaplets.csv')
    #
    # conserved_profiles = [(l.split(',')[0], l.split(',')[1:6]) for l in open(conserved_profiles_file).readlines()[1:]]
    #
    # n_clusters = int(sys.argv[1])
    # # n_clusters = 50
    # cluster_neighborhoods(n_clusters, top_500_profiles, conserved_profiles)
    #
    # cluster_profiles = clustering_postprocess(n_clusters, conserved_profiles)

    neighborhoods_path = os.path.join(gv.project_data_path, 'Archea', 'genes_and_flanks', 'win_10', 'pty')
    print 'starting weighted profile'
    neighborhood_profiles = t.get_weighted_profiles_from_neighborhoods(neighborhoods_path, exclude_target=False)
    limit_to = 1000
    feature_profiles = [k[0] for k in neighborhood_profiles[:limit_to]]
    print 'starting kplets retrieval.'
    kplets = adb.get_archaea_kplets()

    n_clusters = 100
    n_jobs = 5
    cluster_neighborhoods(n_clusters, feature_profiles, kplets[:500000])
    clustering_postprocess(n_clusters, kplets)
