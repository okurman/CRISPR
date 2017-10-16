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
import lib_old.classes as cl


def cluster_neighborhoods(n_clusters, feature_profiles, conserved_profiles):
    print "Generating data matrix for clustering"
    data_matrix = np.zeros((len(conserved_profiles), 500))

    for i in range(len(conserved_profiles)):
        for j in range(len(feature_profiles)):
            if feature_profiles[j] in conserved_profiles[i]:
                data_matrix[i, j] = 1
    data_matrix *= 100
    n_jobs = n_clusters/10
    estimator = KMeans(n_clusters=n_clusters, n_jobs=n_jobs)
    print "Starting clustering job(s):"
    print "Input data:", data_matrix.shape
    print "No of clusters: %d, no of jobs: %d" % (n_clusters, n_jobs)
    predictions = estimator.fit_predict(data_matrix)

    estimator_file = os.path.join(gv.project_data_path, 'clustering/models_predictions/', 'clustering_estimator_%d.p' % n_clusters)
    predictions_file = os.path.join(gv.project_data_path, 'clustering/models_predictions/', 'clustering_predictions_%d.p' % n_clusters)
    print "Clustering finished. Writing the results to files:"
    print "     ", estimator_file
    print "     ", predictions_file
    pickle.dump(estimator, open(estimator_file, 'w'))
    pickle.dump(predictions, open(predictions_file, 'w'))


def clustering_postprocess(n_clusters, conserved_profile_nbrhds, conserved_profiles_map, weight_cutoff):

    predictions_file = os.path.join(gv.project_data_path, 'clustering/models_predictions/clustering_predictions_%d.p' % n_clusters)
    predictions = pickle.load(open(predictions_file))

    all_cluster_profiles = []

    for i in range(n_clusters):
        indcs = np.where(predictions == i)
        cluster_mapped = conserved_profile_nbrhds[indcs]
        cluster_profiles = set()
        cluster_mapped = [k for k in cluster_mapped if conserved_profiles_map[k] > weight_cutoff]
        if len(cluster_mapped) == 0:
            continue

        for nbr in cluster_mapped:
            cluster_profiles = cluster_profiles | set(nbr.split('-'))
        cluster_profiles = sorted(list(cluster_profiles))
        all_cluster_profiles.append(cluster_profiles)

    cluster_path = os.path.join(gv.project_data_path, 'clustering', str(n_clusters))
    if not os.path.exists(cluster_path):
        os.mkdir(cluster_path)

    cluster_profiles_file = os.path.join(cluster_path, "cluster_profiles_cutoff_%f.txt" % weight_cutoff)
    with open(cluster_profiles_file, "w") as f:
        for cl in all_cluster_profiles:
            f.write("\t".join(cl)+"\n")

    return all_cluster_profiles


if __name__=='__main__':

    profile_weights = pickle.load(open('files/profile_weights.p'))
    top_500_profiles = np.asarray([k for k, v in profile_weights[:500]])

    conserved_profiles_file = os.path.join(gv.project_data_path, 'genes_and_flanks/win_10/kmers/500_7.tab')
    conserved_profiles = np.asanyarray([l.split()[0] for l in open(conserved_profiles_file).readlines()[1:]])
    conserved_profiles_map = {l.split()[0]: float(l.split()[1]) for l in open(conserved_profiles_file).readlines()[1:]}

    # n_clusters = int(sys.argv[1])
    n_clusters = 10
    # cluster_neighborhoods(n_clusters, top_500_profiles, conserved_profiles)

    weight_cutoff = 1.2
    cluster_profiles = clustering_postprocess(n_clusters, conserved_profiles, conserved_profiles_map, weight_cutoff)