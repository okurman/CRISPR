__author__ = 'hudaiber'


import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
    sys.path.append('/Users/hudaiber/Projects/lib/BioPy/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')
import global_variables as gv
import os
import numpy as np
from sklearn.cluster import KMeans
import pickle
import lib_old.classes as cl
import lib_old.merger as m
import dm_tools


def cluster_neighborhoods(n_clusters, feature_profiles, conserved_profiles):
    print "Generating data matrix for clustering"
    data_matrix = np.zeros((len(conserved_profiles), 500))

    for i in range(len(conserved_profiles)):
        for j in range(len(feature_profiles)):
            if feature_profiles[j] in conserved_profiles[i][0]:
                data_matrix[i, j] = 1
    data_matrix *= 100
    n_jobs = n_clusters/10
    estimator = KMeans(n_clusters=n_clusters, n_jobs=n_jobs)
    print "Starting clustering job(s):"
    print "Input data:", data_matrix.shape
    print "No of clusters: %d, no of jobs: %d" % (n_clusters, n_jobs)
    predictions = estimator.fit_predict(data_matrix)

    estimator_file = os.path.join('/Users/hudaiber/Projects/NewSystems/data/Bacteria/test/clustering_estimator_%d.p' % n_clusters)
    predictions_file = os.path.join('/Users/hudaiber/Projects/NewSystems/data/Bacteria/test/clustering_predictions_%d.p' % n_clusters)
    print "Clustering finished. Writing the results to files:"
    print "     ", estimator_file
    print "     ", predictions_file
    pickle.dump(estimator, open(estimator_file, 'w'))
    pickle.dump(predictions, open(predictions_file, 'w'))


def clustering_postprocess(n_clusters, conserved_profiles, conserved_profiles_map, weight_cutoff=0):

    predictions_file = os.path.join(gv.project_data_path, '/Users/hudaiber/Projects/NewSystems/data/Bacteria/test/clustering/clustering_predictions_%d.p' % n_clusters)
    predictions = pickle.load(open(predictions_file))

    all_cluster_profiles = []
    all_cluster_files = []

    kplets = np.asanyarray([l[0] for l in conserved_profiles])
    kplet_files = np.asanyarray([l[1] for l in conserved_profiles])

    for i in range(n_clusters):
        indcs = np.where(predictions == i)
        cluster_mapped = kplets[indcs]
        cluster_files = kplet_files[indcs]
        cluster_profiles = set()
        indcs = []
        for j in range(len(cluster_mapped)):
            if conserved_profiles_map[cluster_mapped[j]] >= weight_cutoff:
                indcs.append(j)

        cluster_mapped = cluster_mapped[indcs]
        cluster_files = cluster_files[indcs]

        if len(cluster_mapped) == 0:
            all_cluster_files.append([])
            all_cluster_profiles.append([])
            continue

        tmp_cluster_files = set()
        for l in cluster_files:
            tmp_cluster_files.update(l.split())

        all_cluster_files.append(tmp_cluster_files)
        for kplet in cluster_mapped:
            cluster_profiles = cluster_profiles | set(kplet.split('-'))
        cluster_profiles = sorted(list(cluster_profiles))
        all_cluster_profiles.append(cluster_profiles)

    cluster_path = '/Users/hudaiber/Projects/NewSystems/data/Bacteria/test/clustering/'
    if not os.path.exists(cluster_path):
        os.mkdir(cluster_path)

    cluster_profiles_file = os.path.join(cluster_path, "cluster_profiles_cutoff_%f.txt" % weight_cutoff)
    print "Writing summery file"
    with open(cluster_profiles_file, "w") as f:
        for i in range(n_clusters):
            f.write("Cluster: %d\n" % (i+1))
            f.write("Community:\t"+" ".join(all_cluster_profiles[i])+"\n")
            f.write("Files:\t"+" ".join(all_cluster_files[i])+"\n")

    print "Merging cluster files"
    merge_cluster_files(cluster_profiles_file)
    return all_cluster_profiles, all_cluster_files


def merge_cluster_files(summary_file):

    src_to_org = dm_tools.src2org_map()
    dest_path = '/Users/hudaiber/Projects/NewSystems/data/Bacteria/test/clustering/10/'
    with open(summary_file) as infile:
        lines = infile.readlines()
        cnt = 1
        for i in range(0, 29, 3):
            files = lines[i+2].split()[1:]
            src_path = '/Users/hudaiber/Projects/NewSystems/data/Bacteria/test/genes_and_flanks/pty/'
            src_to_gid_map = m.build_src2gids_map([os.path.join(src_path, f) for f in files])

            cur_dest_path = os.path.join(dest_path, str(cnt))
            if not os.path.exists(cur_dest_path):
                os.mkdir(cur_dest_path)

            m.merge_neighborhoods(src_to_gid_map, src_to_org, cur_dest_path, src_path)
            print 'Done for cluster:', cnt, lines[i]
            cnt += 1


if __name__ == '__main__':

    profile_weights = pickle.load(open('/Users/hudaiber/Projects/NewSystems/data/Bacteria/test/profile_weights.p'))
    top_500_profiles = np.asarray([k for k, v in profile_weights[:500]])

    conserved_profiles_file = os.path.join(gv.project_data_path, '/Users/hudaiber/Projects/NewSystems/data/Bacteria/test/500_5.tab')
    conserved_profiles = [[l.split('\t')[0], l.split('\t')[3].strip()] for l in open(conserved_profiles_file).readlines()[1:]]
    conserved_profiles_map = {l.split('\t')[0]: float(l.split('\t')[1]) for l in open(conserved_profiles_file).readlines()[1:]}

    # n_clusters = int(sys.argv[1])
    n_clusters = 10
    # cluster_neighborhoods(n_clusters, top_500_profiles, conserved_profiles)

    # weight_cutoff = 1.2
    # cluster_profiles = clustering_postprocess(n_clusters, conserved_profiles, conserved_profiles_map, weight_cutoff)
    cluster_profiles, cluster_files = clustering_postprocess(n_clusters, conserved_profiles, conserved_profiles_map)



    # merge_cluster_files('/Users/hudaiber/Projects/NewSystems/data/Bacteria/test/clustering/cluster_profiles_cutoff_0.000000.txt')