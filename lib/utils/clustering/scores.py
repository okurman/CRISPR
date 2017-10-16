__author__ = 'hudaiber'

import os
import sys
if sys.platform == 'darwin':
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))
elif sys.platform == 'linux2':
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))

# import global_variables as gv

import numpy as np


def jackard_weighted_scores(loci, pickle_file=None, save_to_pickle=None):

    if pickle_file:
        return np.load(pickle_file).items()[0][1]

    no_of_loci = len(loci)
    M = np.zeros((no_of_loci, no_of_loci))

    # print "Starting score calculations"

    for i in range(no_of_loci):
        for j in range(i, no_of_loci):
            M[i, j] = loci[i].score(loci[j])

    # print "Score calculations finished"

    if save_to_pickle and pickle_file:
        np.savez_compressed(pickle_file, data=M)

    return M


def generate_jackard_score_matrix(loci, pickle_file=None):

    # pickle_file = os.path.join(gv.project_data_path, 'cas1402/pickle/jw_scores.npz')

    M = jackard_weighted_scores(loci, pickle_file)

    return M


def generate_dot_product_score_matrix(feature_labels, method, loci=None, save_file=None):

    if not method:
        print "Provide a method for score calculations"
        sys.exit()

    feature_weights = [1]*len(feature_labels)

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

    print "Score calculations finished. Sum of matrix", np.sum(np.sum(M))

    if save_file:
        np.savez_compressed(save_file, data=M)

    if method=='dot':
        M = M / np.max(M)

    return M