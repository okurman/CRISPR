#! /usr/bin/env python
__author__ = 'hudaiber'

import sys
import os

if sys.platform == 'darwin':
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))
elif sys.platform == 'linux2':
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))

import global_variables as gv
sys.path.append(gv.project_code_path)

from lib.db.bacteria import db_tools, merged_neighborhoods_path
from lib.utils import tools as t
from lib.utils.merging import main as merging
from lib.db.bacteria.db_tools import file2src_src2org_map
from lib.utils import reporting as r
import dm_tools as dt

from lib.utils.clustering import Locus
from lib.utils.clustering import scores
from lib.utils.clustering import dendrogram
from lib.utils.clustering import reporting
import numpy as np

# target_profiles = [l.strip() for l in open('/Volumes/pan1/patternquest/Projects/NewSystems/data/Archea/arCOG/selected_arcogs.txt').readlines()]
# neighborhood_files_path = '/Volumes/pan1/patternquest/Projects/NewSystems/data/Archea/genes_and_flanks/win_10/pty/'


if __name__ == '__main__':

    print "Loading maps and profiles"

    target_profiles = t.bacteria_target_profiles()

    additional_profiles = ['pfam13542', 'pfam13005', 'pfam13276', 'pfam13817', 'pfam12323', 'pfam14294',
                           'pfam13700', 'pfam05717', 'pfam13011', 'pfam13551', 'pfam13592', 'pfam13481',
                           'COG3311', 'COG3378', 'pfam14261']

    profile2threshold = {
        "COG3311": 2.0,
        "COG3378": 3.7,
        "pfam05717": 3.0,
        "pfam12323": 4.0,
        "pfam13005": 4.0,
        "pfam13011": 5.0,
        "pfam13276": 4.0,
        "pfam13481": 2.9,
        "pfam13542": 5.0,
        "pfam13551": 5.0,
        "pfam13592": 4.0,
        "pfam13700": 4.0,
        "pfam13817": 4.0,
        "pfam14261": 5.0,
        "pfam14294": 5.0
    }

    profile2def = t.map_cdd_profile2def()
    # gid2arcog_cdd = t.map_gid2arcog_cdd()
    neighborhood_files_path = merged_neighborhoods_path()
    org2weight = t.map_genome2weight()

    pan_data_path = '/panfs/pan1/patternquest/Projects/NewSystems/data/Bacteria/'

    pickle_file = os.path.join(pan_data_path, 'pickle/10000/profile2merged_files.p.bz2')
    profile2files = t.load_compressed_pickle(pickle_file)

    baiticity_file = os.path.join(gv.project_data_path, 'baiticity/bacteria/baiticity.tab')
    profile2baiticity = {l.split()[0]: l.split()[4] for l in open(baiticity_file).readlines()[1:] if l.strip()}

    i = 1

    for highlight_profile in additional_profiles:

        _profile_containing_files = profile2files[highlight_profile]
        file_summaries = merging.get_file_summaries(_profile_containing_files, neighborhood_files_path, org2weight)

        print i, highlight_profile, len(file_summaries)

        report_file = os.path.join(gv.project_data_path,'baiticity/bacteria/reports/%s.xlsx' % highlight_profile)

        if len(file_summaries) == 1:

            clusters = []
            singles = [file_summaries[0].file_name]
            continue

        else:

            cl_loci = [Locus(neighborhood_files_path + '/' + f) for f in _profile_containing_files]
            M = scores.jackard_weighted_scores(cl_loci)
            M += np.transpose(M)
            M = -1 * np.log(M)
            M[np.diag_indices_from(M)] = 0
            M[np.where(M == np.inf)] = 50

            singles, clusters = dendrogram.classify_by_scores(M, profile2threshold[highlight_profile], cl_loci)

        args = {}

        args['xls_file_name']       = report_file
        args['file_summaries']      = file_summaries
        args['clusters']            = clusters
        args['singles']             = singles
        args['target_profiles']     = set(target_profiles)
        args['additional_profiles'] = set(additional_profiles)
        args['profile2def']         = profile2def
        args['profile2baiticity']   = profile2baiticity

        r.write_to_xls_baiticity_profile_set(args)
        i += 1


# dendrogram_file_name = os.path.join(
#     gv.project_data_path + '/baiticity/bacteria/clustering/' + highlight_profile + '_dendrogram.png')
# dendrogram.plot_dendrogram_from_score_matrix(M, dendrogram_file_name)
#
# clustering_file_name = os.path.join(
#     gv.project_data_path + '/baiticity/bacteria/clustering/' + highlight_profile + '_clustering.png')
# clustering_values_file_name = os.path.join(
#     gv.project_data_path + '/baiticity/bacteria/clustering/' + highlight_profile + '_clustering_values.txt')
#
# outf = open(clustering_values_file_name, 'w')
#
# for threshold in [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 4, 5, 6, 7, 8, 9, 10]:
#     try:
#         singles, clusters = dendrogram.classify_by_scores(M, threshold, cl_loci)
#         print "             ", threshold
#     except:
#         break
#     outf.write("%f\t%d\t%d\t%f\t%f\n" % (threshold, len(singles), len(clusters), 0, 0))
#     outf.flush()
# outf.close()
# continue


# i = 1
# profile2files = {}
# print "Loading loci from:", neighborhood_files_path
# for f in os.listdir(neighborhood_files_path):
#
#     genes = dt.get_pty_file(neighborhood_files_path + f)
#     for gene in genes:
#         for profile in gene.cogid.split():
#             t.update_dictionary_set(profile2files, profile, f)
#     i += 1
#     if i%1000 == 0:
#         print i
#
# print "No of loci:", i
# t.dump_compressed_pickle(pickle_file, profile2files)
# print "Done"