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

from lib.db.archea import db_tools, merged_neighborhoods_path
from lib.utils import tools as t
from lib.db.archea import duplets as db
from lib.utils.merging import main as merging
from lib.db.archea.db_tools import file2src_src2org_map
from lib.utils import reporting as r
import dm_tools as dt
from lib.utils.clustering import Locus
from lib.utils.clustering import scores
from lib.utils.clustering import dendrogram
from lib.utils.clustering import reporting
import numpy as np
# import matplotlib.pyplot as plt

# target_profiles = [l.strip() for l in open('/Volumes/pan1/patternquest/Projects/NewSystems/data/Archea/arCOG/selected_arcogs.txt').readlines()]
# neighborhood_files_path = '/Volumes/pan1/patternquest/Projects/NewSystems/data/Archea/genes_and_flanks/win_10/pty/'


if __name__ == '__main__':

    print "Loading maps and profiles"

    target_profiles = t.target_profiles()
    additional_profiles = ["arCOG06216", "pfam09820", "pfam08011", "arCOG06699", "pfam14311", "arCOG03352", "arCOG05030",
                           "arCOG03883", "arCOG11984", "arCOG07929", "arCOG10297", "arCOG03884", "arCOG06013", "arCOG08100",
                           "arCOG08494", "arCOG03317", "arCOG08128", "arCOG08578", "arCOG08908", "arCOG04926", "arCOG03423",
                           "arCOG00885", "arCOG10828", "arCOG06279", "arCOG06947", "arCOG05200", "arCOG06154", "arCOG02801",
                           "arCOG02361", "arCOG00108", "arCOG04494", "arCOG10863", "arCOG06914", "arCOG04495"]

    profile2threshold = {"arCOG00108": 2.8,
                         "arCOG00885": 2.0,
                         "arCOG02361": 5,
                         "arCOG02801": 3.0,
                         "arCOG03317": 3.25,
                         "arCOG03352": 2.25,
                         "arCOG03423": 3.0,
                         "arCOG03883": 2.65,
                         "arCOG03884": 2.0,
                         "arCOG04494": 3.5,
                         "arCOG04495": 4.0,
                         "arCOG04926": 3.5,
                         "arCOG05030": 2.25,
                         "arCOG05200": 2.0,
                         "arCOG06013": 4.0,
                         "arCOG06154": 2.75,
                         "arCOG06216": 2.25,
                         "arCOG06279": 3.5,
                         "arCOG06699": 3.0,
                         "arCOG06914": 4.0,
                         "arCOG06947": 3.0,
                         "arCOG07929": 4.0,
                         "arCOG08100": 3.1,
                         "arCOG08128": 2.75,
                         "arCOG08494": 3.0,
                         "arCOG08578": 4.0,
                         "arCOG08908": 2.0,
                         "arCOG10297": 3.4,
                         "arCOG10828": 2.0,
                         "arCOG10863": 4.0,
                         "arCOG11984": 3.5,
                         "pfam08011": 5.0,
                         "pfam09820": 2.0}

    profile2def = t.map_cdd_ar14_profile2def()

    baiticity_file = os.path.join(gv.project_data_path, 'baiticity/archaea/baiticity.tab')
    profile2baiticity = {l.split()[0]: l.split()[4] for l in open(baiticity_file).readlines()[1:]}

    # gid2arcog_cdd = t.map_gid2arcog_cdd()
    neighborhood_files_path = merged_neighborhoods_path()
    org2weight = t.map_archaea_genome2weight()
    pan_data_path = '/panfs/pan1/patternquest/Projects/NewSystems/data/Archea/'

    pickle_file = pickle_file = os.path.join(pan_data_path, 'pickle/10000/profile2merged_files.p.bz2')
    profile2files = t.load_compressed_pickle(pickle_file)

    i = 1

    for highlight_profile in additional_profiles:

        _profile_containing_files = profile2files[highlight_profile]
        file_summaries = merging.get_file_summaries(_profile_containing_files,
                                                    neighborhood_files_path,
                                                    org2weight)

        print i, highlight_profile, len(file_summaries)

        if len(file_summaries) == 1:

            clusters = []
            singles = [file_summaries[0].file_name]

        else:

            cl_loci = [Locus(neighborhood_files_path+'/' + f) for f in _profile_containing_files]
            M = scores.jackard_weighted_scores(cl_loci)
            M += np.transpose(M)
            M = -1 * np.log(M)
            M[np.diag_indices_from(M)] = 0
            M[np.where(M == np.inf)] = 50

            singles, clusters = dendrogram.classify_by_scores(M, profile2threshold[highlight_profile], cl_loci)

        report_file = os.path.join(gv.project_data_path,'baiticity/archaea/reports/%s.xlsx' % highlight_profile)

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

        print report_file
        sys.exit()


# profile2files = {}
# i=0
# for f in os.listdir(neighborhood_files_path):
#
#     genes = dt.get_pty_file(neighborhood_files_path + f)
#     for gene in genes:
#         for profile in gene.cogid.split():
#             t.update_dictionary_set(profile2files, profile, f)
#     i += 1
#
# print "Number of loci,", i
#
# t.dump_compressed_pickle(pickle_file, profile2files)



# dendrogram_file_name = os.path.join(gv.project_data_path+'/baiticity/archaea/clustering/'+highlight_profile+'_dendrogram.png')
# dendrogram.plot_dendrogram_from_score_matrix(M, dendrogram_file_name)
#
# clustering_file_name = os.path.join(gv.project_data_path + '/baiticity/archaea/clustering/' + highlight_profile + '_clustering.png')
# clustering_values_file_name = os.path.join(gv.project_data_path + '/baiticity/archaea/clustering/' + highlight_profile + '_clustering_values.txt')
#
# outf = open(clustering_values_file_name, 'w')
#
# for threshold in [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3,4,5,6,7,8,9,10]:
#     try:
#         singles, clusters = dendrogram.classify_by_scores(M, threshold, cl_loci)
#     except:
#         break
#     outf.write("%f\t%d\t%d\t%f\t%f\n" % (threshold, len(singles), len(clusters), 0, 0))
#     outf.flush()
# outf.close()