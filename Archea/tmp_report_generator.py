#!/sw/bin/python2.7
__author__ = 'Sanjarbek Hudaiberdiev'

import os
import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')

import global_variables as gv
sys.path.append(gv.project_code_path)

from lib.db.archea import pentaplets as p
from lib.db.archea import quadruplets as q
from lib.db.archea import triplets as tr
from lib.db.archea import duplets as d
from lib.db.archea.db_tools import file2src_src2org_map
from lib.db.archea import neighborhoods_path

# import report_generation as r
from lib.utils import reporting as r
import lib.utils.merging as merging
import lib.utils.tools as t
import pickle


if __name__ == '__main__':

    print 'Pre-Loading dictionaries'
    target_profiles = t.target_profiles()
    profile2def = t.map_profile2def()
    gid2arcog_cdd = t.map_gid2arcog_cdd()
    neighborhood_files_path = neighborhoods_path()
    print "\n"

    limit_to = 1000000
    data_path = os.path.join(gv.project_data_path,'Archea/pickle/')
    fname = os.path.join(data_path, str(limit_to), 'pentaplets_merged_across.p.bz2')
    pentaplets = t.load_compressed_pickle(fname)

    report_dir = 'all'
    report_files_dir = os.path.join(gv.project_data_path, 'Archea/reports/merged_across_orders/', report_dir)

    j = 0
    for kplet_sublist in pentaplets:
        cur_reports_folder = os.path.join(report_files_dir, str(5))

        src2org, file_summaries, community, community_count, community_count_with_flanks = \
            merging.merge_into_file_summaries(kplet_sublist,
                                              neighborhood_files_path,
                                              file2src_src2org_map,
                                              'archaea')
        if not src2org:
            continue
        community_classes = merging.arcog_profile_count_into_class_count(community_count)
        community_flank_classes = merging.arcog_profile_count_into_class_count(community_count_with_flanks)

        xls_file_name = os.path.join(cur_reports_folder,  "%d_%d.xls" % (j+1, 5))
        print xls_file_name
        j += 1
        params = dict()
        params['xls_file_name'] = xls_file_name
        params['src2org'] = src2org
        params['file_summaries'] = file_summaries
        params['community'] = community
        params['target_profiles'] = target_profiles
        params['profile2def'] = profile2def
        params['gid2arcog_cdd'] = gid2arcog_cdd
        params['class_counts'] = community_classes
        params['class_flank_counts'] = community_flank_classes
        r.write_to_xls(params)
        break