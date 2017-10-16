#!/usr/bin/env python
__author__ = 'Sanjarbek Hudaiberdiev'

import os
import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')

import global_variables as gv
sys.path.append(gv.project_code_path)

from lib.db.bacteria import pentaplets as p
from lib.db.bacteria import quadruplets as q
from lib.db.bacteria import triplets as tr
from lib.db.bacteria import duplets as d
from lib.db.bacteria.db_tools import file2src_src2org_map
from lib.db.bacteria import neighborhoods_path

# import report_generation as r
from lib.utils import reporting as r
import lib.utils.merging as merging
import lib.utils.tools as t
import lib.utils.distributions as dist



def generate_plots(limit_to, report_dir, target_profiles, profile2def, gid2arcog_cdd, neighborhood_files_path):

    print 'Reading kplets from database'
    pentaplets = p.get_report_kplets(limit_to=limit_to, load_locations=True)
    quadruplets = q.get_report_kplets(limit_to=limit_to, load_locations=True)
    triplets = tr.get_report_kplets(limit_to=limit_to, load_locations=True)
    duplets = d.get_report_kplets(limit_to=limit_to, load_locations=True)

    print 'Merging within order'
    pentaplets = merging.merge_kplets_within_orders(pentaplets, target_profiles)
    quadruplets= merging.merge_kplets_within_orders(quadruplets, target_profiles)
    triplets = merging.merge_kplets_within_orders(triplets, target_profiles)
    duplets = merging.merge_kplets_within_orders(duplets, target_profiles)

    print 'Generationg reports for within orders merged lists'

    report_files_dir = os.path.join(gv.project_data_path, 'Archea/reports/merged_within_orders/', report_dir)
    if not os.path.exists(report_files_dir):
        os.mkdir(report_files_dir)

    for i, kplet_pool in zip([5, 4, 3, 2], [pentaplets, quadruplets, triplets, duplets]):
        for j, kplet_sublist in enumerate(kplet_pool):
            cur_reports_folder = os.path.join(report_files_dir, str(i))
            if not os.path.exists(cur_reports_folder):
                os.mkdir(cur_reports_folder)
            xls_file_name = os.path.join(cur_reports_folder,  "%d_%d.xls" % (j+1, i))
            r.write_to_xls(xls_file_name,kplet_sublist,target_profiles,profile2def,gid2arcog_cdd,neighborhood_files_path,file2src_src2org_map)

    print 'Merging across order'
    triplets, duplets = merging.merge_kplets_across_order(triplets, duplets)
    quadruplets, triplets = merging.merge_kplets_across_order(quadruplets, triplets)
    pentaplets, quadruplets = merging.merge_kplets_across_order(pentaplets, quadruplets)

    print 'Generationg reports for across orders merged lists'

    report_files_dir = os.path.join(gv.project_data_path, 'Archea/reports/merged_across_orders/', report_dir)
    if not os.path.exists(report_files_dir):
        os.mkdir(report_files_dir)

    for i, kplet_pool in zip([5, 4, 3, 2], [pentaplets, quadruplets, triplets, duplets]):
        for j, kplet_sublist in enumerate(kplet_pool):
            cur_reports_folder = os.path.join(report_files_dir, str(i))
            if not os.path.exists(cur_reports_folder):
                os.mkdir(cur_reports_folder)
            xls_file_name = os.path.join(cur_reports_folder,  "%d_%d.xls" % (j+1, i))
            r.write_to_xls(xls_file_name,kplet_sublist,target_profiles,profile2def,gid2arcog_cdd,neighborhood_files_path,file2src_src2org_map)


def generate_plots_from_pickle(limit_to, report_dir, target_profiles, profile2def, gid2arcog_cdd, neighborhood_files_path):

    data_path = os.path.join(gv.project_data_path, 'Bacteria/pickle/')

    fname = os.path.join(data_path, str(limit_to), 'pentaplets_merged_within.p.bz2')
    pentaplets = t.load_compressed_pickle(fname)
    # fname = os.path.join(data_path, str(limit_to), 'quadruplets_merged_within.p.bz2')
    # quadruplets = t.load_compressed_pickle(fname)
    # fname = os.path.join(data_path, str(limit_to), 'triplets_merged_within.p.bz2')
    # triplets = t.load_compressed_pickle(fname)
    # fname = os.path.join(data_path, str(limit_to), 'duplets_merged_within.p.bz2')
    # duplets = t.load_compressed_pickle(fname)

    print 'Generationg reports for within orders merged lists'

    report_files_dir = os.path.join(gv.project_data_path, 'Archea/reports/merged_within_orders/', report_dir)

    if not os.path.exists(report_files_dir):
        os.mkdir(report_files_dir)

    for i, kplet_pool in zip([5, 4, 3, 2], [pentaplets, pentaplets, pentaplets, pentaplets]):
        print 'Starting for', i
        j = 0
        profile2counts_pool = dist.get_flank_distributions(kplet_pool, neighborhood_files_path, target_profiles)

        for k in range(len(kplet_pool)):

            kplet_sublist = kplet_pool[k]
            cur_reports_folder = os.path.join(report_files_dir, str(i))
            if not os.path.exists(cur_reports_folder):
                os.mkdir(cur_reports_folder)

            src2org, file_summaries, community, community_count, community_count_with_flanks = \
                merging.merge_into_file_summaries(kplet_sublist,
                                                  neighborhood_files_path,
                                                  file2src_src2org_map,
                                                  'archaea')
            if not src2org:
                continue

            xls_file_name = os.path.join(cur_reports_folder,  "%d_%d.xls" % (j+1, i))
            class2counts, class2profiles = merging.arcog_profile_count_into_class_count(community_count)
            class2counts_flank, class2profiles_flank = merging.arcog_profile_count_into_class_count(community_count_with_flanks)

            profile2counts = profile2counts_pool[k]

            j += 1
            params = dict()
            params['xls_file_name']      = xls_file_name
            params['src2org']            = src2org
            params['file_summaries']     = file_summaries
            params['community']          = community
            params['target_profiles']    = target_profiles
            params['profile2def']        = profile2def
            params['gid2arcog_cdd']      = gid2arcog_cdd
            params['class2counts']       = class2counts
            params['class2profiles']     = class2profiles
            params['class2counts_flank'] = class2counts_flank
            params['profile2counts']     = profile2counts
            r.write_to_xls(params)


def generate_pickles(save_path, limit_to):
    if not os.path.exists(save_path):
        os.mkdir(save_path)

    duplets = d.get_report_kplets(limit_to=limit_to, load_locations=True)

    print 'Dumping to files'
    dump_file = os.path.join(save_path, 'duplets_raw.p.bz2')
    t.dump_compressed_pickle(dump_file, duplets)

    print 'Done for limit_to:', limit_to
    print
    print


if __name__ == '__main__':

    print 'Pre-Loading dictionaries'
    target_profiles = t.bacteria_target_profiles()
    gid2arcog_cdd = t.map_gid2arcog_cdd()
    profile2def = t.map_cdd_profile2def()
    neighborhood_files_path = neighborhoods_path()
    print "\n"

    data_path = os.path.join(gv.project_data_path, 'Bacteria/pickle/100000/')
    # files = ['pentaplets_basic_merged.p.bz2', 'quadruplets_basic_merged.p.bz2', 'triplets_basic_merged.p.bz2', 'duplets_basic_merged.p.bz2']

    order = int(sys.argv[1])
    if order == 2:
        files = ['duplets_raw.p.bz2']
    elif order == 3:
        files = ['triplets_raw.p.bz2']
    elif order == 4:
        files = ['quadruplets_raw.p.bz2']
    elif order == 5:
        files = ['pentaplets_raw.p.bz2']

    reports_dir_base = os.path.join(gv.project_data_path,'Bacteria/reports/merged_within_orders/top_100000/raw_ranks')
    if not os.path.exists(reports_dir_base):
        os.mkdir(reports_dir_base)

    for k, kplet_file_name in zip([order], files):
    # for k, kplet_file_name in zip([4, 3, 2], files):
    # for k, kplet_file_name in zip([2, 3], files):

        print "Starting for", k
        kplet_list = t.load_compressed_pickle(os.path.join(data_path, kplet_file_name))
        reports_dir = os.path.join(reports_dir_base, str(k))
	
        if not os.path.exists(reports_dir):
            os.mkdir(reports_dir)

        ind = 0
        for kplet in kplet_list[:1000]:

            codes = list(kplet.codes)
            suffix = ['s' if code in target_profiles else 'n' for code in codes]
            if 's' not in suffix:
                continue

            xls_file_name = os.path.join(reports_dir, '%d_%s.xls' % (ind+1, "".join(suffix)))
            ind += 1

            src2org, file_summaries, community, community_count, community_count_with_flanks, weight = \
                    merging.merge_into_file_summaries([kplet],
                                                      neighborhood_files_path,
                                                      file2src_src2org_map)

            # xls_file_name = os.path.join(reports_dir, '%d.xls' % (ind+1))
            params = dict()
            params[     'xls_file_name'] = xls_file_name
            params[           'src2org'] = src2org
            params[    'file_summaries'] = file_summaries
            params[   'target_profiles'] = target_profiles
            params[       'profile2def'] = profile2def
            params[     'gid2arcog_cdd'] = gid2arcog_cdd

            # print "Writing", ind, xls_file_name
            r.write_to_xls_raw_kplet(params)
