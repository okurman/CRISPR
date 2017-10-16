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

from lib.db.archea import pentaplets as p
from lib.db.archea import quadruplets as q
from lib.db.archea import triplets as tr
from lib.db.archea import duplets as d
from lib.db.archea.db_tools import file2src_src2org_map
from lib.db.archea import neighborhoods_path

# import report_generation as r
from lib.utils import reporting as r
from lib.utils import distributions as dist
import lib.utils.merging.main as merging
import lib.utils.tools as t
import pickle
import bz2
from operator import itemgetter


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

    data_path = os.path.join(gv.project_data_path, 'Archea/pickle/')

    fname = os.path.join(data_path, str(limit_to), 'pentaplets_merged_within.p.bz2')
    pentaplets = t.load_compressed_pickle(fname)
    fname = os.path.join(data_path, str(limit_to), 'quadruplets_merged_within.p.bz2')
    quadruplets = t.load_compressed_pickle(fname)
    fname = os.path.join(data_path, str(limit_to), 'triplets_merged_within.p.bz2')
    triplets = t.load_compressed_pickle(fname)
    fname = os.path.join(data_path, str(limit_to), 'duplets_merged_within.p.bz2')
    duplets = t.load_compressed_pickle(fname)

    # pentaplets, quadruplets, triplets = None, None, None

    print 'Generationg reports for within orders merged lists'

    flank_report_fnames = ['pentaplets_flank_counts.xls', 'quadruplets_flank_counts.xls', \
                           'triplets_flank_counts.xls', 'duplets_flank_counts.xls']

    kplets = [pentaplets, quadruplets, triplets, duplets]
    titles = ['Pentaplets', 'Quadruplets', 'Triplets', 'Duplets']

    report_files_dir = os.path.join(gv.project_data_path, 'Archea/reports/merged_within_orders/', report_dir)

    for i in range(len(flank_report_fnames)):

        _fname = flank_report_fnames[i]
        _kplet_list = kplets[i]
        _title = titles[i]
        flank_counts, cog2gids, gid2weight = dist.get_flank_distributions(_kplet_list, neighborhood_files_path, target_profiles)

        # universal_flank_counts = t.merge_dict_list(flank_counts)
        # universal_cog2gids     = t.merge_dict_set_list(cog2gids, gid2weight)

        universal_flank_counts = t.merge_dict_set_list(cog2gids, gid2weight)

        params = dict()
        params[  'xls_file_name'] = os.path.join(report_files_dir, _fname)
        params[    'profile2def'] = profile2def
        params[   'flank_counts'] = universal_flank_counts
        params[          'title'] = _title
        params['target_profiles'] = target_profiles
        r.write_flanking_count_xls(params)

    if not os.path.exists(report_files_dir):
        os.mkdir(report_files_dir)

    for i, kplet_pool in zip([2, 3, 4, 5], [duplets, triplets, quadruplets, pentaplets]):

        j = 0
        # profile2counts_pool = dist.get_flank_distributions(kplet_pool, neighborhood_files_path, target_profiles)

        file_summaries_pool = []

        for kplet_list in kplet_pool:

            summary_terms = merging.merge_into_file_summaries(kplet_list,
                                                      neighborhood_files_path,
                                                      file2src_src2org_map,
                                                      'archaea')
            file_summaries_pool.append(summary_terms)

        for summary_terms in sorted(file_summaries_pool, key=itemgetter(5), reverse=True):

            src2org, file_summaries, community, community_count, community_count_with_flanks, weight = summary_terms

            if not src2org:
                continue

            cur_reports_folder = os.path.join(report_files_dir, str(i))
            if not os.path.exists(cur_reports_folder):
                os.mkdir(cur_reports_folder)

            xls_file_name = os.path.join(cur_reports_folder,  "%d_%d.xls" % (j+1, i))
            class2counts, class2profiles = merging.arcog_profile_count_into_class_count(community_count)
            class2counts_flank, class2profiles_flank = merging.arcog_profile_count_into_class_count(community_count_with_flanks)

            j += 1
            params = dict()
            params[     'xls_file_name'] = xls_file_name
            params[           'src2org'] = src2org
            params[    'file_summaries'] = file_summaries
            params[         'community'] = community
            params[   'target_profiles'] = target_profiles
            params[       'profile2def'] = profile2def
            params[     'gid2arcog_cdd'] = gid2arcog_cdd
            params[      'class2counts'] = class2counts
            params[    'class2profiles'] = class2profiles
            params['class2counts_flank'] = class2counts_flank
            # params['profile2counts'] = profile2counts
            r.write_to_xls(params)


def get_profiles_counts(data_path):
    # file_names = ['pentaplets', 'quadruplets', 'triplets', 'duplets']
    file_names = ['duplets', 'triplets', 'quadruplets', 'pentaplets']
    print 'Reading merged kplet files'

    for file_name in file_names:
        print 'Loading the file:', file_name
        dump_file = bz2.BZ2File(os.path.join(data_path, '%s_merged_across.p.bz2'%file_name))
        kplets_pool = pickle.load(dump_file)

        print 'Counting community'
        community_count_pool = []
        community_count_pool_with_flanks = []
        for kplets in kplets_pool:
            _src2org, _, _, community_count, community_count_with_flanks = merging.merge_into_file_summaries(kplets,
                                                                         neighborhood_files_path,
                                                                         file2src_src2org_map)
            if not _src2org:
                continue
            community_count_pool.append(community_count)
            community_count_pool_with_flanks.append(community_count_with_flanks)

        dump_file_name = os.path.join(data_path, '%s_community_count.p.bz2'%file_name)
        print 'Dumping into', dump_file_name
        t.dump_compressed_pickle(dump_file_name, community_count_pool)

        dump_file_name = os.path.join(data_path, '%s_community_count_with_flanks.p.bz2'%file_name)
        print 'Dumping into', dump_file_name
        t.dump_compressed_pickle(dump_file_name, community_count_pool_with_flanks)

        print
        print


def generate_pickles(save_path, limit_to):

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    # print "Loading kplets from DB"
    # pentaplets  =  p.get_report_kplets(limit_to=limit_to, load_locations=True)
    # quadruplets =  q.get_report_kplets(limit_to=limit_to, load_locations=True)
    # triplets    = tr.get_report_kplets(limit_to=limit_to, load_locations=True)
    # duplets     =  d.get_report_kplets(limit_to=limit_to, load_locations=True)

    # print "Dumping raw kplet data to files"
    # dump_file = bz2.BZ2File(os.path.join(save_path, 'duplets_raw.p.bz2'), 'w')
    # pickle.dump(duplets, dump_file)
    # dump_file = bz2.BZ2File(os.path.join(save_path, 'triplets_raw.p.bz2'), 'w')
    # pickle.dump(triplets, dump_file)
    # dump_file = bz2.BZ2File(os.path.join(save_path, 'quadruplets_raw.p.bz2'), 'w')
    # pickle.dump(quadruplets, dump_file)
    # dump_file = bz2.BZ2File(os.path.join(save_path, 'pentaplets_raw.p.bz2'), 'w')
    # pickle.dump(pentaplets, dump_file)

    print "Loading raw kplets from pickles"
    dump_file = os.path.join(save_path, 'duplets_raw.p.bz2')
    duplets = t.load_compressed_pickle(dump_file)
    dump_file = os.path.join(save_path, 'triplets_raw.p.bz2')
    triplets= t.load_compressed_pickle(dump_file)
    dump_file = os.path.join(save_path, 'quadruplets_raw.p.bz2')
    quadruplets = t.load_compressed_pickle(dump_file)
    # dump_file = os.path.join(save_path, 'pentaplets_raw.p.bz2')
    # pentaplets = t.load_compressed_pickle(dump_file)

    print "Basic within merging"
    # pentaplets = merging.basic_merge_within_orders(pentaplets)
    quadruplets= merging.basic_merge_within_orders(quadruplets)
    triplets = merging.basic_merge_within_orders(triplets)
    duplets = merging.basic_merge_within_orders(duplets)

    print "Dumping basic merges"
    # dump_file = os.path.join(save_path, 'pentaplets_basic_merged.p.bz2')
    # t.dump_compressed_pickle(dump_file, pentaplets)
    dump_file = os.path.join(save_path, 'quadruplets_basic_merged.p.bz2')
    t.dump_compressed_pickle(dump_file, quadruplets)
    dump_file = os.path.join(save_path, 'triplets_basic_merged.p.bz2')
    t.dump_compressed_pickle(dump_file, triplets)
    dump_file = os.path.join(save_path, 'duplets_basic_merged.p.bz2')
    t.dump_compressed_pickle(dump_file, duplets)

    # print "Loading basic merges"
    # dump_file = os.path.join(save_path, 'pentaplets_basic_merged.p.bz2')
    # pentaplets = t.load_compressed_pickle(dump_file)
    # dump_file = os.path.join(save_path, 'quadruplets_basic_merged.p.bz2')
    # quadruplets = t.load_compressed_pickle(dump_file)
    # dump_file = os.path.join(save_path, 'triplets_basic_merged.p.bz2')
    # triplets = t.load_compressed_pickle(dump_file)
    # dump_file = os.path.join(save_path, 'duplets_basic_merged.p.bz2')
    # duplets = t.load_compressed_pickle(dump_file)
    sys.exit()

    print 'Starting iterative within mergings'
    pentaplets  = merging.merge_kplets_within_orders_iterative_2(pentaplets)
    quadruplets = merging.merge_kplets_within_orders_iterative_2(quadruplets)
    triplets    = merging.merge_kplets_within_orders_iterative_2(triplets)
    duplets     = merging.merge_kplets_within_orders_iterative_2(duplets)

    print "Dumping iterative merges"
    dump_file = os.path.join(save_path, 'pentaplets_iterative_merged.p.bz2')
    t.dump_compressed_pickle(dump_file, pentaplets)
    dump_file = os.path.join(save_path, 'quadruplets_iterative_merged.p.bz2')
    t.dump_compressed_pickle(dump_file, quadruplets)
    dump_file = os.path.join(save_path, 'triplets_iterative_merged.p.bz2')
    t.dump_compressed_pickle(dump_file, triplets)
    dump_file = os.path.join(save_path, 'duplets_iterative_merged.p.bz2')
    t.dump_compressed_pickle(dump_file, duplets)
    
    sys.exit()
    print 'Dumping merged kplet lists to files'
    dump_file = bz2.BZ2File(os.path.join(save_path, 'pentaplets_merged_within.p.bz2'), 'w')
    pickle.dump(pentaplets, dump_file)
    dump_file = bz2.BZ2File(os.path.join(save_path, 'quadruplets_merged_within.p.bz2'), 'w')
    pickle.dump(quadruplets, dump_file)
    dump_file = bz2.BZ2File(os.path.join(save_path, 'triplets_merged_within.p.bz2'), 'w')
    pickle.dump(triplets, dump_file)
    dump_file = bz2.BZ2File(os.path.join(save_path, 'duplets_merged_within.p.bz2'), 'w')
    pickle.dump(duplets, dump_file)

    print 'Starting accross mergings'
    triplets, duplets = merging.merge_kplets_across_orders(triplets, duplets)
    quadruplets, triplets = merging.merge_kplets_across_orders(quadruplets, triplets)
    pentaplets, quadruplets = merging.merge_kplets_across_orders(pentaplets, quadruplets)

    print 'Dumping to files'
    dump_file = bz2.BZ2File(os.path.join(save_path, 'pentaplets_merged_across.p.bz2'), 'w')
    pickle.dump(pentaplets, dump_file)
    dump_file = bz2.BZ2File(os.path.join(save_path, 'quadruplets_merged_across.p.bz2'), 'w')
    pickle.dump(quadruplets, dump_file)
    dump_file = bz2.BZ2File(os.path.join(save_path, 'triplets_merged_across.p.bz2'), 'w')
    pickle.dump(triplets, dump_file)
    dump_file = bz2.BZ2File(os.path.join(save_path, 'duplets_merged_across.p.bz2'), 'w')
    pickle.dump(duplets, dump_file)

    print 'Done for limit_to:', limit_to
    print
    print


if __name__ == '__main__':

    print 'Pre-Loading dictionaries'
    # target_profiles = t.target_profiles()
    # profile2def = t.map_profile2def()
    # gid2arcog_cdd = t.map_gid2arcog_cdd()
    # neighborhood_files_path = neighborhoods_path()
    print "\n"

    data_path = os.path.join(gv.project_data_path, 'Archea/pickle/')
    for limit_to in [1000000]:
        print "Limit_to:", limit_to
        print
        cur_path = os.path.join(data_path, str(limit_to))
        generate_pickles(cur_path, limit_to)
        print 'Done'
        print "------------------------"

    # for limit_to, report_dir in zip([1000000], ['all']):
    #     print "Limit_to:", limit_to
    #     print
    #     generate_plots_from_pickle(limit_to, report_dir, target_profiles, profile2def, gid2arcog_cdd, neighborhood_files_path)
    #     print '\t\t Done'
    #     print "------------------------\n\n"
