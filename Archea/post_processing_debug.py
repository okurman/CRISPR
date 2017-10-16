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

    data_path = os.path.join(gv.project_data_path, 'Archea/pickle/')

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

    # pickle.dump(pentaplets, dump_file)
    #
    #
    # pentaplets = p.get_report_kplets(limit_to=limit_to, load_locations=True)
    # quadruplets = q.get_report_kplets(limit_to=limit_to, load_locations=True)
    # triplets = tr.get_report_kplets(limit_to=limit_to, load_locations=True)
    # duplets = d.get_report_kplets(limit_to=limit_to, load_locations=True)
    #
    # print 'Starting within mergings'
    pentaplets = merging.merge_kplets_within_orders_iterative(pentaplets)
    # quadruplets = merging.merge_kplets_within_orders_iterative(quadruplets)
    # triplets = merging.merge_kplets_within_orders_iterative(triplets)
    # duplets = merging.merge_kplets_within_orders_iterative(duplets)
    #
    # print 'Starting accross mergings'
    # triplets, duplets = merging.merge_kplets_across_orders(triplets, duplets)
    # quadruplets, triplets = merging.merge_kplets_across_orders(quadruplets, triplets)
    # pentaplets, quadruplets = merging.merge_kplets_across_orders(pentaplets, quadruplets)
    #
    # print 'Dumping to files'
    # dump_file = bz2.BZ2File(os.path.join(save_path, 'pentaplets_merged_across.p.bz2'), 'w')
    # pickle.dump(pentaplets, dump_file)
    # dump_file = bz2.BZ2File(os.path.join(save_path, 'quadruplets_merged_across.p.bz2'), 'w')
    # pickle.dump(quadruplets, dump_file)
    # dump_file = bz2.BZ2File(os.path.join(save_path, 'triplets_merged_across.p.bz2'), 'w')
    # pickle.dump(triplets, dump_file)
    # dump_file = bz2.BZ2File(os.path.join(save_path, 'duplets_merged_across.p.bz2'), 'w')
    # pickle.dump(duplets, dump_file)
    #
    # print 'Done for limit_to:', limit_to
    # print
    # print


if __name__ == '__main__':

    print 'Pre-Loading dictionaries'
    target_profiles = t.target_profiles()
    gid2arcog_cdd = t.map_gid2arcog_cdd()
    profile2def = t.map_profile2def()
    neighborhood_files_path = neighborhoods_path()
    # print "\n"

    # for limit_to, report_dir in zip([500, 1000, 1000000], ['top_500','top_1000','all']):
    #     print "Limit_to:", limit_to
    #     print
    #     generate_plots(limit_to, report_dir, target_profiles, profile2def, gid2arcog_cdd, neighborhood_files_path)
    #     print '\t\tDone'
    #     print "------------------------\n\n"

    # fname = '/Users/hudaiber/Projects/NewSystems/data/Archea/pickle/1000000/triplets_merged_within.p.bz2'
    # kplets = t.load_compressed_pickle(fname)
    #
    # target_profiles = t.target_profiles()
    # flank_counts = dist.get_flank_distributions(kplets, neighborhood_files_path, target_profiles)
    # universal_flank_counts = t.merge_dict_list(flank_counts)
    #
    # xls_file_name = '/Users/hudaiber/Projects/NewSystems/data/Archea/reports/merged_within_orders/all/triplets_flanks_count.xls'
    # params = dict()
    # params['xls_file_name'] = xls_file_name
    # params['profile2def'] = profile2def
    # params['flank_counts'] = universal_flank_counts
    # params['title'] = 'Triplets'
    #
    # r.write_flanking_count_xls(params)
    #
    # data_path = os.path.join(gv.project_data_path,'Archea/pickle/tmp/')
    #
    # # for limit_to in [10000, 1000000]:
    # for limit_to in [1000000]:
    #     print "Limit_to:", limit_to
    #     print
    #     cur_path = os.path.join(data_path, str(limit_to))
    #     generate_pickles(cur_path, limit_to)
    #     print 'Done'
    #     print "------------------------"
    # sys.exit()

    fname = '/Users/hudaiber/Projects/NewSystems/data/Archea/pickle/1000000/pentaplets_merged_within.p.bz2'
    pentaplets = t.load_compressed_pickle(fname)
    #
    # fname = '/Users/hudaiber/Projects/NewSystems/data/Archea/pickle/1000000/5_community_count.p.bz2'
    # community = t.load_compressed_pickle(fname)
    #
    # fname = '/Users/hudaiber/Projects/NewSystems/data/Archea/pickle/1000000/5_community_count_with_flanks.p.bz2'
    # community_flank = t.load_compressed_pickle(fname)
    #
    # fname = '/Users/hudaiber/Projects/NewSystems/data/Archea/pickle/1000000/5_community_count_classes.p.bz2'
    # community_classes = t.load_compressed_pickle(fname)
    #
    # fname = '/Users/hudaiber/Projects/NewSystems/data/Archea/pickle/1000000/5_community_count_with_flanks_classes.p.bz2'
    # community_flank_classes = t.load_compressed_pickle(fname)
    #
    # fname = '/Users/hudaiber/Projects/NewSystems/data/Archea/pickle/1000000/duplets_merged_within.p.bz2'
    # community_flank_classes = t.load_compressed_pickle(fname)
    #
    report_dir = 'all'
    report_files_dir = os.path.join(gv.project_data_path, 'Archea/reports/merged_within_orders/', report_dir)


    kplets = pentaplets

    pivot = kplets[2]
    pivot_orgs = set([file2organism[f] for f in pivot.files])

    for i in range(2, len(kplets)):
        _kplet = kplets[i]
        if len(pivot.codes.intersection(_kplet.codes)) == 4:
            _orgs = set([file2organism[f] for f in _kplet.files])
            print i+1, len(set(pivot.files).intersection(set(_kplet.files))), len(_kplet.files), len(pivot_orgs), len(_orgs), len(pivot_orgs.intersection(_orgs)), _kplet.codes

    profile2counts_pool = dist.get_flank_distributions(pentaplets, neighborhood_files_path, target_profiles)

    j = 0
    for k in range(len(pentaplets)):

        kplet_sublist = pentaplets[k]

        cur_reports_folder = os.path.join(report_files_dir, str(5))

        src2org, file_summaries, community, community_count, community_count_with_flanks = \
            merging.merge_into_file_summaries(kplet_sublist,
                                              neighborhood_files_path,
                                              file2src_src2org_map,
                                              'archaea')
        if not src2org:
            continue

        class2counts, class2profiles = merging.arcog_profile_count_into_class_count(community_count)
        class2counts_flank, class2profiles_flank = merging.arcog_profile_count_into_class_count(community_count_with_flanks)

        profile2counts = profile2counts_pool[k]

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
        params['class2counts'] = class2counts
        params['class2profiles'] = class2profiles
        params['class2counts_flank'] = class2counts_flank
        params['profile2counts'] = profile2counts
        r.write_to_xls(params)
        break