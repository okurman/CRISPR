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

import xlsxwriter as x
from lib.db.generic import pentaplets as p
from lib.db.generic import quadruplets as q
from lib.db.generic import triplets as tr
from lib.db.generic import duplets as d
from lib.utils.merging.generic import main as merging
import lib.utils.tools as t
import lib.utils.reporting as r
import dm_tools as dt
from lib.utils.arguments import GenericReportingInput
import time
from operator import itemgetter


from lib.db.generic import map_profiles_id2code_code2def
(_, profile_code2def) = map_profiles_id2code_code2def('cas')

def generate_pickles(save_path, limit_to):

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    print "Loading from DB"
    print "pentaplets"
    pentaplets  =  p.get_report_kplets(profile_id2code, limit_to=limit_to)
    print "quadruplets"
    quadruplets =  q.get_report_kplets(profile_id2code, limit_to=limit_to)
    print "triplets"
    triplets    = tr.get_report_kplets(profile_id2code, limit_to=limit_to)
    print "duplets"
    duplets     =  d.get_report_kplets(profile_id2code, limit_to=limit_to)

    print "Dumping to files"

    dump_file = os.path.join(save_path, 'duplets.p.bz2')
    print dump_file
    t.dump_compressed_pickle(dump_file, duplets)

    dump_file = os.path.join(save_path, 'triplets.p.bz2')
    print dump_file
    t.dump_compressed_pickle(dump_file, triplets)

    dump_file = os.path.join(save_path, 'quadruplets.p.bz2')
    print dump_file
    t.dump_compressed_pickle(dump_file, quadruplets)

    dump_file = os.path.join(save_path, 'pentaplets.p.bz2')
    print dump_file
    t.dump_compressed_pickle(dump_file, pentaplets)


def generate_pickle_order(prefix, order, save_path, limit_to):

    print "Loading from DB"
    if order == 2:
        print 'duplets'
        data_file = 'duplets.p.bz2'
        kplets = d.get_report_kplets(prefix, profile_id2code, limit_to=limit_to)
    elif order == 3:
        print 'triplets'
        data_file = 'triplets.p.bz2'
        kplets = tr.get_report_kplets(prefix, profile_id2code, limit_to=limit_to)
    elif order == 4:
        print 'quadruplets'
        data_file = 'quadruplets.p.bz2'
        kplets = q.get_report_kplets(prefix, profile_id2code, limit_to=limit_to)
    elif order == 5:
        print 'pentaplets'
        data_file = 'pentaplets.p.bz2'
        kplets = p.get_report_kplets(prefix, profile_id2code, limit_to=limit_to)

    # # block for work aorund of too bign pentaplet
    # print 'Loading file2genes'
    # neighborhood_files_path = os.path.join(gv.project_data_path,'CRISPR/datasets/crispr/wgs')
    # _file2genes = {}
    # for _f in os.listdir(neighborhood_files_path):
    #     _file2genes[_f] = dt.get_wgs_file(os.path.join(neighborhood_files_path, _f))
    #
    # kplets = filter_seed(kplets, _file2genes)

    # dump_file = os.path.join(save_path, data_file)
    # print "Dumtiping to file", dump_file
    # t.dump_compressed_pickle(kplets, dump_file)
    # print "Finished"
    # sys.exit()

    dump_file = os.path.join(save_path, data_file)
    print "Dumping to file", dump_file
    t.dump_compressed_pickle(dump_file, kplets)
    print "Finished"


def filter_seed(kplets, file2genes=False):

    # if not file2genes:
    #     _file2genes = {}
    #     for _f in os.listdir(neighborhood_files_path):
    #        _file2genes[_f] = dt.get_wgs_file(os.path.join(neighborhood_files_path, _f))
    #
    new_kplets = list()
    for kplet in kplets:

        _code2hit_count = { code:0 for code in kplet.codes}
        len_files = len(kplet.files)

        for _f in kplet.files:

            _genes = file2genes[_f]
            
            for _gene in _genes:
                for _code in _gene.cogid.split(','):
                    if _code in _code2hit_count and _gene.is_seed == "Seed":
                        _code2hit_count[_code] += 1

        if sum([float(v) for v in _code2hit_count.values()])/len_files > 0.9:
            new_kplets.append(kplet)

    return new_kplets


def process_reporting_single_order_old(order):

    if order == 2:
        data_file = 'iterative_merged_duplets.p.bz2'
    elif order == 3:
        data_file = 'iterative_merged_triplets.p.bz2'
    elif order == 4:
        data_file = 'iterative_merged_quadruplets.p.bz2'
    elif order == 5:
        data_file = 'iterative_merged_pentaplets.p.bz2'

    data_file = os.path.join(data_path, data_file)
    print "Loading file", data_file
    merged_lists = t.load_compressed_pickle(data_file)
    reports_dir = os.path.join(reports_path,'merged',str(order))
    print "Generating reports to:", reports_dir

    generate_reports(merged_lists, reports_dir, neighborhood_files_path)


def process_reporting_single_order(order, threshold):

    if order == 2:
        prefix = 'duplets'
    elif order == 3:
        prefix = 'triplets'
    elif order == 4:
        prefix = 'quadruplets'
    elif order == 5:
        prefix = 'pentaplets'

    pickle_file = os.path.join(data_path, "%s.p.bz2"%prefix)
    print "Loading kplets:", pickle_file

    kplets = t.load_compressed_pickle(pickle_file)
    kplets = {kplet.id:kplet for kplet in kplets}

    merged_data_file = os.path.join(data_path, "%s.merged_%s.txt"%(prefix, threshold.split('.')[1]))
    print "Loading merged file:", merged_data_file
    merged_lists = [[kplets[int(id)] for id in cur_list.split(',')[0].split()] for cur_list in open(merged_data_file).readlines() if cur_list.split(',')[0]]

    reports_dir = os.path.join(gv.project_data_path,'cas1402/reports/thr_%s/' % threshold)
    if not os.path.exists(reports_dir):
        os.mkdir(reports_dir)
    reports_dir = os.path.join(reports_dir, str(order))
    print "Generating reports to:", reports_dir
    generate_reports(merged_lists, reports_dir, neighborhood_files_path)


def generate_reports(merged_lists, reports_dir, neighborhood_files_path):

    if not os.path.exists(reports_dir):
        os.mkdir(reports_dir)

    summary_file = os.path.join(reports_dir, 'summary.xls')
    workbook = x.Workbook(summary_file)
    worksheet = workbook.add_worksheet()

    header_format = workbook.add_format()
    header_format.set_font_size(12)
    header_format.set_bold()
    header_format.set_align('center')
    worksheet.set_column(3,3,50)
    worksheet.write_row(0, 0, ["File name", "Weight", "Loci", "CRISPR/Cas systems"], header_format)

    file_summary_list = []
    filter_weak_hits = False

    for i, kplet_list in enumerate(merged_lists):

        ret = merging.kplet_list_to_file_summaries(kplet_list,
                                                   neighborhood_files_path,
                                                   filter_weak_hits)

        if not ret or not ret.file_summaries:
            continue

        file_summary_list.append(ret)

    file_summaries_list = sorted(file_summary_list, key=lambda x: x.weight, reverse=True)

    ind = 0
    for file_summaries_wrapper in file_summaries_list:
        ind += 1
        t.dump_compressed_pickle('file_summaries_wrapper.p', file_summaries_wrapper)

        xls_file_name = os.path.join(reports_dir, '%d.xls' % ind)
        args = GenericReportingInput()

        args.xls_file_name               = xls_file_name
        args.file_summaries              = file_summaries_wrapper.file_summaries
        args.organisms                   = file_summaries_wrapper.organisms
        args.weight                      = file_summaries_wrapper.weight
        args.profile_code2def            = profile_code2def
        args.local_af_kplet2count        = file_summaries_wrapper.kplet2count_af
        args.local_bf_kplet2count        = file_summaries_wrapper.kplet2count_bf
        args.local_profile2count_bf      = file_summaries_wrapper.profile2count_bf
        args.local_profile2count_af      = file_summaries_wrapper.profile2count_af
        args.cas_type2count              = file_summaries_wrapper.cas_type2count

        r.write_to_xls_generic_kplets(args)

        cas_type_summary = ""

        for (cas_type, count) in sorted(file_summaries_wrapper.cas_type2count.items(), key = itemgetter(1), reverse=True):
            cas_type_summary += "%s : %d ; "%(cas_type, count)

        worksheet.write_row(ind+1, 0, ['%d.xls'%ind, len(file_summaries_wrapper.file_summaries), cas_type_summary])


def merging_pipeline_for_order(order, data_path, load_from_db=False):
    limit_to = 1000000000
    print "starting for ", order
    if load_from_db:
        print "Loading kplets from DB"
        if order == 2:
            kplet_file = 'duplets.p.bz2'
        elif order == 3:
            kplet_file = 'triplets.p.bz2'
        elif order == 4:
            kplet_file = 'quadruplets.p.bz2'
            kplets = q.get_report_kplets(profile_id2code, limit_to=limit_to)
        elif order == 5:
            kplet_file = 'pentaplets.p.bz2'
            kplets = p.get_report_kplets(profile_id2code, limit_to=limit_to)
    else:
        print "Loading kplets from pickle file"
        if order == 2:
            kplet_file = 'duplets.p.bz2'
        elif order == 3:
            kplet_file = 'triplets.p.bz2'
        elif order == 4:
            kplet_file = 'quadruplets.p.bz2'
        elif order == 5:
            kplet_file = 'pentaplets.p.bz2'

        kplet_file_full = os.path.join(data_path, kplet_file)
        print "Loading :",kplet_file_full
        kplets = t.load_compressed_pickle(kplet_file_full)

    print "No of kplets:", len(kplets)

    loci_threshold = 0.7

    print "Starting to merge with loci_threshold:", loci_threshold
    #print "Loading file2genes"
    tic = time.time()
    print "Basic merging"
    merged_lists = merging.basic_merge_within_orders(kplets, loci_threshold)
    print "Basic merging done. Merged lists:", len(merged_lists)
    fname = os.path.join(data_path, "basic_merged_"+"%f_"%loci_threshold+kplet_file)
    print "Dumping basic merging: ", fname
    t.dump_compressed_pickle(fname, merged_lists)

    print "Iterative merging"
    merged_lists = merging.merge_kplets_within_orders_iterative(merged_lists, loci_threshold)
    print "Iterative merging done. Merged lists:", len(merged_lists)
    fname = os.path.join(data_path, "iterative_merged_"+"%f_"%loci_threshold+kplet_file)
    print "Dumping Iterative merging: ",fname
    t.dump_compressed_pickle(fname, merged_lists)
    print "Completed in:", time.time()-tic, "(s)"

    loci_threshold = 0.8

    print "Starting to merge with loci_threshold:", loci_threshold
    #print "Loading file2genes"
    tic = time.time()
    print "Basic merging"
    merged_lists = merging.basic_merge_within_orders(kplets, loci_threshold)
    print "Basic merging done. Merged lists:", len(merged_lists)
    fname = os.path.join(data_path, "basic_merged_"+"%f_"%loci_threshold+kplet_file)
    print "Dumping basic merging: ", fname
    t.dump_compressed_pickle(fname, merged_lists)

    print "Iterative merging"
    merged_lists = merging.merge_kplets_within_orders_iterative(merged_lists, loci_threshold)
    print "Iterative merging done. Merged lists:", len(merged_lists)
    fname = os.path.join(data_path, "iterative_merged_"+"%f_"%loci_threshold+kplet_file)
    print "Dumping Iterative merging: ",fname
    t.dump_compressed_pickle(fname, merged_lists)
    print "Completed in:", time.time()-tic, "(s)"

    loci_threshold = 0.9

    print "Starting to merge with loci_threshold:", loci_threshold
    #print "Loading file2genes"
    tic = time.time()
    print "Basic merging"
    merged_lists = merging.basic_merge_within_orders(kplets, loci_threshold)
    print "Basic merging done. Merged lists:", len(merged_lists)
    fname = os.path.join(data_path, "basic_merged_"+"%f_"%loci_threshold+kplet_file)
    print "Dumping basic merging: ", fname
    t.dump_compressed_pickle(fname, merged_lists)

    print "Iterative merging"
    merged_lists = merging.merge_kplets_within_orders_iterative(merged_lists, loci_threshold)
    print "Iterative merging done. Merged lists:", len(merged_lists)
    fname = os.path.join(data_path, "iterative_merged_"+"%f_"%loci_threshold+kplet_file)
    print "Dumping Iterative merging: ",fname
    t.dump_compressed_pickle(fname, merged_lists)
    print "Completed in:", time.time()-tic, "(s)"

    # csv_kplet_file = os.path.join(data_path, kplet_file.split('.')[0]+".csv")
    # csv_merged_lists_file = os.path.join(data_path, "iterative_merged_"+kplet_file.split('.')[0]+".csv")
    # print "Writing kplets to csv file:"
    # print csv_kplet_file
    #
    # t.write_kplets_to_csv(kplets, csv_kplet_file)
    # tic = time.time()
    # print "Starting kpletmerger with params:"
    # print "kpletmerger", csv_kplet_file, csv_merged_lists_file
    # print "\n\n"
    # sp.call(["kpletmerger",csv_kplet_file,csv_merged_lists_file])
    # print "\n\n"
    # print "Completed in:", time.time()-tic, "(s)"


if __name__ == '__main__':


    threshold = sys.argv[1]
    order = int(sys.argv[2])

    data_path = os.path.join(gv.project_data_path,'cas1402/pickle')
    neighborhood_files_path = os.path.join(gv.project_data_path,'cas1402/files')

    process_reporting_single_order(order, threshold)

    # reports_path = os.path.join(gv.project_data_path,'cas1402/reports')
    # Pickle generation
    # generate_pickle_order(prefix, order, data_path, 1000000000)

    # merging_pipeline_for_order(order, data_path)
    # process_reporting_single_order(order)
