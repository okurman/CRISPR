__author__ = 'hudaiber'

import pickle
import sys

import xlsxwriter as x

from lib.utils import tools as t

if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/lib/BioPy/')
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/lib/BioPy/')
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')
import global_variables as gv
import os

from lib.db.bacteria import db_tools

neighborhood_files_path = os.path.join(gv.project_data_path, 'Bacteria/genes_and_flanks/win_10/raw_nbr_files/')

target_profiles = t.bacteria_target_profiles()
profile2def = t.map_cdd_profile2def()
gid2arcog_cdd = t.map_gid2arcog_cdd()


def write_to_xls(xls_file, kplet):

    kplet_cdd_codes = kplet[0]
    kplet_files = kplet[1]
    # kplet_files = kplet_files.split(',')

    _org2src, _src2files = db_tools.bacteria_org2src_src2files_map(kplet_files)

    workbook = x.Workbook(xls_file)
    worksheet = workbook.add_worksheet()

    row_len = 6
    column_names = ['GI', 'From', 'To', 'Strand', 'CDD', 'Definition']

    title_format = workbook.add_format()
    title_format.set_font_size(14)
    title_format.set_bold()
    title_format.set_align('left')

    header_format = workbook.add_format()
    header_format.set_font_size(12)
    header_format.set_bold()
    header_format.set_align('center')

    target_format = workbook.add_format()
    target_format.set_font_color("red")

    kplet_format = workbook.add_format()
    kplet_format.set_font_color("green")

    top_border = 0
    left_border = 0

    worksheet.merge_range(0, 0, 0, 10, 'k-plet: ' + ' '.join(kplet_cdd_codes), title_format)
    top_border += 1
    worksheet.merge_range(top_border, 0, top_border, 10, 'Organisms: ', title_format)
    top_border += 1

    for org in _org2src.keys():
        worksheet.merge_range(top_border, 2, top_border, 10, org, title_format)
        top_border += 1

    top_border += 1
    for org, srcs in _org2src.items():

        for src in srcs:
            _files = _src2files[src]
            neighborhoods = t.load_neighborhoods(neighborhood_files_path, _files)

            for nbr in neighborhoods:

                cur_top_border = top_border

                if not nbr.flank_extension:
                    nbr.extend_flanks(10, os.path.join(gv.pty_data_path, org, "%s.pty" % src), gid2arcog_cdd)

                worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, "%s %s" % (org, src), header_format)
                cur_top_border += 1
                worksheet.write_row(cur_top_border, left_border, column_names, header_format)

                cur_top_border += 2

                for gene in nbr.genes:
                    cur_cogid = gene.cogid
                    if cur_cogid in target_profiles:
                        data_format = target_format
                    elif cur_cogid in kplet_cdd_codes:
                        data_format = kplet_format
                    else:
                        data_format = workbook.add_format()

                    if gene.tag == 'neighborhood':
                        data_format.set_bg_color('#c4bdbd')

                    if cur_cogid in ["", "-", None]:
                        cur_def = ""
                    else:
                        cur_cogid = cur_cogid.split()
                        if len(cur_cogid) > 0:
                            cur_def = []
                            for k in cur_cogid:
                                if k in profile2def:
                                    cur_def.append(profile2def[k])
                                else:
                                    cur_def.append("")
                            cur_def = " | ".join(cur_def)

                            for c in cur_cogid:
                                if c in target_profiles:
                                    data_format = target_format
                                elif c in kplet_cdd_codes:
                                    data_format = kplet_format

                    data_raw = [gene.gid, gene.pFrom, gene.pTo, gene.strand, gene.cogid, cur_def]
                    worksheet.write_row(cur_top_border, left_border, data_raw, data_format)
                    worksheet.write_row(cur_top_border, left_border+row_len, [" "])
                    cur_top_border += 1
                left_border += row_len + 1
    workbook.close()


def generate_reports_for_experiment(weight_cutoff, n_clusters, neighborhoods=None, clustered_profiles=None, target_profiles=None, profile2def=None):

    if not target_profiles:
        target_profiles = t.target_profiles()
    if not profile2def:
        profile2def = t.map_profile2def()
    if not neighborhoods:
        neighborhoods = pickle.load(open('files/neighborhoods.p'))
    if not clustered_profiles:
        clustered_profiles_file = os.path.join(gv.project_data_path, 'clustering/10/cluster_profiles_cutoff_%1.1f.txt' % weight_cutoff)
        clustered_profiles = [l.strip().split() for l in open(clustered_profiles_file).readlines()]

    reports_file_dir = os.path.join('reports', str(n_clusters))
    if not os.path.exists(reports_file_dir):
        os.mkdir(reports_file_dir)

    reports_file_dir = os.path.join(reports_file_dir, 'weight_cutoff_'+str(weight_cutoff))
    if not os.path.exists(reports_file_dir):
        os.mkdir(reports_file_dir)
    print 'Data ready. Starting report generation.'
    cnt = 1
    for cl in clustered_profiles:
        print 'Cluster no:', cnt
        tmp_neighborhoods = []
        for nbr in neighborhoods:
            tmp_profiles = []
            for g in nbr.genes:
                if g.cogid not in ["", None]:
                    tmp_profiles += g.cogid.split()
            prof_cnt = 0
            for prof in set(tmp_profiles):
                if prof in cl:
                    prof_cnt += 1

            if prof_cnt >= 7:
                tmp_neighborhoods.append(nbr)

        xls_file = os.path.join(reports_file_dir, "cl_no_%d.xls" % cnt)
        write_to_xls(xls_file, tmp_neighborhoods, cl, target_profiles, profile2def)
        cnt += 1


if __name__=='__main__':

    n_clusters = 10
    weight_cutoff=1.2
    generate_reports_for_experiment(weight_cutoff, n_clusters)