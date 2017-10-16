__author__ = 'hudaiber'

import sys
import xlsxwriter as x

if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/lib/BioPy/')
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/lib/BioPy/')
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')

import global_variables as gv
sys.path.append(gv.project_code_path)

from lib.db.archea import db_tools, neighborhoods_path
from lib.utils import tools as t
import os

target_profiles = t.target_profiles()

target_profiles = [l.strip() for l in open('/Volumes/pan1/patternquest/Projects/NewSystems/data/Archea/arCOG/selected_arcogs.txt').readlines()]

profile2def = t.map_profile2def()
gid2arcog_cdd = t.map_gid2arcog_cdd()
neighborhood_files_path = neighborhoods_path()

neighborhood_files_path = '/Volumes/pan1/patternquest/Projects/NewSystems/data/Archea/genes_and_flanks/win_10/pty/'


def write_to_xls(xls_file, kplets):

    community = set()
    [community.update(kplet.codes) for kplet in kplets]
    _file2kplets = {}
    for kplet in kplets:
        for f in kplet.files:
            if f in _file2kplets:
                _file2kplets[f].append(kplet)
            else:
                _file2kplets[f] = [kplet]

    kplet_files= _file2kplets.keys()
    _org2src, _src2files = db_tools.archea_org2src_src2files_map(kplet_files)

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

    target_format_neighborhood = workbook.add_format()
    target_format_neighborhood.set_font_color("red")
    target_format_neighborhood.set_bg_color("#c4bdbd")

    kplet_format = workbook.add_format()
    kplet_format.set_font_color("green")

    kplet_format_neighborhood = workbook.add_format()
    kplet_format_neighborhood.set_font_color("green")
    kplet_format_neighborhood.set_bg_color("#c4bdbd")

    top_border = 0
    left_border = 0

    worksheet.merge_range(0, 0, 0, 10, 'Community: ' + ' '.join(community), title_format)
    top_border += 1

    organisms = sorted(_org2src.keys())
    worksheet.merge_range(top_border, 0, top_border, 10, 'Organisms: %d'%len(organisms), title_format)
    top_border += 1

    worksheet.merge_range(top_border, 0, top_border, 10, ' '.join(organisms))
    top_border += 1

    top_border += 1

    outliers = []

    for org in sorted(_org2src.keys()):
        srcs = _org2src[org]

        for src in srcs:
            _files = _src2files[src]
            neighborhoods = t.load_neighborhoods(neighborhood_files_path, _files)

            for nbr in neighborhoods:

                cur_kplets = _file2kplets[os.path.basename(nbr.source_file)]

                if len(cur_kplets) <= 5:
                    outliers.append([org, src, nbr])
                    continue

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
                        data_format = target_format_neighborhood if gene.tag == 'neighborhood' else target_format
                    elif cur_cogid in community:
                        data_format = kplet_format_neighborhood if gene.tag == 'neighborhood' else kplet_format
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
                                    data_format = target_format_neighborhood if gene.tag == 'neighborhood' else target_format
                                    break
                                if c in community:
                                    data_format = kplet_format_neighborhood if gene.tag == 'neighborhood' else kplet_format
                                    break

                    data_raw = [gene.gid, gene.pFrom, gene.pTo, gene.strand, gene.cogid, cur_def]
                    worksheet.write_row(cur_top_border, left_border, data_raw, data_format)
                    worksheet.write_row(cur_top_border, left_border+row_len, [" "])
                    cur_top_border += 1

                cur_top_border += 2
                worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, "Kplets:")
                cur_top_border += 1
                worksheet.write_row(cur_top_border, left_border, ["Id", "Profiles"])
                cur_top_border += 1

                for kplet in cur_kplets:
                    worksheet.write_row(cur_top_border, left_border, [kplet.id, " ".join(kplet.codes)])
                    cur_top_border += 1
                left_border += row_len + 1

    for outlier in outliers:

        [org, src, nbr] = outlier
        cur_kplets = _file2kplets[os.path.basename(nbr.source_file)]

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
                data_format = target_format_neighborhood if gene.tag == 'neighborhood' else target_format
            elif cur_cogid in community:
                data_format = kplet_format_neighborhood if gene.tag == 'neighborhood' else kplet_format
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
                            data_format = target_format_neighborhood if gene.tag == 'neighborhood' else target_format
                            break
                        if c in community:
                            data_format = kplet_format_neighborhood if gene.tag == 'neighborhood' else kplet_format
                            break

            data_raw = [gene.gid, gene.pFrom, gene.pTo, gene.strand, gene.cogid, cur_def]
            worksheet.write_row(cur_top_border, left_border, data_raw, data_format)
            worksheet.write_row(cur_top_border, left_border+row_len, [" "])
            cur_top_border += 1

        cur_top_border += 2
        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, "Kplets:")
        cur_top_border += 1
        worksheet.write_row(cur_top_border, left_border, ["Id", "Profiles"])
        cur_top_border += 1

        for kplet in cur_kplets:
            worksheet.write_row(cur_top_border, left_border, [kplet.id, " ".join(kplet.codes)])
            cur_top_border += 1
        left_border += row_len + 1

    workbook.close()


if __name__=='__main__':


    # print 'Pentaplets'
    # kplets = db_p.get_report_kplets()
    # reports_file_dir = os.path.join('reports', '5')
    #
    # cnt = 0
    # for kplet in kplets:
    #
    #     xls_file_name = os.path.join(reports_file_dir, "%d.xls" % (cnt+1))
    #     write_to_xls_2(xls_file_name, kplet)
    #     print xls_file_name
    #     cnt += 1
    #
    #
    # print 'quadruplets'
    # from lib.db.archea import quadruplets as db
    # kplets = db.get_report_kplets()
    # reports_file_dir = os.path.join('reports', '4')
    #
    # cnt = 0
    # for kplet in kplets:
    #
    #     xls_file_name = os.path.join(reports_file_dir, "%d.xls" % (cnt+1))
    #     write_to_xls_2(xls_file_name, kplet)
    #     # print xls_file_name
    #     cnt += 1
    #
    # print 'triplets'
    # from lib.db.archea import triplets as db
    # kplets = db.get_report_kplets()
    # reports_file_dir = os.path.join('reports', '3')
    #
    # cnt = 0
    # for kplet in kplets:
    #
    #     xls_file_name = os.path.join(reports_file_dir, "%d.xls" % (cnt+1))
    #     write_to_xls_2(xls_file_name, kplet)
    #     # print xls_file_name
    #     cnt += 1

    important_profiles = set(["arCOG06216", "pfam09820", "pfam08011", "arCOG06699", "pfam14311", "arCOG03352", "arCOG05030", "arCOG03883", "arCOG11984", "arCOG07929", "arCOG10297", "arCOG03884", "arCOG06013", "arCOG08100", "arCOG08494", "arCOG03317", "arCOG08128", "arCOG08578", "arCOG08908", "arCOG04926", "arCOG03423", "arCOG00885", "arCOG10828", "arCOG06279", "arCOG06947", "arCOG05200", "arCOG06154", "arCOG02801", "arCOG02361", "arCOG00108", "arCOG04494", "arCOG10863", "arCOG06914", "arCOG04495"])

    print 'duplets'
    from lib.db.archea import duplets as db
    kplets = db.get_report_kplets(limit_to=100000)

    kplets = [kplet for kplet in kplets if len(kplet.codes.intersection(target_profiles)) > 0]

    kplets = [kplet for kplet in kplets if len(kplet.codes.intersection(important_profiles)) > 0]


    reports_file_dir = os.path.join('reports', '2')

    cnt = 0
    for kplet in kplets:

        xls_file_name = os.path.join(reports_file_dir, "%d.xls" % (cnt+1))
        write_to_xls_2(xls_file_name, kplet)
        # print xls_file_name
        cnt += 1