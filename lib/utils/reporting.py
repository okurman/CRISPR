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
from operator import itemgetter
import os
from lib.utils import tools

# bacteria_target_profiles = set(tools.bacteria_target_profiles())
# profile2def = tools.map_cdd_profile2def()
# gnm2weight = tools.map_genome2weight()


def write_to_xls(params):

    xls_file           = params['xls_file_name']
    src2org            = params['src2org']
    file_summaries     = params['file_summaries']
    community          = params['community']
    target_profiles    = params['target_profiles']
    profile2def        = params['profile2def']
    gid2arcog_cdd      = params['gid2arcog_cdd']
    class2counts       = params['class2counts']
    class2profiles     = params['class2profiles']
    class2counts_flank = params['class2counts_flank']
    # profile2counts     = params['profile2counts']

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

    organisms = sorted(set(src2org.values()))
    _org2weight = tools.map_genome2weight()
    _total_weight = sum(_org2weight[file_summary.org] for file_summary in file_summaries)
    worksheet.merge_range(top_border, 0, top_border, 10, 'Organisms: %d,  Total weight: %f'%\
                          (len(organisms), _total_weight), title_format)
    top_border += 1

    worksheet.merge_range(top_border, 0, top_border, 10, ' '.join(organisms))
    top_border += 2

    # Writing the class distribution. The lestmost columns
    cur_top_border = top_border
    worksheet.write_row(cur_top_border, left_border, ['Weight','Class definition'], header_format)
    worksheet.set_column(left_border+1, left_border+1, 30)
    cur_top_border += 1
    worksheet.write_row(cur_top_border, left_border+1, ['Neighborhood'], header_format)
    cur_top_border += 2

    for (class_name, occurence) in sorted(class2counts.items(), key=itemgetter(1), reverse=True):
        worksheet.write_row(cur_top_border, left_border, [occurence, class_name])
        cur_top_border += 1

    cur_top_border += 2
    worksheet.write_row(cur_top_border, left_border+1, ['Neighborhood+flanks'], header_format)
    cur_top_border += 2

    for (class_name, occurence) in sorted(class2counts_flank.items(), key=itemgetter(1), reverse=True):
        worksheet.write_row(cur_top_border, left_border, [occurence, class_name])
        cur_top_border += 1

    # cur_top_border += 2
    # worksheet.write_row(cur_top_border, left_border, ['Code', 'Weight', 'Definition'], header_format)
    # cur_top_border += 2
    #
    # for (k, v) in sorted(profile2counts.items(), key= itemgetter(1), reverse=True):
    #     worksheet.write_row(cur_top_border, 0, [k, v, profile2def[k] if k in profile2def else '-'])
    #     cur_top_border += 1

    left_border = 5
    # Starting to write the data file-wise.
    for file_summary in file_summaries:

        cur_kplets = file_summary.kplets
        nbr = file_summary.neighborhood

        if not nbr.flank_extension:
            nbr.extend_flanks(10, os.path.join(gv.pty_data_path, file_summary.org, "%s.pty" % file_summary.src), gid2arcog_cdd)

        cur_top_border = top_border
        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, "%s %s" % (file_summary.org, file_summary.src), header_format)
        cur_top_border += 1
        worksheet.write_row(cur_top_border, left_border, column_names, header_format)
        cur_top_border += 2

        cur_top_border += 31 - len(nbr.genes)

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
            worksheet.set_column(left_border+row_len-1, left_border+row_len-1, 30)
            cur_top_border += 1
        cur_top_border += 2
        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, "Kplets:")
        cur_top_border += 1

        worksheet.write_row(cur_top_border, left_border, ["Id", "k", "Weight", "Count"])
        worksheet.merge_range(cur_top_border,left_border+4,cur_top_border,left_border+5, "Codes")
        cur_top_border += 1
        worksheet.write_row(cur_top_border, left_border, ["Total", "", file_summary.weight, file_summary.count])
        cur_top_border += 1

        for kplet in sorted(cur_kplets, reverse=True):
            worksheet.write_row(cur_top_border, left_border, [kplet.id, kplet.k, kplet.weight, kplet.count])
            worksheet.merge_range(cur_top_border,left_border+4,cur_top_border,left_border+5, " ".join(kplet.codes))
            cur_top_border += 1
        left_border += row_len + 1

    worksheet_2 = workbook.add_worksheet()

    top_border = 4
    worksheet_2.write_row(top_border, 0, ['Class','','Member arCOGS,CDD'], header_format)
    top_border += 2
    for (class_name, occurence) in sorted(class2counts.items(), key=itemgetter(1), reverse=True):

        profiles = class2profiles[class_name]
        profile_definitions = [profile2def[profile] for profile in profiles if profile in profile2def]

        worksheet_2.write_row(top_border, 0, [class_name, ' ']+profiles)
        top_border += 1
        worksheet_2.write_row(top_border, 0, [' ', ' ']+profile_definitions+[' '])
        top_border += 2
    workbook.close()



def write_flanking_count_xls(params):

    xls_file = params['xls_file_name']
    flank_counts = params['flank_counts']
    profile2def = params['profile2def']
    title_string = params['title']
    target_profiles = params['target_profiles']

    workbook = x.Workbook(xls_file)
    worksheet = workbook.add_worksheet()

    title_format = workbook.add_format()
    title_format.set_font_size(14)
    title_format.set_bold()
    title_format.set_align('left')

    top_border = 0

    worksheet.merge_range(0, 0, 0, 10, 'Flanking genes distribution for: '+title_string, title_format)
    top_border += 1

    worksheet.write_row(top_border, 0, ["Type","Code", "Weight", "Description"], title_format)
    top_border += 1

    for (k, v) in sorted(flank_counts.items(), key= itemgetter(1), reverse=True):
        _type = "seed" if k in target_profiles else "not-seed"
        worksheet.write_row(top_border, 0, [_type, k, v, profile2def[k] if k in profile2def else '-'])
        top_border += 1

    workbook.close()


def write_to_xls_raw_kplet(params):

    xls_file        = params['xls_file_name']
    src2org         = params['src2org']
    file_summaries  = params['file_summaries']
    target_profiles = params['target_profiles']
    profile2def     = params['profile2def']
    gid2arcog_cdd   = params['gid2arcog_cdd']

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

    # worksheet.merge_range(0, 0, 0, 10, 'Community: ' + ' '.join(community), title_format)
    top_border += 1

    organisms = sorted(set(src2org.values()))
    _org2weight = tools.map_genome2weight()
    _total_weight = sum(_org2weight[file_summary.org] for file_summary in file_summaries)
    worksheet.merge_range(top_border, 0, top_border, 10, 'Organisms: %d,  Total weight: %f'%\
                          (len(organisms), _total_weight), title_format)
    top_border += 1

    worksheet.merge_range(top_border, 0, top_border, 10, ' '.join(organisms))
    top_border += 2

    # Starting to write the data file-wise.
    ind = 0
    for file_summary in file_summaries:

        cur_kplets = file_summary.kplets
        nbr = file_summary.neighborhood

        community = set()
        for kplet in cur_kplets:
            community.update(kplet.codes)

        if not nbr.flank_extension:
            nbr.extend_flanks(10, os.path.join(gv.pty_data_path, file_summary.org, "%s.pty" % file_summary.src), gid2arcog_cdd)

        cur_top_border = top_border
        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, "%s %s %f" % (file_summary.org, file_summary.src, file_summary.self_weight), header_format)
        cur_top_border += 1
        worksheet.write_row(cur_top_border, left_border, column_names, header_format)
        cur_top_border += 2

        if nbr.genes[0].tag == 'neighborhood':
            cur_top_border += 31 - len(nbr.genes)

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
            worksheet.set_column(left_border+row_len-1, left_border+row_len-1, 30)
            cur_top_border += 1
        cur_top_border = 41
        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, "Kplets:")
        cur_top_border += 1

        worksheet.write_row(cur_top_border, left_border, ["Id", "k", "Weight", "Count"])
        worksheet.merge_range(cur_top_border,left_border+4,cur_top_border,left_border+5, "Codes")
        cur_top_border += 1
        worksheet.write_row(cur_top_border, left_border, ["Total", "", file_summary.kplets_weight, file_summary.count])
        cur_top_border += 1

        for kplet in sorted(cur_kplets, reverse=True):
            worksheet.write_row(cur_top_border, left_border, [kplet.id, kplet.k, kplet.weight, kplet.count])
            worksheet.merge_range(cur_top_border,left_border+4,cur_top_border,left_border+5, " ".join(kplet.codes))
            cur_top_border += 1
        left_border += row_len + 1

        ind += 1


def write_to_xls_extra_search(xls_file, org2src, src2blocks, kplet):

    kplet_codes = kplet.codes.difference(bacteria_target_profiles)

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

    worksheet.merge_range(top_border, 0, top_border, 10, 'Kplet id: %d'%kplet.id, header_format)
    top_border += 1
    for code in kplet.codes:
        worksheet.write_row(top_border, left_border, [code], header_format)
        worksheet.merge_range(top_border, left_border+1, top_border, 10, profile2def[code])
        top_border += 1

    top_border += 2

    organisms = sorted(org2src.keys())

    num_loci = sum([len(blocks) for blocks in src2blocks.values()])
    num_org = len(org2src)

    worksheet.write_row(top_border, left_border, ["Loci:%d"%num_loci])
    top_border += 1
    worksheet.write_row(top_border, left_border, ["Organisms:%d"%num_org])
    top_border += 2

    for org in organisms:

        cur_top_border = top_border
        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, org, header_format)

        cur_top_border += 2

        for src in org2src[org]:
            worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, src, header_format)
            cur_top_border += 2

            for block in src2blocks[src]:

                worksheet.write_row(cur_top_border, left_border, column_names, header_format)
                cur_top_border += 2

                for gene in block:
                    cur_cogid = gene.cogid
                    if cur_cogid.intersection(bacteria_target_profiles):
                        data_format = target_format_neighborhood if gene.tag == 'neighborhood' else target_format
                    elif cur_cogid.intersection(kplet_codes):
                        data_format = kplet_format_neighborhood if gene.tag == 'neighborhood' else kplet_format
                    else:
                        data_format = workbook.add_format()
                        if gene.tag == 'neighborhood':
                            data_format.set_bg_color('#c4bdbd')

                    if cur_cogid in ["", "-", None, []]:
                        cur_def = ""
                    else:
                        if len(cur_cogid) > 0:
                            cur_def = []
                            for k in cur_cogid:
                                if k in profile2def:
                                    cur_def.append(profile2def[k])
                                else:
                                    cur_def.append("")
                            cur_def = " | ".join(cur_def)

                    data_raw = [gene.gid, gene.pFrom, gene.pTo, gene.strand, " ".join(list(gene.cogid)), cur_def]
                    worksheet.write_row(cur_top_border, left_border, data_raw, data_format)
                    worksheet.write_row(cur_top_border, left_border+row_len, [" "])
                    worksheet.set_column(left_border+row_len-1, left_border+row_len-1, 30)
                    cur_top_border += 1

        left_border += row_len + 1


def write_to_xls_wgs_kplets(args):

    xls_file_name            = args.xls_file_name
    file_summaries           = args.file_summaries
    organisms                = args.organisms
    profile2def              = args.profile_code2def
    crispr_type2files        = args.crispr_type2files
    local_bf_kplet2count     = args.local_bf_kplet2count
    local_af_kplet2count     = args.local_af_kplet2count
    initial_length           = args.initial_length
    wgs_profile2count_bf     = args.wgs_profile2count_bf
    wgs_profile2count_af     = args.wgs_profile2count_af
    local_profile2count_bf   = args.local_profile2count_bf
    local_profile2count_af   = args.local_profile2count_af

    workbook = x.Workbook(xls_file_name)
    worksheet = workbook.add_worksheet()

    column_names = ['GI', 'From', 'To', 'Strand', 'Contig', 'CDD', 'Gene name', 'Definition']
    row_len = len(column_names)

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

    # worksheet.merge_range(0, 0, 0, 10, 'Community: ' + ' '.join(community), title_format)
    top_border += 1

    worksheet.merge_range(top_border, 0, top_border, 10, 'Organisms: %d,  Loci: %d'%\
                          (len(organisms), len(file_summaries)), title_format)
    top_border += 2

    for k in crispr_type2files:
        crispr_type2files[k] = len(crispr_type2files[k])

    crispr_type_summary = ",\t\t\t".join(["%s: %d" % (k,v) for k, v in sorted(crispr_type2files.items(), key=itemgetter(1), reverse=True)])

    worksheet.merge_range(top_border, 0, top_border, 10, 'System types:' + crispr_type_summary, title_format)
    top_border += 1
    worksheet.merge_range(top_border, 0, top_border, 10, 'Filtered out: %d'%(initial_length - len(file_summaries)))
    top_border += 2

    max_len = max([len(fs.genes) for fs in file_summaries])


    ind = 0
    conserved_kplets = {}
    for file_summary in file_summaries:
        [tools.update_dictionary(conserved_kplets, kplet.id, 1) for kplet in file_summary.kplets]
    conserved_kplet_ids = [k for (k,v) in conserved_kplets.items() if v > len(file_summaries)/2]

    # Starting to write the data file-wise.
    for file_summary in file_summaries:

        kplet_codes = set()
        [kplet_codes.update(kplet.codes) for kplet in file_summary.kplets if kplet.id in conserved_kplet_ids]

        cur_top_border = top_border
        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, "%s\t%s\t%s" % (file_summary.org, file_summary.src, file_summary.file_name), header_format)
        cur_top_border += 1
        worksheet.write_row(cur_top_border, left_border, column_names, header_format)
        cur_top_border += 2

        cur_top_border += max_len/2 - len(file_summary.genes)/2

        for gene in file_summary.genes:
            cur_cogid = set(gene.cogid.split(','))

            if gene.is_seed:
                data_format = target_format
            elif cur_cogid.intersection(kplet_codes):
            # if cur_cogid.intersection(kplet_codes):
                data_format = kplet_format
            else:
                data_format = workbook.add_format()

            if cur_cogid in [[""], ["-"], [], None]:
                cur_def = ""
            else:
                if len(cur_cogid) > 0:
                    cur_def = []
                    for k in cur_cogid:
                        if k in profile2def:
                            cur_def.append(profile2def[k])
                        else:
                            cur_def.append("")
                    cur_def = " | ".join(cur_def)

            data_raw = [gene.gid, gene.pFrom, gene.pTo, gene.strand, gene.src, gene.cogid, gene.gene_name, cur_def]
            worksheet.write_row(cur_top_border, left_border, data_raw, data_format)
            worksheet.write_row(cur_top_border, left_border+row_len, [" "])
            worksheet.set_column(left_border+row_len-1, left_border+row_len-1, 30)
            cur_top_border += 1


        cur_top_border = top_border + max_len + 4

        if file_summary.cluster == 'singleton':
            cluster_raw = "Represents: Singleton"
        else:
            cluster_raw = "Represents: Cluster size %d, local %d" % (len(file_summary.cluster.files), file_summary.cluster_local_count)
        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, cluster_raw)
        cur_top_border += 2

        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, "Kplets:")
        cur_top_border += 1

        worksheet.write_row(cur_top_border, left_border, ["Id", "Gl cnt", "Lc cnt bf", "Lc cnt af"])
        worksheet.merge_range(cur_top_border,left_border+4,cur_top_border,left_border+5, "Codes")
        cur_top_border += 1

        fs_kplets = sorted(file_summary.kplets, key = lambda x: x.count, reverse=True)
        fs_kplets = [kplet for kplet in fs_kplets if kplet.id in conserved_kplet_ids]

        for kplet in fs_kplets[:10]:

            if kplet.id not in conserved_kplet_ids:
                continue

            worksheet.write_row(cur_top_border, left_border, [kplet.id, kplet.count, local_bf_kplet2count[kplet.id], local_af_kplet2count[kplet.id]])
            worksheet.merge_range(cur_top_border,left_border+4,cur_top_border,left_border+5, " ".join(kplet.codes))
            cur_top_border += 1

        if len(fs_kplets) > 10:

            worksheet.merge_range(cur_top_border,left_border,cur_top_border,left_border+5, " ... ")
            cur_top_border += 1
            worksheet.merge_range(cur_top_border,left_border,cur_top_border,left_border+5, " ... ")
            cur_top_border += 1
            worksheet.merge_range(cur_top_border,left_border,cur_top_border,left_border+5, " ... ")
            cur_top_border += 1
            worksheet.merge_range(cur_top_border,left_border,cur_top_border,left_border+5, "%d kplets more"%(len(file_summary.kplets) - 10))
            cur_top_border += 1

        left_border += row_len + 1

        ind += 1

    local_profile_count = dict()
    local_kplet_member_count = dict()

    for file_summary in file_summaries:
        _codes = set()
        [_codes.update(kplet.codes) for kplet in file_summary.kplets]
        for _code in _codes:
            tools.update_dictionary(local_kplet_member_count, _code, 1)
        for gene in file_summary.genes:
            for code in gene.cogid.split(','):
                tools.update_dictionary(local_profile_count, code, 1)

    worksheet2 = workbook.add_worksheet()
    cur_top_border = 0

    worksheet2.write_row(cur_top_border, 0, ["Profile", "Local af", "Loacl bf", "Global af", "Global bf", "Description"], header_format)
    cur_top_border += 1

    worksheet2.set_column(0, 4, 15)
    worksheet2.set_column(5, 5, 50)

    for code, local_count in sorted(local_profile2count_af.items(), key=itemgetter(1), reverse=True):

        if code == "Unidentified":continue

        row = [code,
               local_count,
               local_profile2count_bf[code],
               wgs_profile2count_af[code] if code in wgs_profile2count_af else "-",
               wgs_profile2count_bf[code],
               profile2def[code],
               " "]

        worksheet2.write_row(cur_top_border, 0, row)
        cur_top_border += 1

    return crispr_type_summary


def write_to_summary_file(args):

    worksheet                = args.worksheet
    #cur_worksheet            = args.cur_worksheet
    kplet_list               = args.kplet_list
    ind                      = args.ind
    local_af_kplet2count     = args.local_af_kplet2count
    local_bf_kplet2count     = args.local_bf_kplet2count
    crispr_type_summary      = args.crispr_type_summary

    _id2kplet = {kplet.id: kplet for kplet in kplet_list.kplets}
    sorted_kplet2count = sorted(local_af_kplet2count.items(), key=itemgetter(1), reverse=True)
    kplets_summary = ";\t".join(["-".join(_id2kplet[kplet_id].codes) + ":" +str(cnt) for (kplet_id, cnt) in sorted_kplet2count])

    worksheet.write_row(ind+1, 0, ['%d.xls' % ind, crispr_type_summary, kplets_summary])

    #cur_worksheet.write_row(0, 0, ["Id", "Gl cnt", "Lc cnt bf", "Lc cnt af", "Codes"])
    #tmp_ind = 1

    #for kplet_id, cnt in sorted_kplet2count:
    #    codes = " ".join(_id2kplet[kplet_id].codes)
    #    row = [kplet_id, _id2kplet[kplet_id].count, local_bf_kplet2count[kplet_id], local_af_kplet2count[kplet_id], codes]
    #    cur_worksheet.write_row(tmp_ind, 0, row)
    #    tmp_ind += 1


def write_to_xls_generic_kplets(args):

    xls_file_name            = args.xls_file_name
    file_summaries           = args.file_summaries
    organisms                = args.organisms
    total_weight             = args.weight
    profile2def              = args.profile_code2def
    local_bf_kplet2count     = args.local_bf_kplet2count
    local_af_kplet2count     = args.local_af_kplet2count
    local_profile2count_bf   = args.local_profile2count_bf
    local_profile2count_af   = args.local_profile2count_af

    workbook = x.Workbook(xls_file_name)
    worksheet = workbook.add_worksheet()

    column_names = ['GI', 'From', 'To', 'Strand', 'Contig', 'CDD', 'Definition']
    row_len = len(column_names)

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

    # worksheet.merge_range(0, 0, 0, 10, 'Community: ' + ' '.join(community), title_format)
    top_border += 1

    worksheet.merge_range(top_border, 0, top_border, 10, 'Organisms: weight: %f count: %d,  Loci: %d'%\
                          (total_weight ,len(organisms), len(file_summaries)), title_format)
    top_border += 2

    max_len = max([len(fs.genes) for fs in file_summaries])

    ind = 0
    conserved_kplets = {}
    for file_summary in file_summaries:
        [tools.update_dictionary(conserved_kplets, kplet.id, 1) for kplet in file_summary.kplets]

    if len(file_summaries)>20:
        conserved_kplet_ids = [k for (k,v) in conserved_kplets.items() if v > len(file_summaries)/10]
    else:
        conserved_kplet_ids = [k for (k,v) in conserved_kplets.items() if v > 2]

    tmp = set()
    # Starting to write the data file-wise.
    for file_summary in file_summaries:
        tmp.update([file_summary.org])
        kplet_codes = set()
        [kplet_codes.update(kplet.codes) for kplet in file_summary.kplets if kplet.id in conserved_kplet_ids]

        cur_top_border = top_border
        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, "%s\t%s\t%s" % (file_summary.org, file_summary.src, file_summary.file_name), header_format)
        cur_top_border += 1
        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, "type: %s" % (file_summary.cas_type), header_format)
        cur_top_border += 1
        worksheet.write_row(cur_top_border, left_border, column_names, header_format)
        cur_top_border += 2

        cur_top_border += max_len/2 - len(file_summary.genes)/2

        for gene in file_summary.genes:
            cur_cogid = set(gene.cogid.split(','))

            if hasattr(gene,'is_seed') and gene.is_seed:
                data_format = target_format
            elif cur_cogid.intersection(kplet_codes):
                data_format = kplet_format
            else:
                data_format = workbook.add_format()

            if cur_cogid in [[""], ["-"], [], None]:
                cur_def = ""
            else:
                if len(cur_cogid) > 0:
                    cur_def = []
                    for k in cur_cogid:
                        if k in profile2def:
                            cur_def.append(profile2def[k])
                        else:
                            cur_def.append("")
                    cur_def = " | ".join(cur_def)

            data_raw = [gene.gid, gene.pFrom, gene.pTo, gene.strand, gene.src, gene.cogid, cur_def]
            worksheet.write_row(cur_top_border, left_border, data_raw, data_format)
            worksheet.write_row(cur_top_border, left_border+row_len, [" "])
            worksheet.set_column(left_border+row_len-1, left_border+row_len-1, 30)
            cur_top_border += 1

        cur_top_border = top_border + max_len + 4

        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, "Kplets:")
        cur_top_border += 1

        worksheet.write_row(cur_top_border, left_border, ["Id", "Gl cnt", "Lc cnt bf", "Lc cnt af"])
        worksheet.merge_range(cur_top_border,left_border+4,cur_top_border,left_border+5, "Codes")
        cur_top_border += 1

        fs_kplets = sorted(file_summary.kplets, key = lambda x: x.count, reverse=True)
        # fs_kplets = [kplet for kplet in fs_kplets if kplet.id in conserved_kplet_ids]

        for kplet in fs_kplets[:10]:

            # if kplet.id not in conserved_kplet_ids:
            #     continue

            worksheet.write_row(cur_top_border, left_border, [kplet.id, kplet.count, local_bf_kplet2count[kplet.id], local_af_kplet2count[kplet.id]])
            worksheet.merge_range(cur_top_border,left_border+4,cur_top_border,left_border+5, " ".join(kplet.codes))
            cur_top_border += 1

        if len(fs_kplets) > 10:

            worksheet.merge_range(cur_top_border,left_border,cur_top_border,left_border+5, " ... ")
            cur_top_border += 1
            worksheet.merge_range(cur_top_border,left_border,cur_top_border,left_border+5, " ... ")
            cur_top_border += 1
            worksheet.merge_range(cur_top_border,left_border,cur_top_border,left_border+5, " ... ")
            cur_top_border += 1
            worksheet.merge_range(cur_top_border,left_border,cur_top_border,left_border+5, "%d more kplets "%(len(file_summary.kplets) - 10))
            cur_top_border += 1

        left_border += row_len + 1

        ind += 1

    local_profile_count = dict()
    local_kplet_member_count = dict()

    for file_summary in file_summaries:
        _codes = set()
        [_codes.update(kplet.codes) for kplet in file_summary.kplets]
        for _code in _codes:
            tools.update_dictionary(local_kplet_member_count, _code, 1)
        for gene in file_summary.genes:
            for code in gene.cogid.split(','):
                tools.update_dictionary(local_profile_count, code, 1)

    worksheet2 = workbook.add_worksheet()
    cur_top_border = 0

    # worksheet2.write_row(cur_top_border, 0, ["Profile", "Local af", "Loacl bf", "Global af", "Global bf", "Description"], header_format)
    worksheet2.write_row(cur_top_border, 0, ["Profile", "Local af", "Description"], header_format)
    cur_top_border += 1

    worksheet2.set_column(0, 4, 15)
    worksheet2.set_column(5, 5, 50)

    for code, local_count in sorted(local_profile2count_af.items(), key=itemgetter(1), reverse=True):

        if code == "Unidentified":continue

        row = [code,
               local_count,
               # local_profile2count_bf[code],
               # wgs_profile2count_af[code] if code in wgs_profile2count_af else "-",
               # wgs_profile2count_bf[code],
               profile2def[code] if code in profile2def else ""]

        worksheet2.write_row(cur_top_border, 0, row)
        cur_top_border += 1


def write_to_xls_generic_loci(args):

    xls_file_name  = args['xls_file_name']
    loci           = args['loci']
    total_weight   = args['weight']
    profile2def    = args['profile_code2def']
    file2crispr_type = args['file2crispr_type']
    feature_labels = args['feature_labels']
    local_profile2weight = args['local_profile2weight']
    global_profile2weight = args['global_profile2weight']

    workbook = x.Workbook(xls_file_name)
    worksheet = workbook.add_worksheet()

    column_names = ['GI', 'From', 'To', 'Strand', 'Contig', 'CDD', 'Definition']
    row_len = len(column_names)

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

    # worksheet.merge_range(0, 0, 0, 10, 'Community: ' + ' '.join(community), title_format)
    top_border += 1

    organisms = set([locus.organism for locus in loci])

    worksheet.merge_range(top_border, 0, top_border, 10, 'Total weight: %f count: %d. Unique organisms: %d '%\
                          (total_weight ,len(loci), len(organisms)), title_format)
    top_border += 2

    max_len = max([len(locus.genes) for locus in loci])

    # Starting to write the data file-wise.
    for locus in loci:
        file_name = os.path.basename(locus.file_name)
        cur_top_border = top_border
        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, "%s\t%s\t%s" %
                              (locus.organism, locus.source, file_name),
                              header_format)
        cur_top_border += 1

        cas_type = ";".join(file2crispr_type[file_name]) if file_name in file2crispr_type else "NA"

        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1,
                              "Weight: %f, type: %s" % (gnm2weight[locus.organism],cas_type),
                              header_format)
        cur_top_border += 1
        worksheet.write_row(cur_top_border, left_border, column_names, header_format)
        cur_top_border += 2

        cur_top_border += max_len/2 - len(locus.genes)/2

        for gene in locus.genes:
            cur_cogid = set(gene.cogid.split(',')) if ',' in gene.cogid else set(gene.cogid.split())

            if cur_cogid.intersection(feature_labels):
                data_format = kplet_format
            else:
                data_format = workbook.add_format()

            cur_def = ""
            if cur_cogid not in [[""], ["-"], [], None] and len(cur_cogid) > 0:
                cur_def = []
                for k in cur_cogid:
                    if k in profile2def:
                        cur_def.append(profile2def[k])
                    else:
                        cur_def.append("")
                cur_def = " | ".join(cur_def)

            data_raw = [gene.gid, gene.pFrom, gene.pTo, gene.strand, gene.src, gene.cogid, cur_def]
            worksheet.write_row(cur_top_border, left_border, data_raw, data_format)
            worksheet.write_row(cur_top_border, left_border+row_len, [" "])
            worksheet.set_column(left_border+row_len-1, left_border+row_len-1, 30)
            cur_top_border += 1

        left_border += row_len + 1

    worksheet2 = workbook.add_worksheet()

    worksheet2.write_row(0, 0, ["Profile", "Local weight", "Global weight", "Desctiption"], header_format)

    cnt = 1
    for profile, local in sorted(local_profile2weight.items(), key=itemgetter(1), reverse=True):

        if profile.isdigit():
            continue

        global_weight = global_profile2weight[profile] if profile in global_profile2weight else "-"
        profile_def = profile2def[profile] if profile in profile2def else "-"
        worksheet2.write_row(3+cnt, 0, [profile, local, global_weight, profile_def])
        cnt += 1


def write_to_xls_generic_loci_cas4(args):

    xls_file_name  = args['xls_file_name']
    loci           = args['loci']
    cluster_packs  = args['clusters']
    singles        = args['singles']
    profile2def    = args['profile_code2def']
    feature_labels = args['feature_labels']

    workbook = x.Workbook(xls_file_name)
    worksheet = workbook.add_worksheet()

    column_names = ['GI', 'From', 'To', 'Strand', 'Contig', 'Cluster id', 'Genes', 'CDD', 'Definition']
    row_len = len(column_names)

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

    # worksheet.merge_range(0, 0, 0, 10, 'Community: ' + ' '.join(community), title_format)
    top_border += 1

    organisms = set([locus.organism for locus in loci])

    worksheet.merge_range(top_border, 0, top_border, 10, 'count: %d. Unique organisms: %d '%\
                          (len(loci), len(organisms)), title_format)
    top_border += 2

    max_len = max([len(locus.genes) for locus in loci])

    file_name2locus = {}
    for locus in loci:
        _file_name = os.path.basename(locus.file_name)
        file_name2locus[_file_name] = locus

    # This part should be paid attention to. The content of cluster_packs can change
    # depending on the function used to produce cluster_packs
    # clusters = []
    # for cluster_pack in cluster_packs:
    #     clusters.append(cluster_pack[1])

    clusters = cluster_packs

    clusters.sort(key=lambda x: len(x), reverse=True)

    for cluster in clusters:

        cluster_size = len(cluster)
        cluster_ind = 1

        for cl_ind in cluster:

            # locus = file_name2locus[cur_file]
            locus = loci[cl_ind]
            cur_file = os.path.basename(locus.file_name)

            cur_top_border = top_border
            worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1,
                                  "%s\t%s\t%s Cluster: %d / %d" %
                                  (locus.organism, locus.source, cur_file, cluster_ind, cluster_size),
                                  header_format)
            cur_top_border += 1

            worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1,
                                  "type: %s" % locus.crispr_type,
                                  header_format)
            cur_top_border += 1
            worksheet.write_row(cur_top_border, left_border, column_names, header_format)
            cur_top_border += 2

            cur_top_border += max_len/2 - len(locus.genes)/2

            for gene in locus.genes:
                # cur_cogid = set(gene.cogid.split(',')) if ',' in gene.cogid else set(gene.cogid.split())

                if gene.is_seed:
                    data_format = target_format
                elif gene.cluster_id in feature_labels:
                    data_format = kplet_format
                else:
                    data_format = workbook.add_format()

                cdd_profiles = gene.cogid
                gene_names = gene.gene_name
                first_profile = cdd_profiles.split(',')[0].strip()
                cur_def = profile2def[first_profile] if first_profile in profile2def else ""

                data_raw = [gene.gid, gene.pFrom, gene.pTo, gene.strand, gene.src, gene.cogid, gene_names, cdd_profiles, cur_def]
                worksheet.write_row(cur_top_border, left_border, data_raw, data_format)
                worksheet.write_row(cur_top_border, left_border+row_len, [" "])
                worksheet.set_column(left_border+row_len-1, left_border+row_len-1, 30)
                cur_top_border += 1

            left_border += row_len + 1
            cluster_ind += 1


    for ind in singles:

        # locus = file_name2locus[os.path.basename(cur_file)]
        locus = loci[ind]
        cur_file = os.path.basename(locus.file_name)

        cur_top_border = top_border
        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1,
                              "%s\t%s\t%s Cluster: singleton" %
                              (locus.organism, locus.source, cur_file),
                              header_format)
        cur_top_border += 1

        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1,
                              "type: %s" % locus.crispr_type,
                              header_format)
        cur_top_border += 1
        worksheet.write_row(cur_top_border, left_border, column_names, header_format)
        cur_top_border += 2

        cur_top_border += max_len/2 - len(locus.genes)/2

        for gene in locus.genes:
            cur_cogid = set(gene.cogid.split(',')) if ',' in gene.cogid else set(gene.cogid.split())

            if gene.is_seed:
                data_format = target_format
            elif cur_cogid.intersection(feature_labels):
                data_format = kplet_format
            else:
                data_format = workbook.add_format()

            cdd_profiles = gene.cogid
            gene_names = gene.gene_name
            first_profile = cdd_profiles.split(':')[0].strip()
            cur_def = profile2def[first_profile] if first_profile in profile2def else ""

            data_raw = [gene.gid, gene.pFrom, gene.pTo, gene.strand, gene.src, gene.cogid, gene_names, cdd_profiles, cur_def]
            worksheet.write_row(cur_top_border, left_border, data_raw, data_format)
            worksheet.write_row(cur_top_border, left_border+row_len, [" "])
            worksheet.set_column(left_border+row_len-1, left_border+row_len-1, 30)
            cur_top_border += 1

        left_border += row_len + 1



def write_to_xls_generic_loci_cas4_vertical(args):

    xls_file_name  = args['xls_file_name']
    loci           = args['loci']
    cluster_packs  = args['clusters']
    singles        = args['singles']
    profile2def    = args['profile_code2def']
    feature_labels = args['feature_labels']

    workbook = x.Workbook(xls_file_name)
    worksheet = workbook.add_worksheet()

    column_names = ['GI', 'From', 'To', 'Strand', 'Contig', 'Cluster id', 'Genes', 'CDD', 'Definition']
    row_len = len(column_names)

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

    # worksheet.merge_range(0, 0, 0, 10, 'Community: ' + ' '.join(community), title_format)
    top_border += 1

    organisms = set([locus.organism for locus in loci])

    worksheet.merge_range(top_border, 0, top_border, 10, 'count: %d. Unique organisms: %d '%\
                          (len(loci), len(organisms)), title_format)
    top_border += 2

    # clusters = []
    clusters = cluster_packs

    # for cluster_pack in cluster_packs:
    #
    #     clusters.append(cluster_pack[0])

    clusters.sort(key=lambda x: len(x), reverse=True)

    top_border += 3

    for cluster in clusters:

        cluster_size = len(cluster)
        cluster_ind = 1

        max_len = max([len(loci[cl_ind].genes) for cl_ind in cluster])
        
        worksheet.merge_range(top_border, left_border, top_border, left_border + row_len-1,
                                  "Group size: %d" % cluster_size, header_format)
        
        for cl_ind in cluster:

            cur_top_border = top_border + 2

            locus = loci[cl_ind]
            cur_file = os.path.basename(locus.file_name)

            worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1,
                                  "%s\t%s\t%s, index:%d" % (locus.organism, locus.source, cur_file, cluster_ind),
                                  header_format)
            cur_top_border += 1

            worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1,
                                  "type: %s" % locus.crispr_type,
                                  header_format)
            cur_top_border += 1
            worksheet.write_row(cur_top_border, left_border, column_names, header_format)
            cur_top_border += 2

            cur_top_border += max_len/2 - len(locus.genes)/2

            for gene in locus.genes:
                
                if gene.is_seed:
                    data_format = target_format
                elif gene.cluster_id in feature_labels:
                    data_format = kplet_format
                else:
                    data_format = workbook.add_format()

                cdd_profiles = gene.cogid
                gene_names = gene.gene_name
                first_profile = cdd_profiles.split(',')[0].strip()
                cur_def = profile2def[first_profile] if first_profile in profile2def else ""

                data_raw = [gene.gid, gene.pFrom, gene.pTo, gene.strand, gene.src, gene.cluster_id, gene_names, cdd_profiles, cur_def]
                worksheet.write_row(cur_top_border, left_border, data_raw, data_format)
                worksheet.write_row(cur_top_border, left_border+row_len, [" "])
                worksheet.set_column(left_border+row_len-1, left_border+row_len-1, 30)
                cur_top_border += 1

            left_border += row_len + 1
            cluster_ind += 1

        left_border = 0
        top_border += max_len + 10

    cur_top_border = top_border + 2
    worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1,
                                  "Sinlges: %d " % len(singles), header_format)

    top_border += 4

    max_len = max([len(loci[cl_ind].genes) for cl_ind in singles]) if singles else 0

    for ind in singles:

        locus = loci[ind]
        cur_file = os.path.basename(locus.file_name)

        cur_top_border = top_border 

        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1,
                              "%s\t%s\t%s Cluster: singleton" %
                              (locus.organism, locus.source, cur_file),
                              header_format)
        cur_top_border += 1

        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1,
                              "type: %s" % locus.crispr_type,
                              header_format)
        cur_top_border += 1
        worksheet.write_row(cur_top_border, left_border, column_names, header_format)
        cur_top_border += 2

        cur_top_border += max_len/2 - len(locus.genes)/2

        for gene in locus.genes:
            cur_cogid = set(gene.cogid.split(',')) if ',' in gene.cogid else set(gene.cogid.split())

            if gene.is_seed:
                data_format = target_format
            elif gene.cluster_id in feature_labels:
                data_format = kplet_format
            else:
                data_format = workbook.add_format()

            cdd_profiles = gene.cogid
            gene_names = gene.gene_name
            first_profile = cdd_profiles.split(',')[0].strip()
            cur_def = profile2def[first_profile] if first_profile in profile2def else ""

            data_raw = [gene.gid, gene.pFrom, gene.pTo, gene.strand, gene.src, gene.cluster_id, gene_names, cdd_profiles, cur_def]
            worksheet.write_row(cur_top_border, left_border, data_raw, data_format)
            worksheet.write_row(cur_top_border, left_border+row_len, [" "])
            worksheet.set_column(left_border+row_len-1, left_border+row_len-1, 30)
            cur_top_border += 1

        left_border += row_len + 1


def write_loci_to_tab_file_cas4(args):

    tab_file_name  = args['tab_file_name']
    loci           = args['loci']
    profile2def    = args['profile_code2def']

    outf = open(tab_file_name, 'w')

    column_names = ['GI', 'From', 'To', 'Strand', 'Contig', 'Cluster id', 'Genes', 'CDD']

    outf.write("#"+"\t".join(column_names)+"\n")

    for locus in loci:
        outf.write(locus.summary_line)
        for g in locus.genes:
            out_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
            outf.write(out_line % (g.gid, g.pFrom, g.pTo, g.strand, g.src, g.cluster_id, g.gene_name, g.cogid))

    outf.close()


def write_to_xls_generic_singleton_loci(args):

    xls_file_name  = args['xls_file_name']
    loci           = args['loci']
    profile2def    = args['profile_code2def']
    feature_labels = args['feature_labels']
    file2crispr_type = args['file2crispr_type']

    workbook = x.Workbook(xls_file_name)
    worksheet = workbook.add_worksheet()

    column_names = ['GI', 'From', 'To', 'Strand', 'Contig', 'CDD', 'Definition']
    row_len = len(column_names)

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

    # worksheet.merge_range(0, 0, 0, 10, 'Community: ' + ' '.join(community), title_format)
    top_border += 1

    max_len = max([len(locus.genes) for locus in loci])

    # Starting to write the data file-wise.
    for locus in loci:
        file_name = os.path.basename(locus.file_name)
        cur_top_border = top_border
        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, "%s\t%s\t%s" %
                              (locus.organism, locus.source, file_name),
                              header_format)
        cur_top_border += 1

        cas_type = ";".join(file2crispr_type[file_name]) if file_name in file2crispr_type else "NA"

        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1,
                              "Weight: %f, type: %s" % (gnm2weight[locus.organism],cas_type),
                              header_format)
        cur_top_border += 1
        worksheet.write_row(cur_top_border, left_border, column_names, header_format)
        cur_top_border += 2

        cur_top_border += max_len/2 - len(locus.genes)/2

        for gene in locus.genes:
            cur_cogid = set(gene.cogid.split(',')) if ',' in gene.cogid else set(gene.cogid.split())

            if cur_cogid.intersection(feature_labels):
                data_format = kplet_format
            else:
                data_format = workbook.add_format()

            cur_def = ""
            if cur_cogid not in [[""], ["-"], [], None] and len(cur_cogid) > 0:
                cur_def = []
                for k in cur_cogid:
                    if k in profile2def:
                        cur_def.append(profile2def[k])
                    else:
                        cur_def.append("")
                cur_def = " | ".join(cur_def)

            data_raw = [gene.gid, gene.pFrom, gene.pTo, gene.strand, gene.src, gene.cogid, cur_def]
            worksheet.write_row(cur_top_border, left_border, data_raw, data_format)
            worksheet.write_row(cur_top_border, left_border+row_len, [" "])
            worksheet.set_column(left_border+row_len-1, left_border+row_len-1, 30)
            cur_top_border += 1

        left_border += row_len + 1


def write_to_xls_baiticity_profile_set(args):

    xls_file            = args['xls_file_name']
    file_summaries      = args['file_summaries']
    clusters            = args['clusters']
    singles             = args['singles']
    target_profiles     = args['target_profiles']
    additional_profiles = args['additional_profiles']
    profile2def         = args['profile2def']
    profile2baiticity   = args['profile2baiticity']

    workbook = x.Workbook(xls_file)
    worksheet = workbook.add_worksheet()

    column_names = ['GI', 'From', 'To', 'Strand', 'CDD', 'Baiticity', 'Definition']
    row_len = len(column_names)

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

    default_format = workbook.add_format()

    default_format_neigborhood = workbook.add_format()
    default_format_neigborhood.set_bg_color('#c4bdbd')

    top_border = 0
    left_border = 0

    # worksheet.merge_range(0, 0, 0, 10, 'Community: ' + ' '.join(community), title_format)
    top_border += 1

    organisms = set()
    _unique_weight = 0

    file_name2file_summary = {}

    for file_summary in file_summaries:

        _file_name = os.path.basename(file_summary.file_name)

        file_name2file_summary[_file_name] = file_summary

        if not file_summary.org in organisms:
            organisms.update([file_summary.org])
            _unique_weight += file_summary.weight

    _total_weight = sum([file_summary.weight for file_summary in file_summaries])
    worksheet.merge_range(top_border, 0, top_border, 10, 'Organisms: %d, Unique weight: %f  Total weight: %f' % \
                          (len(organisms), _unique_weight, _total_weight), title_format)
    top_border += 1

    worksheet.merge_range(top_border, 0, top_border, 10, ' '.join(organisms))
    top_border += 2

    # Starting to write the data file-wise.
    ind = 0
    # for file_summary in file_summaries:
    for cluster in clusters:

        cluster_size = len(cluster)
        cluster_ind = 1

        for cur_file in cluster:

            file_summary = file_name2file_summary[cur_file]

            nbr = file_summary.neighborhood

            cur_top_border = top_border
            worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len - 1,
                                  "%s %s %f Cluster: %d / %d" %
                                  (file_summary.org, file_summary.src, file_summary.weight, cluster_ind, cluster_size),
                                  header_format)

            cur_top_border += 1
            worksheet.write_row(cur_top_border, left_border, column_names, header_format)
            cur_top_border += 2

            if nbr.genes[0].tag == 'neighborhood':
                cur_top_border += 31 - len(nbr.genes)

            if cur_top_border < 7:
                cur_top_border = 7

            for gene in nbr.genes:

                data_format = default_format

                _baiticity_terms = []

                if gene.tag == 'neighborhood':
                    data_format = default_format_neigborhood

                for cur_cogid in gene.cogid.split():

                    if cur_cogid in additional_profiles:
                        data_format = kplet_format_neighborhood if gene.tag == 'neighborhood' else kplet_format
                        break
                    elif cur_cogid in target_profiles:
                        data_format = target_format_neighborhood if gene.tag == 'neighborhood' else target_format
                        break

                for cur_cogid in gene.cogid.split():

                    if cur_cogid in target_profiles:
                        _baiticity_terms.append('bait')
                    elif cur_cogid in additional_profiles:
                        _str_baiticity = "%.2f" % float(profile2baiticity[cur_cogid])
                        _baiticity_terms.append(_str_baiticity)

                if gene.cogid in ["", "-", None]:
                    cur_def = ""
                else:
                    cur_cogid = gene.cogid.split()
                    if len(cur_cogid) > 0:
                        cur_def = []
                        for k in cur_cogid:
                            if k in profile2def:
                                cur_def.append(profile2def[k])
                            else:
                                cur_def.append("")
                        cur_def = " | ".join(cur_def)

                _baiticity = "|".join(_baiticity_terms) if _baiticity_terms else ' '

                data_raw = [gene.gid, gene.pFrom, gene.pTo, gene.strand, gene.cogid, _baiticity, cur_def]
                worksheet.write_row(cur_top_border, left_border, data_raw, data_format)
                worksheet.write_row(cur_top_border, left_border + row_len, [" "])
                worksheet.set_column(left_border + row_len - 1, left_border + row_len - 1, 30)
                cur_top_border += 1

            left_border += row_len + 1
            cluster_ind += 1


    for cur_file in singles:

        file_summary = file_name2file_summary[cur_file]

        nbr = file_summary.neighborhood

        cur_top_border = top_border
        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len - 1,
                              "%s %s %f Cluster: singleton" %
                              (file_summary.org, file_summary.src, file_summary.weight),
                              header_format)

        cur_top_border += 1
        worksheet.write_row(cur_top_border, left_border, column_names, header_format)
        cur_top_border += 2

        if nbr.genes[0].tag == 'neighborhood':
            cur_top_border += 31 - len(nbr.genes)

        if cur_top_border < 7:
            cur_top_border = 7

        for gene in nbr.genes:

            data_format = default_format

            _baiticity_terms = []

            if gene.tag == 'neighborhood':
                data_format = default_format_neigborhood

            for cur_cogid in gene.cogid.split():

                if cur_cogid in additional_profiles:
                    data_format = kplet_format_neighborhood if gene.tag == 'neighborhood' else kplet_format
                    break
                elif cur_cogid in target_profiles:
                    data_format = target_format_neighborhood if gene.tag == 'neighborhood' else target_format
                    break

            for cur_cogid in gene.cogid.split():

                if cur_cogid in target_profiles:
                    _baiticity_terms.append('bait')
                elif cur_cogid in additional_profiles:
                    _str_baiticity = "%.2f" % float(profile2baiticity[cur_cogid])
                    _baiticity_terms.append(_str_baiticity)

            if gene.cogid in ["", "-", None]:
                cur_def = ""
            else:
                cur_cogid = gene.cogid.split()
                if len(cur_cogid) > 0:
                    cur_def = []
                    for k in cur_cogid:
                        if k in profile2def:
                            cur_def.append(profile2def[k])
                        else:
                            cur_def.append("")
                    cur_def = " | ".join(cur_def)

            _baiticity = "|".join(_baiticity_terms) if _baiticity_terms else ' '

            data_raw = [gene.gid, gene.pFrom, gene.pTo, gene.strand, gene.cogid, _baiticity, cur_def]
            worksheet.write_row(cur_top_border, left_border, data_raw, data_format)
            worksheet.write_row(cur_top_border, left_border + row_len, [" "])
            worksheet.set_column(left_border + row_len - 1, left_border + row_len - 1, 30)
            cur_top_border += 1

        left_border += row_len + 1

        ind += 1



def write_to_xls_loci_plain(args):

    xls_file_name  = args['xls_file_name']
    loci           = args['loci']
    profile2def    = args['profile_code2def']
    feature_labels = args['feature_labels']

    workbook = x.Workbook(xls_file_name)
    worksheet = workbook.add_worksheet()

    column_names = ['GI', 'From', 'To', 'Strand', 'Contig', 'Cluster id', 'Genes', 'CDD', 'Definition']
    row_len = len(column_names)

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

    # worksheet.merge_range(0, 0, 0, 10, 'Community: ' + ' '.join(community), title_format)
    top_border += 1
    organisms = set([locus.organism for locus in loci])
    worksheet.merge_range(top_border, 0, top_border, 10, 'count: %d. Unique organisms: %d ' % \
                          (len(loci), len(organisms)), title_format)
    top_border += 2

    max_len = max([len(locus.genes) for locus in loci])

    # Starting to write the data file-wise.
    for locus in loci:
        file_name = os.path.basename(locus.file_name)
        cur_top_border = top_border
        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, "%s\t%s\t%s" %
                              (locus.organism, locus.source, file_name),
                              header_format)
        cur_top_border += 1

        worksheet.write_row(cur_top_border, left_border, column_names, header_format)
        cur_top_border += 2

        cur_top_border += max_len/2 - len(locus.genes)/2

        for gene in locus.genes:
            if gene.is_seed:
                data_format = target_format
            elif gene.cluster_id in feature_labels:
                data_format = kplet_format
            else:
                data_format = workbook.add_format()

            cdd_profiles = gene.cogid
            gene_names = gene.gene_name
            first_profile = cdd_profiles.split(',')[0].strip()
            cur_def = profile2def[first_profile] if first_profile in profile2def else ""

            data_raw = [gene.gid, gene.pFrom, gene.pTo, gene.strand, gene.src, gene.cluster_id, gene_names,
                        cdd_profiles, cur_def]
            worksheet.write_row(cur_top_border, left_border, data_raw, data_format)
            worksheet.write_row(cur_top_border, left_border + row_len, [" "])
            worksheet.set_column(left_border + row_len - 1, left_border + row_len - 1, 30)
            cur_top_border += 1

        left_border += row_len + 1
