import os
import sys
if sys.platform == 'darwin':
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))
elif sys.platform == 'linux2':
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))
import global_variables as gv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import shutil as sh
import xlsxwriter as x
from operator import itemgetter

import scores
import dendrogram

from lib.utils import tools as t
gnm2weight = t.map_genome2weight()
# file2org = {l.split()[0]:l.strip().split()[1] for l in open(os.path.join(gv.project_data_path,'cas1402/file2org.txt')).readlines()}
# file2crispr_type = {l.split('\t')[0]:l.strip().split('\t')[1].split(';') for l in open(os.path.join(gv.project_data_path,'cas1402/file2type.tab'))}

import lib.utils.reporting as r

from lib.db.generic import map_profiles_id2code_code2def
(_, profile_code2def) = map_profiles_id2code_code2def('cas')


def plot_block(block):

    _fname = block[0].strip()

    thresholds = []
    singles = []
    clusters = []
    mean_errors = []
    entropies = []

    for l in block[1:]:

        terms = l.strip().split('\t')

        thresholds.append(terms[0])
        singles.append(terms[1])
        clusters.append(terms[2])
        entropies.append(terms[3])

    thresholds = np.asarray(thresholds,dtype=np.float)
    singles = np.asarray(singles, dtype=np.int)
    clusters = np.asarray(clusters, dtype=np.int)
    entropies = 1000*np.asarray(entropies, dtype=np.float)

    plt.plot(thresholds, singles)
    plt.plot(thresholds, clusters)
    plt.plot(thresholds, entropies)

    _cnt = _fname.split('_')[1]
    _crispricity = _fname.split()[0].split('_')[2][:-4]
    _profiles = _fname.split()[1]

    plt.title("Occurrence:%s, CRISPRicity:%s, profiles: %s"%(_cnt, _crispricity, _profiles))
    plt.grid(True)
    plt.legend(['Singletons', 'Clusters', r'$ \left<\textit{I}\right> x(10^3) $'],loc='upper left')
    plt.xticks(thresholds, [str(t) for t in thresholds], rotation='vertical')
    plt.xlabel('Clustering thresholds')


def plot_results(data_file_name = 'results.txt', image_file_name='results.png'):

    plt.figure(figsize=(60,20))
    plt.rc('text', usetex=True)

    font = {'family': 'serif',
            'weight': 'bold',
            'size': 22}

    plt.rc('font', **font)

    block = []
    i = 1

    for l in open(data_file_name).readlines():
        if not block:
            block.append(l)
            continue

        if l.startswith('#'):
            plt.subplot(1,3,i)
            plt.tight_layout()
            plot_block(block)
            i += 1
            block = [l]
            continue

        block.append(l)

    plt.subplot(1,3,i)
    plot_block(block)

    plt.savefig(image_file_name)


def report_clustering_dot_product(loci, thresholds_pack, method, feature_labels):

    thr_occ, thr_crisp, cluster_thresholds = thresholds_pack

    M = scores.generate_dot_product_score_matrix(feature_labels, method, loci=loci)
    M += np.transpose(M)
    M = -1 * np.log(M)
    M[np.diag_indices_from(M)] = 0
    M[np.where(M==np.inf)] = 100

    reports_dir_base = os.path.join(gv.project_data_path, 'cas4/reports/')

    cluster2summary_file_path = os.path.join(gv.project_data_path, 'cas4/reports/cluster_summary.tab')

    for threshold in cluster_thresholds:

        repors_dir = reports_dir_base + 'dot_%s_%d_%.2f_%.2f'%(method, thr_occ, thr_crisp, threshold)
        # print "Thresholds:", thr_occ, thr_crisp, threshold
        # print repors_dir
        # if os.path.exists(repors_dir):
        #     sh.rmtree(repors_dir)
        # os.mkdir(repors_dir)

        singles, cluster_packs, entropies = dendrogram.classify_by_scores_cas4(M, threshold, loci)

        _local_thresholds_pack = (thr_occ, thr_crisp, threshold)

        generate_cluster_reports_cas4(cluster_packs,
                                      loci,
                                      repors_dir,
                                      feature_labels,
                                      method,
                                      _local_thresholds_pack)

        generate_cas4_gi_summary_file(singles, cluster_packs, loci, repors_dir, cluster2summary_file_path)



def report_clustering_jw(loci, thresholds):

    M = scores.generate_jackard_score_matrix(loci)

    M += np.transpose(M)
    M = -1 * np.log(M)
    M[np.diag_indices_from(M)] = 0
    M[np.where(M==np.inf)] = 50

    reports_dir_base = os.path.join(gv.project_data_path, 'cas1402/reports/')

    for threshold in thresholds:

        repors_dir = reports_dir_base + 'jw_%.2f' % threshold
        print repors_dir
        if os.path.exists(repors_dir):
            sh.rmtree(repors_dir)
        os.mkdir(repors_dir)

        singles, cluster_packs, sum_errors, entropies = dendrogram.classify_by_scores(M, threshold, loci)

        # generate_cluster_reports(cluster_packs, loci, repors_dir, feature_labels, method, _local_thresholds_pack)
        generate_jw_cluster_reports(cluster_packs, loci, repors_dir, threshold)


def generate_cluster_reports(cluster_packs, loci, reports_dir, feature_labels, method, thresholds_pack):

    if not feature_labels:
        local_features = True
    else:
        local_features = False

    thr_occ, thr_crisp, cluster_threshold = thresholds_pack

    summary_file = os.path.join(reports_dir,
                                'summary_%s_%d_%.2f_%.2f.xls' % (method, thr_occ, thr_crisp, cluster_threshold))

    workbook = x.Workbook(summary_file)
    worksheet = workbook.add_worksheet()

    header_format = workbook.add_format()
    header_format.set_font_size(12)
    header_format.set_bold()
    header_format.set_align('center')
    worksheet.set_column(4,5,50)
    worksheet.write_row(0, 0, ["File name", "Weight", "Loci", "Entropy", "systems weight", "systems count"], header_format)

    print "Generating report files"
    ind = 0

    weights = np.zeros(len(cluster_packs))
    entropies = np.zeros(len(cluster_packs))

    for outer_i in range(len(cluster_packs)):

        (cluster, type2count, type2weight, entropy) = cluster_packs[outer_i]

        ind += 1
        cl_files = [os.path.basename(loci[i].file_name) for i in cluster]

        weight = sum([gnm2weight[file2org[file]] for file in cl_files])

        weights[outer_i] = weight
        entropies[outer_i] = entropy

        crispr_cas_types_count = " ; ".join([k+":"+str(v) for (k,v) in sorted(type2count.items(), key=itemgetter(1), reverse=True)])
        crispr_cas_types_weight = " ; ".join([k+":"+str(v) for (k,v) in sorted(type2weight.items(), key=itemgetter(1), reverse=True)])

        xls_file_name = os.path.join(reports_dir, '%d.xls' % ind)

        worksheet.write_row(ind+1, 0, ['%d.xls'%ind,
                                       weight,
                                       len(cl_files),
                                       entropy,
                                       crispr_cas_types_weight,
                                       crispr_cas_types_count,
                                       " "])

        cl_loci = sorted([loci[_i] for _i in cluster], key = lambda x: gnm2weight[x.organism], reverse=True)

        local_profile2weight = {}
        for locus in cl_loci:
            for gene in locus.genes:
                for profile in gene.cogid.split(','):
                    t.update_dictionary(local_profile2weight, profile, gnm2weight[locus.organism])

        global_profile2weight = t.map_global_cdd_profile_count()

        if local_features:
            feature_labels = [ k for k,v in local_profile2weight.items() if v/weight >= 0.5 ]

        params = {}

        params['xls_file_name']         = xls_file_name
        params['loci']                  = cl_loci
        params['weight']                = weight
        params['profile_code2def']      = profile_code2def
        params['gnm2weight']            = gnm2weight
        params['feature_labels']        = feature_labels
        params['file2crispr_type']      = file2crispr_type
        params['local_profile2weight']  = local_profile2weight
        params['global_profile2weight'] = global_profile2weight

        r.write_to_xls_generic_loci(params)

    worksheet.write_row(ind+3, 0, ['Average entropy'], header_format)
    worksheet.write_row(ind+3, 1, [np.sum(weights*entropies)/np.sum(weights)])

    worksheet.write_row(ind + 4, 0, ['Exp(Average entropy)'], header_format)
    worksheet.write_row(ind + 4, 1, [np.exp(np.sum(weights * entropies) / np.sum(weights))])


def generate_cluster_reports_cas4(cluster_packs,
                                  loci,
                                  reports_dir,
                                  feature_labels,
                                  method=None,
                                  thresholds_pack=None):

    if method and thresholds_pack:

        thr_occ, thr_crisp, cluster_threshold = thresholds_pack

        summary_file = os.path.join(reports_dir,
                                    'summary_%s_%d_%.2f_%.2f.xlsx' % (method, thr_occ, thr_crisp, cluster_threshold))
    else:
        summary_file = os.path.join(reports_dir, 'summary.xlsx')

    workbook = x.Workbook(summary_file)
    worksheet = workbook.add_worksheet()

    header_format = workbook.add_format()
    header_format.set_font_size(12)
    header_format.set_bold()
    header_format.set_align('center')

    worksheet.set_column(2,3,20)
    worksheet.set_column(4,4,50)
    worksheet.set_column(5,5,50)

    worksheet.write_row(0, 0, ["File name", "Loci", "Effective size", "Entropy", "System types", "Genes", "Clusters"], header_format)

    print "Generating report files"
    ind = 0

    entropies = np.zeros(len(cluster_packs))

    for outer_i in range(len(cluster_packs)):

        (cluster, type2count, entropy, _) = cluster_packs[outer_i]

        ind += 1

        entropies[outer_i] = entropy

        crispr_cas_types_count = " ; ".join([k+":"+str(v) for (k,v) in sorted(type2count.items(), key=itemgetter(1), reverse=True)])

        xls_file_name = os.path.join(reports_dir, '%d.xlsx' % ind)
        tab_file_name = os.path.join(reports_dir, '%d.tab' % ind)
        dendrogram_file_name = os.path.join(reports_dir, '%d.png' % ind)

        cl_loci = [loci[_i] for _i in cluster]

        M = scores.jackard_weighted_scores(cl_loci)

        threshold = dendrogram.plot_dendrogram_from_score_matrix(M, dendrogram_file_name)

        # M gets modified in the previous procedure. Didn't have time to debug. 
        M = scores.jackard_weighted_scores(cl_loci)

        sub_singles, sub_clusters, gene2count = dendrogram.sub_classify_by_scores_cas4(M, threshold, cl_loci)
        
        params = {}

        params['xls_file_name']         = xls_file_name
        params['tab_file_name']         = tab_file_name
        params['loci']                  = cl_loci
        params['clusters']              = sub_clusters
        params['singles']               = sub_singles
        params['profile_code2def']      = profile_code2def
        params['feature_labels']        = feature_labels

        r.write_to_xls_generic_loci_cas4_vertical(params)
        # r.write_to_xls_generic_loci_cas4(params)
        r.write_loci_to_tab_file_cas4(params)

        sorted_gene2count = sorted(gene2count.items(), key=lambda x: x[1], reverse=True)

        gene_counts = ";". join( [ "%s:%.2f"% (gene_name, count) for (gene_name, count) in sorted_gene2count[:5]] )

        worksheet.write_row(ind, 0, ['%d.xlsx' % ind,
                                     len(cluster),
                                     len(sub_clusters) + len(sub_singles),
                                     entropy,
                                     crispr_cas_types_count,
                                     gene_counts,
                                     len(sub_clusters)])

    worksheet.write_row(ind+3, 0, ['Average entropy'], header_format)
    worksheet.write_row(ind+3, 1, [np.average(entropies)])

    worksheet.write_row(ind + 4, 0, ['Exp(Average entropy)'], header_format)
    worksheet.write_row(ind + 4, 1, [np.exp(np.average(entropies))])


def generate_cas4_gi_summary_file(singles,
                                  cluster_packs,
                                  loci,
                                  reports_dir,
                                  cluster2summary_file_name):

        cluster2summary = {int(l.split('\t')[0]):l.split('\t')[1].strip() for l in open(cluster2summary_file_name)}
        summary_file = open(os.path.join(reports_dir, 'cas4_gi_summary.tab'), 'w')

        cas_filter = set(["cas3", "cas5", "cas8c", "cas7", "cas4", "cas1", "cas2"])

        gi2crispr_type = {}

        ind = 1
        for outer_i in range(len(cluster_packs)):
            (cluster, type2count, entropy, gene2count) = cluster_packs[outer_i]

            sorted_gene2count = sorted(gene2count.items(), key=lambda x: x[1], reverse=True)
            top_genes = set(k for k,v in sorted_gene2count[:5])

            cl_loci = [loci[_i] for _i in cluster]

            cl_summary = cluster2summary[ind]

            for locus in cl_loci:

                cas4_genes = [g for g in locus.genes if g.is_seed]

                cl_genes = set(gene_name for gene in locus.genes for gene_name in gene.gene_name.split(','))

                # cas_genes = ",".join(cas_filter.intersection(cl_genes))
                cas_genes_summary = []

                buffer = []
                for gene in locus.genes:
                    for gene_name in set(gene.gene_name.split(",")):

                        if gene_name in top_genes:

                            buffer.append(gene_name)
                            cas_genes_summary += buffer
                            buffer = []
                            continue

                        buffer.append(gene_name if gene_name else "?")

                cas_genes_summary = "+".join(cas_genes_summary)

                for gene in cas4_genes:
                    summary_file.write("%s\t%s\t%s\t%s\t%s\n" % (gene.gid,
                                                             locus.organism,
                                                             cl_summary,
                                                             cas_genes_summary,
                                                             "%d.xlsx" % ind))

                    gi2crispr_type[gene.gid] = locus.crispr_type

            ind += 1

        for single_ind in singles:

            locus = loci[single_ind]

            cas4_genes = [g for g in locus.genes if g.is_seed]
            cas_genes_summary = []

            buffer = []
            for gene in locus.genes:
                for gene_name in gene.gene_name.split(","):

                    if gene_name in cas_filter:
                        buffer.append(gene_name)
                        cas_genes_summary += buffer
                        buffer = []

                    buffer.append(gene_name if gene_name else "?")

            cas_genes_summary = "+".join(cas_genes_summary)

            for gene in cas4_genes:
                summary_file.write("%s\t%s\t%s\t%s\n" % (gene.gid, locus.organism, "singleton", cas_genes_summary))
                gi2crispr_type[gene.gid] = locus.crispr_type

        summary_file.close()

        fname = os.path.join(gv.project_data_path, 'cas4/pickle/cas4_gi2crispr_type.p.bz2')

        t.dump_compressed_pickle(fname, gi2crispr_type)

        return gi2crispr_type



def generate_jw_cluster_reports(cluster_packs, loci, reports_dir, threshold):

    # if not feature_labels:
    #     local_features = True
    # else:
    #     local_features = False

    # thr_occ, thr_crisp, cluster_threshold = thresholds_pack

    summary_file = os.path.join(reports_dir,
                                'summary_jw_%.2f.xls' % threshold)

    workbook = x.Workbook(summary_file)
    worksheet = workbook.add_worksheet()

    header_format = workbook.add_format()
    header_format.set_font_size(12)
    header_format.set_bold()
    header_format.set_align('center')
    worksheet.set_column(4,5,50)
    worksheet.write_row(0, 0, ["File name", "Weight", "Loci", "Entropy", "systems weight", "systems count"], header_format)

    print "Generating report files"
    ind = 0

    weights = np.zeros(len(cluster_packs))
    entropies = np.zeros(len(cluster_packs))

    for outer_i in range(len(cluster_packs)):

        (cluster, type2count, type2weight, entropy) = cluster_packs[outer_i]

        ind += 1
        cl_files = [os.path.basename(loci[i].file_name) for i in cluster]

        weight = sum([gnm2weight[file2org[file]] for file in cl_files])

        weights[outer_i] = weight
        entropies[outer_i] = entropy

        crispr_cas_types_count = " ; ".join([k+":"+str(v) for (k,v) in sorted(type2count.items(), key=itemgetter(1), reverse=True)])
        crispr_cas_types_weight = " ; ".join([k+":"+str(v) for (k,v) in sorted(type2weight.items(), key=itemgetter(1), reverse=True)])

        xls_file_name = os.path.join(reports_dir, '%d.xls' % ind)

        worksheet.write_row(ind+1, 0, ['%d.xls'%ind,
                                       weight,
                                       len(cl_files),
                                       entropy,
                                       crispr_cas_types_weight,
                                       crispr_cas_types_count,
                                       " "])

        cl_loci = sorted([loci[_i] for _i in cluster], key = lambda x: gnm2weight[x.organism], reverse=True)

        local_profile2weight = {}
        for locus in cl_loci:
            for gene in locus.genes:
                for profile in gene.cogid.split(','):
                    t.update_dictionary(local_profile2weight, profile, gnm2weight[locus.organism])

        global_profile2weight = t.map_global_cdd_profile_count()

        # if local_features:
        #     feature_labels = [ k for k,v in local_profile2weight.items() if v/weight >= 0.5 ]

        params = {}

        params['xls_file_name']         = xls_file_name
        params['loci']                  = cl_loci
        params['weight']                = weight
        params['profile_code2def']      = profile_code2def
        params['gnm2weight']            = gnm2weight
        # params['feature_labels']        = feature_labels
        params['feature_labels']        = []
        params['file2crispr_type']      = file2crispr_type
        params['local_profile2weight']  = local_profile2weight
        params['global_profile2weight'] = global_profile2weight

        r.write_to_xls_generic_loci(params)

    worksheet.write_row(ind+3, 0, ['Average entropy'], header_format)
    worksheet.write_row(ind+3, 1, [np.sum(weights*entropies)/np.sum(weights)])

    worksheet.write_row(ind + 4, 0, ['Exp(Average entropy)'], header_format)
    worksheet.write_row(ind + 4, 1, [np.exp(np.sum(weights * entropies) / np.sum(weights))])



def generate_community_reports(nodes_pool,
                               reports_dir,
                               locus2weight,
                               file2locus,
                               profile2def,
                               feature_profiles_file=None):

    # if not feature_labels:
    #     local_features = True
    # else:
    #     local_features = False

    # thr_occ, thr_crisp, cluster_threshold = thresholds_pack

    summary_file = os.path.join(reports_dir, 'summary.xlsx' )

    workbook = x.Workbook(summary_file)
    worksheet = workbook.add_worksheet()

    header_format = workbook.add_format()
    header_format.set_font_size(12)
    header_format.set_bold()
    header_format.set_align('center')
    # worksheet.set_column(4,5,50)
    worksheet.write_row(0, 0, ["File name", "Size", "Effective size", "Genes"], header_format)

    print "Generating report files"
    ind = 1

    for nodes in nodes_pool:

        loci_size = len([node for node in nodes if node.type == 2])
        loci_esize = sum(node.weight for node in nodes if node.type == 2)

        # if loci_esize < 5:
        #     continue

        loci = [file2locus[node.file_name] for node in nodes if node.type == 2]

        xls_file_name = os.path.join(reports_dir, '%d.xlsx' % ind)
        loci_file_name = os.path.join(reports_dir, '%d.tab' % ind)

        with open(loci_file_name, 'w') as outf:
            loci_files = ",".join(os.path.basename(locus.file_name) for locus in loci)
            outf.write(loci_files + "\n")

        gene2cnt = {}
        profile2cnt = {}
        for locus in loci:
            weight = locus2weight[os.path.basename(locus.file_name)]
            for gene_name in locus.gene_names:
                t.update_dictionary(gene2cnt, gene_name, weight)
            for cl in locus.clusters:
                t.update_dictionary(profile2cnt, cl, weight)

        sorted_gene2count = sorted(gene2cnt.items(), key=lambda x: x[1], reverse=True)
        gene_counts = ";".join(["%s:%.2f" % (gene_name, count) for (gene_name, count) in sorted_gene2count[:10]])

        worksheet.write_row(ind + 1, 0, ['%d.xlsx' % ind,
                                         loci_size,
                                         loci_esize,
                                         gene_counts])
        args = {}

        args['xls_file_name'] = xls_file_name
        args['loci'] = loci
        args['profile_code2def'] = profile2def
        if not feature_profiles_file:
            args['feature_labels'] = [ k for k,v in profile2cnt.items() if v >= loci_esize/2 ]
        else:
            args['feature_labels'] = [l.strip() for l in open(feature_profiles_file)]

        try:
            r.write_to_xls_loci_plain(args)
        except:
            sys.exit()
        ind += 1
