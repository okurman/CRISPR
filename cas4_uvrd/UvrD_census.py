#!/home/hudaiber/env2.7/bin/python
__author__ = 'hudaiber'

import os
import sys
if sys.platform == 'darwin':
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))
elif sys.platform == 'linux2':
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))

import global_variables as gv
sys.path.append(gv.project_code_path)

from lib.utils import BasicLocus, Locus
import lib.utils.clustering.dendrogram as dnd
import lib.utils.clustering.reporting as r
from lib.utils import tools as t
from lib.utils.clustering import scores
import numpy as np
import cPickle
import time
import subprocess as sp
import glob
import dm_tools as dt
import random
from lib.utils.prok1603 import parser


system_gene_names = set(['uvra','uvrb','uvrc','uvrd','sbcc','sbcd', 'cas4'])

def cas4_dataset():

    print "Loading loci"
    def_file = os.path.join(gv.project_data_path, 'cas4/profiles/defenseProfiles.tab')
    profile2gene = {}

    for l in open(def_file):
        terms = l.split('\t')
        profile = terms[0]
        gene_names = terms[3].split(',')
        if len(gene_names) > 1:
            profile2gene[profile] = gene_names[1]
        else:
            profile2gene[profile] = gene_names[0]

    cdd_profile2gene = t.map_cdd_profile2gene_name()
    cdd_profile2gene.update(profile2gene)

    loci = [Locus(os.path.join(work_dir, 'files', f), file_format='generic', profile2gene=cdd_profile2gene)
            for f in os.listdir(os.path.join(work_dir, 'files'))]

    # dendrogram_file = os.path.join(work_dir, 'loci_dendrogram.pdf')

    singles, cluster_packs, _ = dnd.classify_loci_hierarchically(loci, threshold=2.6)

    print "Clusters: %d, singles %d" % (len(cluster_packs), len(singles))

    reports_dir = os.path.join(work_dir, 'reports/cas4')

    feature_profiles = [l.strip() for l in open(os.path.join(work_dir, 'feature_profiles.txt')).readlines()]

    r.generate_cluster_reports_cas4(cluster_packs, loci, reports_dir, feature_profiles)


def prok1603_extract_dendrogram():

    work_dir = os.path.join(gv.project_data_path, 'UvrD/')

    files_path = os.path.join(work_dir, 'prok1603/merged_files/')

    print "Loading loci"

    def_file = os.path.join(gv.project_data_path, 'cas4/profiles/defenseProfiles.tab')
    profile2gene = {}

    for l in open(def_file):
        terms = l.split('\t')
        profile = terms[0]
        gene_names = terms[3].split(',')
        if len(gene_names) > 1:
            profile2gene[profile] = gene_names[1]
        else:
            profile2gene[profile] = gene_names[0]

    cdd_profile2gene = t.map_cdd_profile2gene_name()
    cdd_profile2gene.update(profile2gene)

    loci = [BasicLocus(os.path.join(files_path, f), profile2gene=cdd_profile2gene) for f in
            os.listdir(files_path)]

    prok1603_loci_file = os.path.join(work_dir, 'prok1603/prok1603_loci.p.bz2')
    # loci = t.load_compressed_pickle(prok1603_loci_file)
    print "Loci:", len(loci)
    print "Dumping loci to:", prok1603_loci_file

    t.dump_compressed_pickle(prok1603_loci_file, loci)
    sys.exit()

    tic = time.time()
    print "Generating score matrix"
    M = scores.generate_jackard_score_matrix(loci)

    tac = time.time() - tic
    print "Elapsed time:", float(tac)/60/60, float(tac)/60, float(tac)

    tic = time.time()
    prok1603_jw_file = os.path.join(work_dir, 'prok1603_jw_scores.p')
    print "Dumping JW scores to:", prok1603_jw_file
    with open(prok1603_jw_file, 'wb') as outf:
        cPickle.dump(M, outf, protocol=cPickle.HIGHEST_PROTOCOL)
    tac = time.time() - tic
    print "Elapsed time:", float(tac) / 60 / 60, float(tac) / 60, float(tac)

    tic = time.time()
    prok1603_jw_file = os.path.join(work_dir, 'prok1603_jw_scores.p.bz2')
    print "Dumping JW scores to:", prok1603_jw_file
    t.dump_compressed_pickle(prok1603_jw_file, M)
    tac = time.time() - tic
    print "Elapsed time:", float(tac) / 60 / 60, float(tac) / 60, float(tac)

    # print "Loading JW scores from:", prok1603_jw_file
    # M = t.load_compressed_pickle(prok1603_jw_file)

    tic = time.time()
    prok1603_jw_file = os.path.join(work_dir, 'prok1603/prok1603_jw_scores.npz')
    # prok1603_jw_file = os.path.join('/Users/hudaiber/Projects/NewSystems/data/UvrD/prok1603/prok1603_jw_scores.npz')
    print "Dumping JW scores to:", prok1603_jw_file
    np.savez_compressed(prok1603_jw_file, M)
    tac = time.time() - tic
    print "Elapsed time:", float(tac) / 60 / 60, float(tac) / 60, float(tac)

    # print "Loading JW scores from:", prok1603_jw_file
    # M = t.load_compressed_pickle(prok1603_jw_file)

    prok1603_tree_file = os.path.join(work_dir, 'prok1603/prok1603_upgma.tre')
    print "Generating tree:", prok1603_tree_file
    dnd.convert_score_to_newick(M, [os.path.basename(l.file_name) for l in loci], prok1603_tree_file)


def cas4_extract_dendrogram():

    work_dir = os.path.join(gv.project_data_path, 'cas4')

    print "Loading loci"

    def_file = os.path.join(gv.project_data_path, 'cas4/profiles/defenseProfiles.tab')
    profile2gene = {}

    for l in open(def_file):
        terms = l.split('\t')
        profile = terms[0]
        gene_names = terms[3].split(',')
        if len(gene_names) > 1:
            profile2gene[profile] = gene_names[1]
        else:
            profile2gene[profile] = gene_names[0]

    cdd_profile2gene = t.map_cdd_profile2gene_name()
    cdd_profile2gene.update(profile2gene)

    files_path = os.path.join(work_dir, 'files')

    loci = [Locus(os.path.join(files_path, f), file_format='generic', profile2gene=cdd_profile2gene) for f in
            os.listdir(files_path)]

    tic = time.time()
    print "Generating score matrix"
    M = scores.generate_jackard_score_matrix(loci)

    tac = time.time() - tic
    print "Elapsed time:", float(tac) / 60 / 60, float(tac) / 60, float(tac)

    tic = time.time()
    jw_file = os.path.join(work_dir, 'pickle/jw_scores.p.bz2')
    print "Dumping JW scores to:", jw_file
    t.dump_compressed_pickle(jw_file, M)
    tac = time.time() - tic
    print "Elapsed time:", float(tac) / 60 / 60, float(tac) / 60, float(tac)

    # print "Loading JW scores from:", prok1603_jw_file
    # M = t.load_compressed_pickle(prok1603_jw_file)

    tree_file = os.path.join(work_dir, 'jw_upgma.tre')
    print "Generating tree:", tree_file
    dnd.convert_score_to_newick(M, [os.path.basename(l.file_name) for l in loci], tree_file)


def prok1603_architecture_frequencies():

    work_dir = os.path.join(gv.project_data_path, 'UvrD/')

    map_file = os.path.join(work_dir, 'prok1603/prok1603_weights.txt')

    locus2weight = {l.split()[0]:float(l.split()[1]) for l in open(map_file)}

    def_file = os.path.join(gv.project_data_path, 'cas4/profiles/defenseProfiles.tab')
    profile2gene={}
    profile2def = {}

    for l in open(def_file):
        terms = l.strip().split('\t')
        profile = terms[0]
        gene_names = terms[3].split(',')
        if len(gene_names)>1:
            profile2gene[profile] = gene_names[1]
        else:
            profile2gene[profile] = gene_names[0]

        profile2def[profile] = terms[4]

    cdd_profile2gene = t.map_cdd_profile2gene_name()
    cdd_profile2gene.update(profile2gene)

    cdd_profile2def = t.map_cdd_profile2def()
    cdd_profile2def.update(profile2def)

    prok1603_loci_file = os.path.join(work_dir, 'prok1603The CRISPR/prok1603_loci.p.bz2')
    loci = t.load_compressed_pickle(prok1603_loci_file)

    profile2loci = {}

    for locus in loci:
        for _profile in locus.profiles:
            if _profile in profile2loci:
                profile2loci[_profile].append(locus)
            else:
                profile2loci[_profile] = [locus]

    for (profile, loci) in sorted(profile2loci.items(), key=lambda x: len(x[1]), reverse=True):
        _weight = sum([locus2weight[locus.base_file_name] for locus in loci])
        print "%s\t%s\t%d\t%f\t%s" % (profile,
                                      cdd_profile2gene[profile] if profile in cdd_profile2gene else "",
                                      len(loci),
                                      _weight,
                                      cdd_profile2def[profile] if profile in cdd_profile2def else "")


def get_prok1603_weights():

    work_dir = os.path.join(gv.project_data_path, 'UvrD/prok1603')
    files_dir = os.path.join(work_dir, 'merged_files')




def get_system_genes(gi, files_dir, profile2gene):

    pty_file = os.path.join(files_dir, "%s.pty" % gi)

    if not os.path.exists(pty_file):
        return None

    _genes = dt.get_pty_file(pty_file, profile2gene=profile2gene)

    res_genes = set([])

    for _gene in _genes:

        for _gene_name in _gene.gene_name.split(','):

            if not _gene_name.strip():
                continue

            for system_gene in system_gene_names:

                if system_gene in _gene_name.lower():
                    res_genes.update([_gene_name])
                    break

    return "|".join(res_genes)


def tree_leaves():

    work_dir = os.path.join(gv.project_data_path, 'UvrD/prok1603')
    tree_dir = os.path.join(work_dir, 'clust_tree/')
    files_dir = os.path.join(work_dir, 'merged_files')

    profile2gene = t.map_cdd_profile2gene_name()

    gi2org = {l.split()[0]: l.rstrip().split()[1] for l in open(work_dir + '/gi_org.txt')}

    gi2weight = {l.split()[0].split('.')[0]: float(l.split()[1]) for l in open(work_dir + '/prok1603_weights.txt')}

    cl2size, cl2gis, cl2weight = {}, {}, {}

    for l in open(tree_dir + 'uvrd.cls'):
        terms = l.rstrip().split()
        cl2size[terms[1]] = terms[0]
        cl2gis[terms[1]] = terms[2:]
        cl2weight[terms[1]] = sum([ gi2weight[gi] if gi in gi2weight else 0 for gi in terms[2:]])

    tree_string = open(tree_dir + 'uvrd.up.tre').readline()

    leave_file_names = [os.path.basename(l) for l in glob.glob(tree_dir + '*.sr')]

    for leave_file_name in leave_file_names:

        leave_file_gis = [ l.split()[0] for l in open(os.path.join(tree_dir, leave_file_name))]

        system_gene_pool = []

        sgp_count = {}

        for gi in leave_file_gis:
            system_genes = get_system_genes(gi, files_dir, profile2gene)

            if not system_genes:
                continue
            system_gene_pool.append(system_genes)

            t.update_dictionary(sgp_count, system_genes, gi2weight[gi])

        sorted_sgp_count = sorted(sgp_count.items(), key=lambda x: x[1], reverse=True)

        leaf_name = os.path.splitext(leave_file_name)[0]
        gene_names = sorted_sgp_count[0][0] if sorted_sgp_count else ""

        representative = gi2org[random.choice(leave_file_gis)]

        total_weight = sum([v for k,v in sgp_count.items()])
        leaf_prefix = "%s|" % int(total_weight) if total_weight else "-"

        has_genes = False

        for _gene_name in ["Cas4", "UvrA", "UvrB", "UvrC", "SbcS", "SbcD"]:
        # for _gene_name in ["Cas4"]:
            _weight = sum([v for k, v in sgp_count.items() if _gene_name.lower() in k.lower()])
            if _weight:
                leaf_prefix += "%s=%d|" % (_gene_name, _weight)
                has_genes = True

        if has_genes:
            new_leaf_name = leaf_prefix + representative + "|" + leaf_name.split('.')[1]
        else:
            new_leaf_name = leaf_prefix + leaf_name.split('.')[1]

        # new_leaf_name = "cas4=%s/%s|%s|%s" % (int(cas4_weight) if total_weight else "-",
        #                                  int(total_weight) if total_weight else "-",
        #                                  representative,
        #                                  leaf_name.split('.')[1])

        print leaf_name, new_leaf_name

        tree_string = tree_string.replace(leaf_name + ":", new_leaf_name + ":")

        # new_file_name = os.path.join(tree_dir, os.path.splitext(leave_file_name)[0] + '.def')
        # with open(new_file_name, 'w') as new_file:
        #
        #     for k, v in sorted_sgp_count:
        #         new_file.write("#%s\t%f\n" % (k, v))
        #
        #     new_file.write("\n")
        #
        #     [new_file.write("%s\t%s\n" % (gi, gi2org[gi])) for gi in leave_file_gis]

    with open(tree_dir + 'uvrd.up_all_genes.tree', 'w') as outf:
        outf.write(tree_string)


def prok1603_cas4_crispr_uvrd():

    loci_file = '/panfs/pan1.be-md.ncbi.nlm.nih.gov/prokdata/CRISPRicity/All1603_pty_loci/Islands.ann_full'
    work_dir = os.path.join(gv.project_data_path, 'UvrD/prok1603')

    prok1603_genomes = os.listdir(os.path.join(os.path.expanduser('~'), 'data/Prok1603/pty/'))
    org2weight = parser.map_org2weight()

    all_crispr_cas4_gis = []
    all_fusion_gis = []
    all_cas4_uvrd_gis = []
    all_solo_cas4_gis = []

    loci = t.parse_crispr_loci(loci_file, filter_set=prok1603_genomes)

    org2loci = {}

    for locus in loci:

        organism = locus.organism

        if organism not in org2loci:
            org2loci[organism] = [locus]
        else:
            org2loci[organism].append(locus)

    header = ["Organism", "Weight", "Fusion", "cas4-UvrD", "Total", "Solo Cas4", "Solo UvrD", "Total", "Known", "Unknown", "Cas4+", "Cas4-"]

    subtypes = ['CAS-I',
                'CAS-I-A',
                'CAS-I-B',
                'CAS-I-C',
                'CAS-I-D',
                'CAS-I-E',
                'CAS-I-F',
                'CAS-I-U',
                'CAS-II-A',
                'CAS-II-B',
                'CAS-II-C',
                'CAS-III',
                'CAS-III-A',
                'CAS-III-B',
                'CAS-III-C',
                'CAS-III-D',
                'CAS-IV-A',
                'CAS-V-A',
                'CAS-V-B',
                'CAS-V-D',
                'CAS-V-E',
                'CAS-VI-A',
                'CAS-VI-B',
                'Unknown']

    outf = open(os.path.join(work_dir, 'census_cas4_uvrd.tsv'), 'w')
    outf.write("\t".join(header + subtypes) + "\n")

    for cnt, _org in enumerate(prok1603_genomes):

        subtypes_count = [0] * len(subtypes)

        if _org in org2loci:
            w_cas4 = len([l for l in org2loci[_org] if len(l.cas4_genes)>0])
            wo_cas4 = len([l for l in org2loci[_org] if len(l.cas4_genes)==0])
            crispr_cas4_gis = set([gene.gid for locus in org2loci[_org] for gene in locus.cas4_genes])
            crispr_gis = set([gid for locus in org2loci[_org] for gid in locus.gis])

            known_loci = len([locus for locus in org2loci[_org] if locus.status == "Known"])
            unknown_loci = len([locus for locus in org2loci[_org] if locus.status != "Known"])

            for locus in org2loci[_org]:

                if "unk" in locus.type or "Unidentified" in locus.type:
                    _type="Unknown"
                    subtypes_count[subtypes.index(_type)] += 1
                    continue
                    
                for _type in locus.type.split(";"):
                    if "part" in _type:
                        _type = _type.split()[1]
                    subtypes_count[subtypes.index(_type)] += 1
        else:

            w_cas4, wo_cas4 = 0, 0
            crispr_cas4_gis = set([])
            crispr_gis = set([])

            known_loci, unknown_loci = 0, 0

        solo_cas4 = set()
        solo_uvrd = set()

        joint_cas4 = set()
        joint_uvrd = set()

        crispr_uvrd = set()

        fusion = set([])

        for source_genes in parser.get_pty_organism(_org):

            for i in range(len(source_genes)):
                _gene = source_genes[i]
                _gene_name = _gene.gene_name.lower()

                if "cas4" in _gene_name:

                    if _gene.gid in crispr_cas4_gis:
                        continue

                    if "uvrd" in _gene_name:
                        fusion.update([_gene.gid])
                        continue

                    for _sub_gene in source_genes[i-5:i+5]:
                        if "uvrd" in _sub_gene.gene_name.lower():
                            joint_cas4.update([_gene.gid])
                            print _org, _gene.gid
                            continue

                    solo_cas4.update([_gene.gid])

                if "uvrd" in _gene_name:

                    if _gene.gid in crispr_gis:
                        crispr_uvrd.update([_gene.gid])
                        continue

                    for _sub_gene in source_genes[i-5:i+5]:

                        if _sub_gene.gid in fusion or _sub_gene.gid in crispr_gis:
                            continue

                        if "cas4" in _sub_gene.gene_name.lower():

                            # if the UvrD comes before the fusion, skip it, let the cas4 section above catch it as a fusion
                            if "uvrd" in _sub_gene.gene_name.lower():
                                continue

                            joint_uvrd.update([_gene.gid])
                            continue

                    solo_uvrd.update([_gene.gid])

        out_terms = [_org, org2weight[_org], len(fusion), len(joint_cas4), len(fusion)+len(joint_cas4), len(solo_cas4),
                     len(solo_uvrd), known_loci + unknown_loci, known_loci, unknown_loci, w_cas4, wo_cas4]

        out_terms += subtypes_count

        for i in range(2,len(out_terms)):
            out_terms[i] *= float(org2weight[_org]) if org2weight[_org] else 1

        outf.write("\t".join([str(_) for _ in out_terms]) + "\n")

        all_cas4_uvrd_gis += joint_cas4
        all_fusion_gis += fusion
        all_crispr_cas4_gis += crispr_cas4_gis
        all_solo_cas4_gis += solo_cas4

    outf.close()

    for i in range(4):
        _fnames= ['cas4_uvrd_gis.txt', 'fusion_gis.txt', 'crispr_cas4_gis.txt', 'solo_cas4_gis.txt']
        _gi_lists = [all_cas4_uvrd_gis, all_fusion_gis, all_crispr_cas4_gis, all_solo_cas4_gis]

        with open(os.path.join(work_dir, 'dnds', _fnames[i]), 'w') as outf:
            [outf.write(gi + "\n") for gi in _gi_lists[i]]




if __name__=='__main__':

    # prok1603_extract_dendrogram()
    # cas4_extract_dendrogram()

    # prok1603_architecture_frequencies()

    # tree_leaves()

    # prok1603_cas4_crispr_uvrd()


    loci_file = '/panfs/pan1.be-md.ncbi.nlm.nih.gov/prokdata/CRISPRicity/All1603_pty_loci/Islands.ann_full'

    prok1603_genomes = os.listdir(os.path.join(os.path.expanduser('~'), 'data/Prok1603/pty/'))

    loci = t.parse_crispr_loci(loci_file, filter_set=prok1603_genomes)

    print "Complete loci:\n%s\n"," ".join([gene.gid for locus in loci if locus.complete for gene in locus.cas4_genes])
    print "InComplete loci:\n%s\n", " ".join([gene.gid for locus in loci if not locus.complete for gene in locus.cas4_genes])
