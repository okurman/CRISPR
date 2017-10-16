#!/usr/bin/env python
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

from lib.utils.clustering import Locus
import lib.utils.tools as t


def get_feature_labels():

    feature_definition_file = os.path.join(gv.project_data_path, 'cas4/profiles/profiles_10_0.50.tab')

    feature_labels = []

    cluster_description_file = os.path.join(gv.project_data_path, 'cas4/baiticity/cluster_id_size.tsv')
    cluster_id2effective_size = {l.split()[0]: l.strip().split()[1] for l in open(cluster_description_file)}

    for id in open(feature_definition_file).readlines():
        id = id.strip()
        feature_labels.append("CLUSTER_%s_%s" % (id, cluster_id2effective_size[id]))

    return set(feature_labels)


if __name__=='__main__':

    files_path = os.path.join(gv.project_data_path, 'cas4/files/')

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

    loci = [Locus(os.path.join(files_path, f), file_format='generic', profile2gene=cdd_profile2gene)
            for f in os.listdir(files_path)]

    locus2id = {}
    cluster2id = {}
    feature_clusters = get_feature_labels()

    print "No of feature clusters", len(feature_clusters)

    vertices, edges = [], []

    id = 1
    for locus in loci:
        locus2id[locus.base_file_name] = id

        _crispr_type = locus.crispr_type
        _type_code = None

        if _crispr_type.startswith("part unk") or _crispr_type.startswith("unk"):

            _type_code = 0

        elif _crispr_type.startswith("CAS-I-") or _crispr_type.strip()=="CAS-I" \
                or _crispr_type.strip()=="part CAS-I" \
                or _crispr_type.startswith("part CAS-I-"):

            _type_code = 1

        elif _crispr_type.startswith("CAS-II-"):

            _type_code = 2

        elif _crispr_type.startswith("CAS-III")\
                or _crispr_type.strip() == "part CAS-III" \
                or _crispr_type.startswith("part CAS-III-"):

            _type_code = 3

        elif _crispr_type.startswith("CAS-IV-") \
                or _crispr_type.startswith("part CAS-IV-"):

            _type_code = 4

        elif _crispr_type.startswith("CAS-V-"):

            _type_code = 5

        else:

            _type_code = 6
            print "Not parsed:", _crispr_type

        vertices.append("%d\t\"%s %s\"\t1\t%d\n" % (id, locus.crispr_type, locus.base_file_name, _type_code))
        id += 1

    for cluster in feature_clusters:
        vertices.append("%d\t\"%s\"\t2\t7\n" % (id, cluster))
        cluster2id[cluster] = id
        id += 1

    ########################################################################
    ################## This part is for preparing pajek file ###############
    ########################################################################
    # net_file = open(os.path.join(gv.project_data_path, 'cas4/infomap/cas4.net'), 'w')
    # for locus in loci:
    #     locus_id = locus2id[locus.base_file_name]
    #     for cluster in locus.clusters:
    #         if cluster in cluster2id:
    #             cluster_id = cluster2id[cluster]
    #             edges.append("%d %d 1\n" % (locus_id, cluster_id))
    #
    # net_file.write("*Vertices %d\n" % len(vertices))
    # [net_file.write(l) for l in vertices]
    #
    # net_file.write("*Edges %d\n"%len(edges))
    # [net_file.write(l) for l in edges]
    # net_file.close()

    ########################################################################
    ############ This part is for preparing link file for cytoscape ########
    ########################################################################
    for locus in loci:
        locus_id = locus2id[locus.base_file_name]

        for cluster in locus.clusters:
            if cluster in cluster2id:
                cluster_id = cluster2id[cluster]
                edges.append("%d %d\n" % (locus_id, cluster_id))

    # net_file.write("*Vertices %d\n" % len(vertices))
    # [net_file.write(l) for l in vertices]
    #
    # net_file.write("*Edges %d\n"%len(edges))
    # [net_file.write(l) for l in edges]
    # net_file.close()

    vertex_description_file = open(os.path.join(gv.project_data_path, 'cas4/infomap/vertices.txt'), 'w')
    vertex_description_file.write("#Id\tDescription\tType\tCRISPR type\n")
    [vertex_description_file.write(l) for l in vertices]
    vertex_description_file.close()

    net_file = open(os.path.join(gv.project_data_path, 'cas4/infomap/cas4_cytoscape.net'), 'w')
    [net_file.write(l) for l in edges]
    net_file.close()