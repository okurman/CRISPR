#!/home/hudaiber/env2.7/bin/python

import os
import sys
import configparser
from collections import defaultdict, Counter

###############################################################
config_file = os.path.join(os.path.expanduser('~'),'paths.cfg')
cfg=configparser.ConfigParser()
cfg.read(config_file)

code_path = cfg.get('NewSystems','code_path')
data_path = cfg.get('NewSystems','data_path')
sys.path.append(code_path)
###############################################################

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import powerlaw
import lib.utils.graph.tools as gt
import lib.utils.tools as t
from lib.utils import infomap

PROK1603_WEIGHT = 4930.0
RAND_ADJ_PROB = 0.0451314274706
work_dir = os.path.join(data_path, "prok1603/graph/")



def cas_synonyms():

    _f = os.path.join(work_dir, "infomap/merged/cas_synonyms.txt")
    cl2name = {}
    for l in open(_f):
        if l.startswith("#"):
            continue
        parts = l.split("\t")
        cl2name[parts[0]] = parts[1]

    return cl2name


def extract_cas_communities(tree_file, node_file, G=None):

    module2members = infomap.parse_tree_file(tree_file, node_names=True, node_desc_file=node_file)

    modules = sorted([(k,v) for (k,v) in module2members.items()], key=lambda x: len(x[1]), reverse=True)

    cl2name = cas_synonyms()
    profile2gene = t.map_cdd_defense2gene_name()

    def_profs_file = "/panfs/pan1/patternquest/Projects/NewSystems/data/profiles/defenseProfiles.tab"
    cas_profiles = [l.split()[0] for l in open(def_profs_file) if l.split()[1] == "CRISPR"]

    # uvrd_profiles = ["COG0210", "PRK11773", "TIGR01075", "cl22977", "cl26202", "cl28312", "pfam00580", "pfam13361",
    #                  "pfam13538"]

    cas_profiles = cas_profiles + list(cl2name.keys())

    save_dir = os.path.join(work_dir, "cases/crispr_graph")

    for module in modules:
        nodes = set(module[1])
        cas_nodes = nodes.intersection(cas_profiles)

        if cas_nodes:

            # if "crispr" in nodes:
            #     print "CRISPR"

            cas_ratio = len(cas_nodes)/float(len(nodes))
            
            if cas_ratio < 0.3:
                continue

            print module[0], len(nodes), cas_ratio

            continue

            graph_file = os.path.join(save_dir, "%s_%d.net" % (module[0], len(nodes)))
            node_file = os.path.join(save_dir, "%s_%d.nodes.txt" % (module[0], len(nodes)))

            subgraph = gt.subgraph_crosslink(G, list(nodes))
            gt.write_edgelist(subgraph, graph_file)

            with open(node_file, "w") as outf:
                outf.write("Node\tgene_name\t")
                for node in subgraph:
                    if node in cl2name:
                        _name = cl2name[node]
                    elif node in profile2gene:
                        _name = profile2gene[node]
                    else:
                        _name = node

                    outf.write("%s\t%s\n" % (node, _name))
            break


def extract_extract_cas_communities_from_oslom(G):

    oslom_file = "/panfs/pan1/patternquest/Projects/NewSystems/data/prok1603/graph/infomap/merged/consensus_clustering/multi_level/norm_min/3/oslom/tp"
    
    module2pvalue = {}
    module2members = {}
    with open(oslom_file) as inf:
        l = inf.readline()

        while l:
            assert l.startswith("#")

            parts = l.split()
            name = parts[1]
            size = int(parts[3])
            if size ==1:
                l = inf.readline()
                l = inf.readline()
                continue

            pvalue = float(parts[5])

            module2pvalue[name] = pvalue
            module2members[name] = inf.readline().strip().split()
            l = inf.readline()

    # modules = [l.strip().split() for l in open(oslom_file) if not l.startswith("#")]

    node_file = os.path.join(work_dir, "infomap/merged/nodes.txt")

    id2node = {l.split()[0]:l.split()[1] for l in open(node_file)}

    def_profs_file = "/panfs/pan1/patternquest/Projects/NewSystems/data/profiles/defenseProfiles.tab"
    cas_profiles = [l.split()[0] for l in open(def_profs_file) if l.split()[1] == "CRISPR"]
    cl2name = cas_synonyms()
    cas_profiles = cas_profiles + list(cl2name.keys())
    save_dir = os.path.join(work_dir, "cases/crispr_graph_2/")

    profile2gene = t.map_cdd_defense2gene_name()

    for (name, members) in module2members.items():
        # nodes = set(module[1])
        nodes = set([id2node[id] for id in members])
        cas_nodes = nodes.intersection(cas_profiles)

        if cas_nodes:

            cas_ratio = len(cas_nodes) / float(len(nodes))

            if cas_ratio < 0.3:
                continue

            if "crispr" in nodes:
                print "CRISPR"

            # print len(nodes), cas_ratio, " ".join(nodes)
            # print len(nodes), cas_ratio, Counter([profile2gene[n] for n in cas_nodes])
            vals = [name, str(module2pvalue[name]), str(len(nodes)), str(cas_ratio), " ".join(nodes)]
            # print "\t".join(vals)

            graph_file = os.path.join(save_dir, "%s_%d.net" % (name, len(nodes)))
            node_file_s = os.path.join(save_dir, "%s_%d.nodes.txt" % (name, len(nodes)))

            subgraph = gt.subgraph_crosslink(G, list(nodes))
            gt.write_edgelist(subgraph, graph_file)

            with open(node_file_s, "w") as outf:
                outf.write("Node\tgene_name\tcas\n")
                for node in subgraph:
                    if node in cl2name:
                        _name = cl2name[node]
                    elif node in profile2gene:
                        _name = profile2gene[node]
                    else:
                        _name = node

                    _type = 0
                    if node in cas_nodes:
                        _type = 1

                    outf.write("%s\t%s\t%d\n" % (node, _name, _type))




if __name__=="__main__":

    inf_file = os.path.join(work_dir, "infomap/merged/out/4/link_list.tree")
    node_file = os.path.join(work_dir, "infomap/merged/nodes.txt")

    # communities_in_merged_graph(inf_file, node_file)
    G = nx.read_gpickle(os.path.join(work_dir, "adj_graph_merged_prp.p"))
    extract_extract_cas_communities_from_oslom(G)
