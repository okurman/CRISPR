#!/home/hudaiber/env2.7/bin/python

import os
import sys
import configparser
import networkx as nx

###############################################################
config_file = os.path.join(os.path.expanduser('~'),'paths.cfg')
cfg=configparser.ConfigParser()
cfg.read(config_file)

code_path = cfg.get('NewSystems','code_path')
data_path = cfg.get('NewSystems','data_path')
sys.path.append(code_path)
###############################################################

import lib.utils.tools as t
from lib.utils.classes import Kplet
import lib.db.prok1402.db_tools as dt
from lib.db.prok1402.duplets import extract_baited_duplet_aggregates
import numpy as np
import matplotlib.pyplot as plt

file_name2file_id, file_id2file_name = dt.map_baited_file2id()
genome2weight = t.map_genome2weight()


def load_adjacent_duplets_graph():

    duplet_rows = extract_baited_duplet_aggregates()

    duplets = []
    graph = nx.DiGraph()

    for row in duplet_rows:
        kplet_id = row[0]
        profiles = (row[1], row[2])
        number_of_loci = int(row[3])
        weight_of_loci = float(row[4])
        locus_file_ids = row[5].split(",")

        assert(len(locus_file_ids) == number_of_loci)

        graph.add_edge(row[1], row[2], weight=weight_of_loci)

        _kplet = Kplet(id=kplet_id,
                       profiles=profiles,
                       count=number_of_loci,
                       weight=weight_of_loci,
                       files=locus_file_ids)

        duplets.append(_kplet)

    return graph, duplets


def load_duplets_graph():

    duplet_rows = extract_baited_duplet_aggregates("all")

    duplets = []
    graph = nx.DiGraph()

    for row in duplet_rows:
        kplet_id = row[0]
        profiles = (row[1], row[2])
        number_of_loci = int(row[3])
        weight_of_loci = float(row[4])
        locus_file_ids = row[5].split(",")

        assert(len(locus_file_ids) == number_of_loci)
        graph.add_edge(row[1], row[2], weight=weight_of_loci)
        _kplet = Kplet(id=kplet_id,
                       profiles=profiles,
                       count=number_of_loci,
                       weight=weight_of_loci,
                       files=locus_file_ids)

        duplets.append(_kplet)

    return graph, duplets


def generate_graphs_from_baited_loci():
    
    # Adjacent duplets
    print "Loading adjacent duplets"
    graph, duplets = load_adjacent_duplets_graph()
    linked_list_file = os.path.join(work_dir, "baited_adj_duplets.txt")
    with open(linked_list_file, "w") as outf:
        outf.write("#From\tTo\tWeight\n")
        for edge in graph.edges_iter(data=True):
            _from = edge[0]
            _to = edge[1]
            _weight = edge[2]["weight"]
            outf.write("%s\t%s\t%s\n" % (_from, _to, _weight))
    print "Writing adjacent duplets"
    pajek_file = os.path.join(work_dir, "baited_adj_duplets.net")
    nx.write_pajek(graph, pajek_file)
    
    # duplets
    print "Loading duplets"
    graph, duplets = load_duplets_graph()
    print "Writing duplets"
    linked_list_file = os.path.join(work_dir, "baited_duplets.txt")
    with open(linked_list_file, "w") as outf:
        outf.write("#From\tTo\tWeight\n")
        for edge in graph.edges_iter(data=True):
            _from = edge[0]
            _to = edge[1]
            _weight = edge[2]["weight"]
            outf.write("%s\t%s\t%s\n" % (_from, _to, _weight))
    pajek_file = os.path.join(work_dir, "baited_duplets.net")
    nx.write_pajek(graph, pajek_file)


if __name__=="__main__":

    work_dir = os.path.join(data_path, 'prok1402/graph/graph_files')
    generate_graphs_from_baited_loci()
