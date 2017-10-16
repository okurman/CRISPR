# coding: utf-8
get_ipython().magic(u'pwd ')
get_ipython().magic(u'll ')
import os
import sys
import configparser
from collections import defaultdict

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

uvrd_profiles = [ "COG0210", "PRK11773", "TIGR01075", "cl22977", "cl26202", "cl28312", "pfam00580", "pfam13361",
                  "pfam13538" ]
cas4_profiles = ["cd09637", "cd09659", "cls000170", "COG1468", "COG4343", "pfam01930", "pfam06023"]

profile2gene = t.map_cdd_defense2gene_name()
work_dir = os.path.join(data_path, "prok1603/graph/")
preprocessed_file = os.path.join(work_dir, "adj_graph_merged_prp.p")
G = nx.read_gpickle(preprocessed_file)

def_profs_file = "/panfs/pan1/patternquest/Projects/NewSystems/data/profiles/defenseProfiles.tab"
cas_profiles = [l.split()[0] for l in open(def_profs_file) if l.split()[1]=="CRISPR" ]
