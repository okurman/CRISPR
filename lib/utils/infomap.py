__author__ = 'Sanjarbek Hudaiberdiev'

import os
import sys
import configparser
from collections import defaultdict

###############################################################
config_file = os.path.join(os.path.expanduser('~'), 'paths.cfg')
cfg = configparser.ConfigParser()
cfg.read(config_file)

code_path = cfg.get('NewSystems', 'code_path')
data_path = cfg.get('NewSystems', 'data_path')
sys.path.append(code_path)

###############################################################


def parse_tree_file(file_name, node_names=True, node_desc_file=None):

    if not node_desc_file:
        node_desc_file = "/panfs/pan1/patternquest/Projects/NewSystems/data/prok1603/graph/infomap/merged/nodes.txt"

    module2members = {}
    member2modules = defaultdict(set)

    for line in open(file_name):

        if line.startswith("#"):
            continue

        parts = line.split()

        cluster_desc = parts[0].split(":")
        member = parts[2][1:-1]

        assert member == parts[3]

        module_name = "-".join(cluster_desc[:-1])

        if module_name in module2members:
            module2members[module_name].append(member)
        else:
            module2members[module_name] = [member]

        member2modules[member].update([module_name])
        # if member in member2modules:
        #     member2modules[member].update([module_name])
        # else:
        #     member2modules[member] = set([module_name])

    if node_names:

        new_member2modules = {}

        node_id2node_name = {}
        if not node_desc_file:
            raise ValueError("Provide the file with node descriptions")

        for l in open(node_desc_file):
            parts = l.rstrip().split()
            node_id2node_name[parts[0]] = parts[1]

        for k in module2members.keys():
            module2members[k] = [node_id2node_name[v] for v in module2members[k]]

        for (k,v) in member2modules.items():
            new_member2modules[node_id2node_name[k]] = v

        member2modules = new_member2modules

    return module2members, member2modules


def tree2tp(tree_file, tp_file):

    module2members, member2modules = parse_tree_file(tree_file, node_names=False)

    with open(tp_file, "w") as outf:
        for (community_name, members) in sorted([(k,v) for (k,v) in module2members.items()], key=lambda x: len(x[1]), reverse=True):
            outf.write("#Community %s Size: %d\n" % (community_name, len(members)))
            outf.write("%s\n" % " ".join(members))



if __name__=="__main__":

    f = "/panfs/pan1/patternquest/Projects/NewSystems/data/prok1603/graph/infomap/merged/out/1/link_list.tree"
    node_desc_file = "/panfs/pan1/patternquest/Projects/NewSystems/data/prok1603/graph/infomap/merged/nodes.txt"
    g2m, sg2m = parse_tree_file(f, node_desc_file)


