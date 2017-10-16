#!/home/hudaiber/env2.7/bin/python
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

import lib.utils.tools as t
from itertools import combinations, product
import parser
import networkx as nx

PROK1603_WEIGHT = 4930.0


def add_to_graph_old_version(pair2weight, pair2count, profile2weight, genes, weight):

    previous_profiles = genes[0].profiles

    if previous_profiles and len(previous_profiles) > 1:
        domain_duplets = list(combinations(previous_profiles, 2))

        for duplet in domain_duplets:
            [kplet_1, kplet_2] = sorted(duplet)
            key = "%s-%s" % (kplet_1, kplet_2)
            pair2weight[key] += weight
            pair2count[key] += 1

    if previous_profiles:
        for profile in previous_profiles:
            profile2weight[profile] += weight

    for gene in genes[1:]:
        cur_profiles = gene.profiles

        if not previous_profiles:
            previous_profiles = cur_profiles

            if cur_profiles:
                for profile in cur_profiles:
                    profile2weight[profile] += weight

            continue

        if not cur_profiles:
            continue

        if len(cur_profiles) > 1:
            domain_duplets = list(combinations(previous_profiles, 2))

            for duplet in domain_duplets:
                [kplet_1, kplet_2] = sorted(duplet)
                key = "%s-%s" % (kplet_1, kplet_2)
                pair2weight[key] += weight
                pair2count[key] += 1

        for profile in cur_profiles:
            profile2weight[profile] += weight

        adjacent_duplets = list(product(previous_profiles, cur_profiles))

        for duplet in adjacent_duplets:
            [kplet_1, kplet_2] = sorted(duplet)
            key = "%s-%s" % (kplet_1, kplet_2)
            pair2weight[key] += weight
            pair2count[key] += 1

        previous_profiles = cur_profiles


def add_to_graph(pair2weight, pair2count, profile2weight, genes, weight):
    """Edit: The earlier version of this function was corrected with this one after discussion with Yuri.
        The difference is in the treatment of multidomain proteins (genes). My reasoning was, that as long as inside one protein,
        the order of domains should not play role, hence the the procedure took every combination of domains as adjacent. As a consequence,
        neighboring multidomain genes produced adjecent links produced by cartesian product of domains (in addition to
        within-gene combinations)
        This function treats all such cases as one plain string.
        Edit 2:
        if the multidomain protein genes are in opposite strands with the previous neighbor, revert the order of domains.
        """

    all_profiles = ""

    reference_strand = genes[0].strand

    for gene in genes:

        if gene.profiles:

            if type(gene.profiles) == str:

                all_profiles += " %s" % gene.profiles

            elif type(gene.profiles) == list:

                _profiles = gene.profiles if gene.strand == reference_strand else reversed(gene.profiles)
                all_profiles += " %s" % " ".join(_profiles)

            else:

                raise Exception("The gene.profiles is neither list nor a string")

    all_profiles = all_profiles.strip().split()

    kplet_2 = None

    for i in range(0, len(all_profiles)-1):

        duplet = all_profiles[i:i+2]
        [kplet_1, kplet_2] = sorted(duplet)

        assert kplet_2 and kplet_1

        key = "%s-%s" % (kplet_1, kplet_2)
        pair2weight[key] += weight
        pair2count[key] += 1

        profile2weight[kplet_1] += weight

    if kplet_2:
        profile2weight[kplet_2] += weight


def graph_from_prok1603():

    genome2weight = parser.map_org2weight()

    pty_path = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/patternquest/data/Prok1603/pty"
    work_dir = os.path.join(data_path, 'prok1603/graph/')

    crisprs = parser.get_crispr_annotations()

    pair2weight = defaultdict(float)
    pair2count = defaultdict(int)
    profile2weight = defaultdict(float)

    cnt = 1

    weighted_sum = 0
    total_weight = 0

    for dir in os.listdir(pty_path):

        if dir not in genome2weight:
            continue

        print cnt, dir
        _weight = genome2weight[dir]

        source_files = [f for f in os.listdir(os.path.join(pty_path, dir)) if f.endswith(".pty2")]

        for source_file in source_files:

            source = os.path.splitext(source_file)[0]

            genes = t.parse_pty_file(os.path.join(pty_path, dir, source_file))

            if not genes:
                continue

            genes += crisprs[source]
            genes = sorted(genes)

            add_to_graph(pair2weight, pair2count, profile2weight, genes, _weight)

            profiles_num = 0

            for gene in genes:
                _p = gene.profiles
                if _p:
                    profiles_num += len(_p)
                else:
                    profiles_num += 1

            if profiles_num > 0:
                weighted_sum += 2 * _weight / profiles_num
                total_weight += _weight

        cnt += 1
        # if cnt == 10:
        #     break

    print "Weighted average of 2/N", weighted_sum/total_weight
    # Result of the above is: 0.0451314274706

    G = nx.Graph()

    for k,v in pair2weight.items():
        p1, p2 = k.split("-")
        G.add_edge(p1, p2, weight=v)

        _count = pair2count[k]
        G.add_edge(p1, p2, count=_count)

    graph_file = os.path.join(work_dir, "adj_graph.p")
    print "Writing to file:", graph_file
    nx.write_gpickle(G, graph_file)

    graph_file = os.path.join(work_dir, "adj_graph.p.bz2")
    print "Writing to file:", graph_file
    t.dump_compressed_pickle(graph_file, G)
    # nx.write_gpickle(G, graph_file)

    # profile_weight_file = os.path.join(work_dir, 'profile_weight.txt')
    # print "Writing to file:", profile_weight_file
    # with open(profile_weight_file, 'w') as outf:
    #     for (k,v) in sorted(profile2weight.items(), key= lambda x: x[1], reverse=True):
    #         outf.write("%s\t%f\n" % (k, v))


if __name__=="__main__":

    graph_from_prok1603()