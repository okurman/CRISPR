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


import networkx as nx
import matplotlib.pyplot as plt
from scipy.cluster._hierarchy import linkage
from numpy import linalg
import numpy as np
# import powerlaw
import lib.utils.graph.tools as gt
import lib.utils.tools as t
from lib.utils import infomap

PROK1603_WEIGHT = 4930.0
RAND_ADJ_PROB = 0.0451314274706

uvrd_profiles = [ "COG0210", "PRK11773", "TIGR01075", "cl22977", "cl26202", "cl28312", "pfam00580", "pfam13361",
                  "pfam13538" ]
cas4_profiles = ["cd09637", "cd09659", "cls000170", "COG1468", "COG4343", "pfam01930", "pfam06023"]

profile2gene = t.map_cdd_defense2gene_name()
work_dir = os.path.join(data_path, "prok1603/graph")

def graph_preprocessing_with_prob_expectation(work_dir, save_file=None):

    graph_file = os.path.join(work_dir, "adj_graph.p")

    G = nx.read_gpickle(graph_file)
    print "Raw graph size:", G.size()

    profile2prob = {l.split()[0]: float(l.split()[1]) for l in open(os.path.join(work_dir, 'profile_weight.txt'))}

    for k in profile2prob:
        profile2prob[k] /= PROK1603_WEIGHT

    total_edges = 0
    removed_edges = 0

    for edge in G.edges(data=True):
        nodes = edge[:2]
        _weight = edge[2]['weight']/PROK1603_WEIGHT
        _expected_weight = profile2prob[nodes[0]] * profile2prob[nodes[1]] * RAND_ADJ_PROB

        total_edges += 1
        if _weight < _expected_weight:
            G.remove_edge(*nodes)
            removed_edges +=1
        else:
            G[nodes[0]][nodes[1]]['weight'] = _weight

    print total_edges
    print removed_edges

    for node in G.nodes():
        if G.degree(node) == 0:
            print G[node]

    print "Pre-processed graph size", G.size()

    if save_file:
        print "Saving to", save_file
        nx.write_gpickle(G,save_file)

    return G


def graph_preprocessing_with_counts(G_input, weight_thr=2, save_file=None, small_components_dir=None):

    G = G_input.copy()

    print "Raw graph size:", G.size()
    print "Raw graph nodes", G.number_of_nodes()

    for edge in G.edges(data=True):

        nodes = edge[:2]
        _weight = edge[2]['weight']
        # _count = edge[2]['count']

        if _weight < weight_thr or nodes[0] == nodes[1]:
            G.remove_edge(*nodes)

    print "Pre-processed graph size", G.size(), G.number_of_edges()
    print "Pre-processed graph nodes", G.number_of_nodes()
    # print nx.isolates(G)
    print "Removing the isolates and small components"

    nodes_to_remove = set()
    [nodes_to_remove.update(c) for c in nx.connected_components(G) if len(c) < 2]

    if small_components_dir:
        i = 1
        for c in nx.connected_components(G):
            if len(c) < 10:
                if len(c) == 1:
                    continue

                c_dir = os.path.join(small_components_dir, str(i+1))
                if not os.path.exists(c_dir):
                    os.mkdir(c_dir)
                subgraph = G.subgraph(c)

                gt.write_edgelist(subgraph, c_dir)

                i += 1

    for node in nodes_to_remove:
        G.remove_node(node)

    print "Pre-processed graph size", G.size(), G.number_of_edges()
    print "Pre-processed graph nodes", G.number_of_nodes()

    # assert len(list(nx.connected_components(G))) == 1

    if save_file:
        print "Saving to", save_file
        nx.write_gpickle(G, save_file)

    return G


def uvrd_cas4_links(G, graph_save_dir, radius=1):

    all_profiles = uvrd_profiles + cas4_profiles
    profile2gene = t.map_cdd_defense2gene_name()

    print "Extracting subgraphs"
    cas4_g = gt.subgraph(G, cas4_profiles, radius=radius)
    uvrd_g = gt.subgraph(G, uvrd_profiles, radius=radius)
    print "Joining subgraphs"
    joint_g = nx.compose(cas4_g, uvrd_g)
    print "Writing graph to files"

    nodes_to_write= set()

    with open(os.path.join(graph_save_dir, "edgelist.txt"),"w") as edge_outf:

        for edge in joint_g.edges(data=True):
            (p1, p2, data) = edge

            # if p1 in all_profiles and p2 in all_profiles:
            edge_outf.write("%s\t%s\t%f\n" % (p1, p2, data['weight']))
            nodes_to_write.update([p1, p2])

    with open(os.path.join(graph_save_dir, "nodes.txt"),"w") as node_outf:

        for p in nodes_to_write:
            if p in cas4_profiles:
                _type = 1
            elif p in uvrd_profiles:
                _type = 2
            else:
                _type = 0

            node_outf.write("%s\t%s\t%d\n" % (p, profile2gene[p], _type))

    return joint_g


def crispr_graph(G, graph_save_dir=None, radius=1):

    def_profs_file = "/panfs/pan1/patternquest/Projects/NewSystems/data/profiles/defenseProfiles.tab"
    cas_profiles = set()

    profile2cluster = {}

    for l in open(os.path.join(work_dir, "synonyms_in_graph.txt")):
        parts = l.strip().split("\t")
        cl_name = parts[0]
        members = parts[1].split()

        for member in members:
            profile2cluster[member] = cl_name

    for l in open(def_profs_file):

        if l.split()[1] != "CRISPR":
            continue

        profile = l.split()[0]

        if profile in profile2cluster:
            cas_profiles.update([profile2cluster[profile]])
        else:
            cas_profiles.update([profile])

    profile2gene = t.map_cdd_defense2gene_name()
    cas_g = gt.subgraph(G, list(cas_profiles), radius=radius)
    nodes_to_write = set()

    if graph_save_dir:

        with open(os.path.join(graph_save_dir, "edgelist.txt"), "w") as edge_outf:

            for edge in cas_g.edges(data=True):
                (p1, p2, data) = edge

                # if p1 in all_profiles and p2 in all_profiles:
                edge_outf.write("%s\t%s\t%f\n" % (p1, p2, data['weight']))
                nodes_to_write.update([p1, p2])

        with open(os.path.join(graph_save_dir, "nodes.txt"), "w") as node_outf:

            for p in nodes_to_write:
                if p in cas4_profiles:
                    _type = 1
                elif p in cas_profiles:
                    _type = 2
                else:
                    _type = 0

                node_outf.write("%s\t%s\t%d\n" % (p, profile2gene[p], _type))

    return cas_g


def merge_nodes(G, nodes, new_node):

    G.add_node(new_node)

    for node in nodes:
        neighbors = G[node].keys()
        for neighbor in neighbors:
            edge_to_rewire = G[node][neighbor]
            if neighbor in G[new_node]:
                for _key, _value in edge_to_rewire.items():
                    G[new_node][neighbor][_key] += _value
            else:
                G.add_edge(new_node, neighbor, attr_dict=edge_to_rewire)

        G.remove_node(node)


def merge_synonyms(G, synonym_rep_file=None):

    synonyms_file = os.path.join(data_path, "profiles/synonyms.txt")
    print "reading from:", synonyms_file

    graph_nodes = set(G.nodes())

    group2profiles = {}
    profile2group = {}

    merged_nodes = 0
    with open(synonyms_file) as inf:
        l = inf.readline()
        while l:
            profiles = l.split("\t")[1]
            profiles = profiles.split()
            
            profiles = list(graph_nodes.intersection(profiles))

            if len(profiles) < 2:
                l = inf.readline()
                continue

            size = len(profiles)
            new_node = "s_%d_%s" % (size, profiles[0])

            group2profiles[new_node] = profiles

            for pr in profiles:
                profile2group[pr] = new_node
            merged_nodes += size
            l = inf.readline()

    if synonym_rep_file:
        print "Writing synonym lookup to:", synonym_rep_file
        with open(synonym_rep_file, "w") as outf:
            for k, v in group2profiles.items():
                outf.write("%s\t%s\n" % (k, " ".join(v)))

    new_G = nx.Graph()

    print "Total number of nodes: ", len(graph_nodes)
    unmerged_nodes = graph_nodes.difference(profile2group.keys())
    print "%d nodes are unmerged." % len(unmerged_nodes)
    for node in unmerged_nodes:
        new_G.add_node(node)
    print "Nodes merged to groupds %d => %d" % (merged_nodes, len(group2profiles))
    for group in group2profiles:
        new_G.add_node(group)
    print "Number of nodes in new graph: ", new_G.number_of_nodes()

    groups = group2profiles.keys()

    rewired_edges = 0

    # Wiring the groups among each other
    for i in range(0, len(groups)-1):
        outer_group = groups[i]
        outer_group_profiles = group2profiles[outer_group]
        outer_group_neighbors = set([v for node in group2profiles[outer_group] for v in G.neighbors(node)])

        for j in range(i+1, len(groups)):
            inner_group = groups[j]
            connected_nodes = outer_group_neighbors.intersection(group2profiles[inner_group])

            if not connected_nodes:
                continue

            total_weight = 0
            total_count = 0

            for node1 in connected_nodes:
                for node2 in group2profiles[outer_group]:

                    if G.has_edge(node1, node2):
                        data = G[node1][node2]
                        total_weight += data["weight"]
                        total_count += data["count"]

                        rewired_edges += 1

            new_G.add_edge(outer_group, inner_group, weight=total_weight, count=total_count)

        for i in range(len(outer_group_profiles)-1):
            for j in range(i,len(outer_group_profiles)):
                if G.has_edge(outer_group_profiles[i], outer_group_profiles[j]):
                    rewired_edges += 1

    outer_group_profiles = group2profiles[groups[-1]]
    for i in range(len(outer_group_profiles) - 1):
        for j in range(i, len(outer_group_profiles)):
            if G.has_edge(outer_group_profiles[i], outer_group_profiles[j]):
                rewired_edges += 1

    # wiring the unmerged nodes

    for unmerged_node in unmerged_nodes:

        to_be_grouped = []

        for (_, neighbor, data) in G.edges(unmerged_node, data=True):
            if neighbor in profile2group:
                to_be_grouped.append(neighbor)
            else:
                new_G.add_edge(unmerged_node, neighbor, weight=data["weight"], count=data["count"])
                rewired_edges += 1
                continue

        linked_groups = set([profile2group[n] for n in to_be_grouped])

        for group in linked_groups:
            total_count = 0
            total_weight = 0

            linked_nodes = set(to_be_grouped).intersection(group2profiles[group])
            for node in linked_nodes:
                data = G[unmerged_node][node]
                total_count += data["count"]
                total_weight += data["weight"]
                rewired_edges += 1

            new_G.add_edge(unmerged_node, group, weight=total_weight, count=total_count)

    print "Size of original graph: %d, number of rewired edges: %d" % (G.size(), rewired_edges)

    return new_G


def write_linked_list_profile_nodes(G, save_dir):

    node2id = {node: id+1 for (id, node) in enumerate(G)}

    print "Saving link lists with node IDs"
    f = os.path.join(save_dir, "nodes.txt")
    print f
    with open(f, "w") as outf:
        for node, id in node2id.items():
            outf.write("%d\t%s\t%s\n" % (id, node, profile2gene[node]))

    f = os.path.join(save_dir, "link_list.txt")
    print f
    with open(f, "w") as outf:
        for edge in G.edges(data=True):
            node1, node2, data = edge

            outf.write("%d\t%d\t%f\n" % (node2id[node1], node2id[node2], data["weight"]))


def write_linked_list_id_nodes(G, save_dir):

    f = os.path.join(save_dir, "link_list.txt")

    with open(f, "w") as outf:
        for edge in G.edges(data=True):
            node1, node2, data = edge
            outf.write("%d\t%d\t%f\n" % (node1+1, node2+1, data["weight"]))


def process_plain_graph():

    print "Loading the graph"
    graph_file = os.path.join(work_dir, "adj_graph.p")
    G = nx.read_gpickle(graph_file)

    prp_save_file = os.path.join(work_dir, "adj_graph_raw_prp.p")
    small_components_dir = os.path.join(work_dir, "small_components_raw")
    G = graph_preprocessing_with_counts(G_input=G, save_file=prp_save_file, small_components_dir=small_components_dir)

    link_list_dir = os.path.join(work_dir, "infomap/plain/")
    write_linked_list(G, link_list_dir)


def process_merged_synonym_graph():

    print "Loading the graph"
    graph_file = os.path.join(work_dir, "adj_graph.p")
    G = nx.read_gpickle(graph_file)

    link_list_dir = os.path.join(work_dir, "infomap/merged/")
    synonym_file = os.path.join(link_list_dir, "synonyms_in_graph.txt")
    prp_save_file = os.path.join(work_dir, "adj_graph_merged_prp.p")
    small_components_dir = os.path.join(work_dir, "small_components_merged")

    print "Merging the synonyms"
    G = merge_synonyms(G, synonym_file)

    G = graph_preprocessing_with_counts(G_input=G, save_file=prp_save_file, small_components_dir=small_components_dir)

    write_linked_list(G, link_list_dir)

    return G


def clustering_score_matrix(clusters_dir, m_file):

    #TODO: Get it from the nodes_file
    N = 12105
    M = np.zeros((N,N))

    cnt = 1
    for cl_dir in os.listdir(clusters_dir):
        tree_file = os.path.join(clusters_dir, cl_dir, "link_list.tree")

        _, member2modules = infomap.parse_tree_file(tree_file, node_names=False)

        print cnt, tree_file

        for i in range(N):
            for j in range(i, N):
                i1 = str(i+1)
                i2 = str(j+1)
                common_modules = member2modules[i1].intersection(member2modules[i2])
                M[i,j] += len(common_modules)

        cnt += 1

    np.savez_compressed(m_file, matrix=M)


def clustering_score_matrix_single_file(clusters_dir, m_file, score_normalizer="min"):

    #TODO: Get it from the nodes_file
    N = 12105
    M = np.zeros((N,N))

    tree_file = os.path.join(clusters_dir, "link_list.tree")
    _, member2modules = infomap.parse_tree_file(tree_file, node_names=False)

    for i in range(N):
        for j in range(i, N):
            i1 = str(i+1)
            i2 = str(j+1)

            common_modules = len(member2modules[i1].intersection(member2modules[i2]))

            if common_modules == 0:
                _score = 0
            else:
                if score_normalizer == "min":
                    _normalizer = min(len(member2modules[i1]), len(member2modules[i2]))
                elif score_normalizer == "max":
                    _normalizer = max(len(member2modules[i1]), len(member2modules[i2]))

                _score = float(common_modules) / _normalizer

            M[i,j] += _score

    print "Saving to file:", m_file
    t.dump_compressed_pickle(m_file, M)
    # np.savez_compressed(m_file, matrix=M)


def run_single_matrix_for_sge():

    if len(sys.argv) < 3:
        print "graph_analysis.py iteration normalizer run_no"
        sys.exit()

    iteration = sys.argv[1]
    normalizer = sys.argv[2]
    run_no = sys.argv[3]

    communities_dir = os.path.join(data_path,
                    "prok1603/graph/infomap/merged/consensus_clustering/multi_level/norm_%s/%s/infomap_runs/%s" %
                                   (normalizer, iteration, run_no))

    if not os.path.exists(communities_dir):
        print "Path doesn't exist: "
        print communities_dir
        sys.exit()

    m_file = os.path.join(communities_dir, "m.p.bz2")

    print "Normalizer: ", normalizer
    print "Iteration #: ", iteration
    print "Run #: ", run_no
    print "Dir: ", communities_dir
    print "Save file: ", m_file

    clustering_score_matrix_single_file(communities_dir, m_file, score_normalizer=normalizer)


def run_single_matrix_for_sge_tmp():

    if len(sys.argv) < 3:
        print "graph_analysis.py iteration normalizer run_no"
        sys.exit()

    iteration = sys.argv[1]
    normalizer = sys.argv[2]
    run_no = sys.argv[3]

    communities_dir = os.path.join(data_path,
                    "prok1603/graph/infomap/merged/consensus_clustering/multi_level/norm_%s/%s/runs/%s" %
                                   (normalizer, iteration, run_no))

    if not os.path.exists(communities_dir):
        print "Path doesn't exist: "
        print communities_dir
        sys.exit()

    m_file = os.path.join(communities_dir, "m.p.bz2")

    print "Normalizer: ", normalizer
    print "Iteration #: ", iteration
    print "Run #: ", run_no
    print "Dir: ", communities_dir
    print "Save file: ", m_file

    clustering_score_matrix_single_file(communities_dir, m_file, score_normalizer=normalizer)



def merge_score_matrices():

    if len(sys.argv) < 3:
        print "graph_analysis.py iteration normalizer "
        sys.exit()

    iteration = sys.argv[1]
    normalizer = sys.argv[2]

    current_iter_dir = os.path.join(data_path,
                        "prok1603/graph/infomap/merged/consensus_clustering/multi_level/norm_%s/%s/" %
                        (normalizer, iteration))
    # source_dir = os.path.join(current_iter_dir, "infomap_runs")

    # if not os.path.exists(source_dir):
    #     print "Path doesn't exist: ", source_dir
    #     sys.exit()
    #
    # N = 12105
    # M = np.zeros((N, N))
    # print "Processing folder ", source_dir
    # for i, d in enumerate(os.listdir(source_dir)):
    #     f = os.path.join(source_dir, d, "m.p.bz2")
    #
    #     if not os.path.exists(f):
    #         print "Not found: ", f
    #         continue
    #
    #     print i, f
    #     _m = t.load_compressed_pickle(f)
    #     M += _m

    save_file = os.path.join(data_path,
                 "prok1603/graph/infomap/merged/consensus_clustering/multi_level/norm_%s/%s/" % (normalizer, iteration),
                 "consensus_matrix.p.bz2")

    # print "Saving to ", save_file
    # t.dump_compressed_pickle(save_file, M)

    # next_iter_dir =  os.path.join(data_path,
    #                 "prok1603/graph/infomap/merged/consensus_clustering/multi_level/norm_%s/%s/" %
    #                                (normalizer, str(int(iteration)+1)))

    # G, M = extract_consensus_scheme(save_file, t_thr=0.5)

    # print "Writing the concensus graph for next iteration:"
    # print next_iter_dir
    # write_linked_list_id_nodes(G, next_iter_dir)

    print "Calculating the eigenvalues"
    M = t.load_compressed_pickle(save_file)
    m = M / np.sum(np.diag(M))
    w, _ = linalg.eig(m)

    f = os.path.join(current_iter_dir, "eigenvalues.p.bz2")
    print "Saving eigenvalues to ", f
    t.dump_compressed_pickle(f, w)

    return M


def merge_score_matrices_tmp():

    if len(sys.argv) < 3:
        print "graph_analysis.py iteration normalizer "
        sys.exit()

    iteration = sys.argv[1]
    normalizer = sys.argv[2]

    current_iter_dir = os.path.join(data_path,
                        "prok1603/graph/infomap/merged/consensus_clustering/multi_level/norm_%s/%s/" %
                        (normalizer, iteration))
    source_dir = os.path.join(current_iter_dir, "runs")

    if not os.path.exists(source_dir):
        print "Path doesn't exist: ", source_dir
        sys.exit()

    # N = 12105
    # M = np.zeros((N, N))
    # print "Processing folder ", source_dir
    # for i, d in enumerate(os.listdir(source_dir)):
    #     f = os.path.join(source_dir, d, "m.p.bz2")
    #
    #     if not os.path.exists(f):
    #         print "Not found: ", f
    #         continue
    #
    #     print i, f
    #     _m = t.load_compressed_pickle(f)
    #     M += _m

    save_file = os.path.join(data_path,
                 "prok1603/graph/infomap/merged/consensus_clustering/multi_level/norm_%s/%s/" % (normalizer, iteration),
                 "consensus_matrix.p.bz2")

    # print "Saving to ", save_file
    # t.dump_compressed_pickle(save_file, M)
    M = t.load_compressed_pickle(save_file)
    next_iter_dir =  os.path.join(data_path,
                    "prok1603/graph/infomap/merged/consensus_clustering/multi_level/norm_%s/%s/" %
                                   (normalizer, str(int(iteration)+1)))

    G, M = extract_consensus_scheme(save_file, t_thr=0.5)

    print "Writing the concensus graph for next iteration:"
    print next_iter_dir
    write_linked_list_id_nodes(G, next_iter_dir)

    # print "Calculating the eigenvalues"
    #
    # m = M / np.sum(np.diag(M))
    # w, _ = linalg.eig(m)
    #
    # f = os.path.join(current_iter_dir, "eigenvalues.p.bz2")
    # print "Saving eigenvalues to ", f
    # t.dump_compressed_pickle(f, w)

    return M


def extract_consensus_scheme(score_matrix_file, t_thr=0.3):

    M = t.load_compressed_pickle(score_matrix_file)

    if M.diagonal().any():
        assert len(set(M.diagonal())) == 1
        np.fill_diagonal(M, 0)

    if not np.all(M == M.T):
        M += M.T

    if np.max(M) != 1:
        M /= np.max(M)

    M[M <= t_thr] = 0
    G = nx.from_numpy_matrix(M)

    return G, M


def cluster_with_upgma(M):

    inf_default = 50

    M += np.transpose(M)
    M /= np.max(M)
    M = np.negative(np.log(M))

    np.fill_diagonal(M, 0)

    inf_idx = np.isinf(M)
    M[inf_idx] = inf_default

    # M_array = ssd.squareform(M)
    # Z = linkage(M_array, method='average')

    return M


def iterate_consensus_clustering():

    source_dir = sys.argv[1]
    starting_ind = sys.argv[2]

    next_ind = str(int(starting_ind)+1)
    iter_dir = os.path.join(source_dir, next_ind)

    if not os.path.exists(iter_dir):
        os.mkdir(iter_dir)


def plot_iteration_eigenvalues():

    n = 5
    source_dir = "/panfs/pan1/patternquest/Projects/NewSystems/data/prok1603/graph/infomap/merged/consensus_clustering/multi_level/norm_max/"

    _f = os.path.join(source_dir, "1/eigenvalues.p.bz2")

    w = t.load_compressed_pickle(_f)
    w = np.sort(w)[::-1]
    plt.plot(w)
    plt.xlabel(r'$i$')
    plt.ylabel(r'$\lambda_i$')

    # for i in range(1,5):


def trace_iterations():

    n = 4
    source_dir = "/panfs/pan1/patternquest/Projects/NewSystems/data/prok1603/graph/infomap/merged/consensus_clustering/multi_level/norm_min/"
    t_thr=0.5

    nodes = np.zeros(n)
    edges = np.zeros(n)

    for i in range(n):
        curdir = os.path.join(source_dir, str(i+1))
        cons_file = os.path.join(curdir, "consensus_matrix.p.bz2")
        link_file = os.path.join(curdir, "link_list.txt")

        G_l = gt.read_edgelist(link_file)
        nodes[i] = G_l.number_of_nodes()
        edges[i] = G_l.number_of_edges()

    print nodes
    print edges


def generate_block_diagonal(score_file, tree_file):

    group2members, member2groups = infomap.parse_tree_file(tree_file)


def processing_order_track():

    work_dir = os.path.join(data_path, "prok1603/graph/")
    graph_file = os.path.join(work_dir, "adj_graph.p.bz2")
    synonym_def_file = os.path.join(work_dir, "synonyms_in_graph.txt")
    merged_file = os.path.join(work_dir, "adj_graph_merged.p.bz2")
    preprocessed_file = os.path.join(work_dir, "adj_graph_prp.p.bz2")

    # G = t.load_compressed_pickle(graph_file)
    # G = merge_synonyms(G, synonym_def_file)
    # t.dump_compressed_pickle(merged_file, G)
    # graph_preprocessing_with_counts(G_input=G, save_file=preprocessed_file)

    G = t.load_compressed_pickle(preprocessed_file)





if __name__ == "__main__":


    processing_order_track()

    # synonym_file = os.path.join(work_dir, "synonyms_in_graph.txt")
    # print "Merging the synonyms"
    # G = merge_synonyms(G, synonym_file)

    # preprocessed_file = os.path.join(work_dir, "adj_graph_merged_prp.p")
    # G = graph_preprocessing_with_counts(G_input=G, save_file=preprocessed_file, )
    # G = nx.read_gpickle(preprocessed_file)

    # save_dir = os.path.join(work_dir, "crispr_graph")
    # crispr_graph(G, save_dir)

    # link_list_dir = os.path.join(work_dir, "infomap/plain/")
    # write_linked_list(G, link_list_dir)
    # nx.write_pajek(G, os.path.join(link_list_dir, "graph.net"))
    # gt.degree_distributions(G)
    # gt.clustering_coefficients(G)

    # process_plain_graph()
    # process_merged_synonym_graph()

    # run_single_matrix_for_sge_tmp()
    # merge_score_matrices_tmp()
    # f = os.path.join(work_dir, 'infomap/merged/consensus_clustering/multi_level/norm_max/2/consensus_matrix.p.bz2')
    # save_dir = os.path.join(work_dir, "infomap/merged/consensus_clustering/multi_level/norm_max/3/")
    # G, M = extract_consensus_scheme(f, t_thr=0.5, save_dir=save_dir)

    # trace_iterations()