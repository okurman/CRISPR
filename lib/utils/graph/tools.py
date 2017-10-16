import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import os

def log_binning(data, bin_count=50):

    data = np.array(data)
    ix = data > 0

    max_exp = np.log10(max(data))
    min_exp = np.log10(min(data[ix]))

    bins = np.logspace(min_exp, max_exp, num=bin_count)

    H = np.histogram(data, bins)
    h_y = H[0]
    h_x = H[1]

    h_x_size = h_x[1:] - h_x[:-1]
    h_x = h_x[:-1] + h_x_size / 2
    h_y = h_y / h_x_size
    ix = h_y > 0

    return h_x[ix], h_y[ix]


def linear_binning(data, bin_count=50):
    H = np.histogram(data, bin_count)
    h_y = H[0]
    h_x = H[1]
    bin_size = h_x[1] - h_x[0]
    h_x = h_x[:-1] + bin_size / 2
    ix = h_y > 0

    return h_x[ix], h_y[ix]


def degree_distributions(G, work_dir=None):

    bin_count = 100
    plt.figure()
    degrees = np.array(G.degree(weight='weight').values())

    h_x, h_y = linear_binning(degrees, bin_count=bin_count)
    # h_x, h_y = log_binning(degrees, bin_count=bin_count)
    plt.scatter(h_x, h_y, label=r'$f(d)$')

    plt.xscale('log')
    plt.yscale('log')

    x = np.log10(h_x)
    y = np.log10(h_y)
    [a,b] = np.polyfit(x, y, 1)
    Y = 10**b * h_x**a

    plt.plot(h_x, Y, c='r', linestyle='--', label=r'$\alpha = %.2f$' % a)

    plt.xlabel("Degree")
    plt.ylabel("Frequency")
    plt.title("Degree distribution")
    plt.legend()

    if work_dir:
        image_file = os.path.join(work_dir, 'degree_distribution.png')
        plt.savefig(image_file)
    else:
        plt.show()


def clustering_coefficients(G, work_dir=None):

    node_names = []
    degrees = []
    cc = []
    cnt = 1
    for node in G.nodes():

        if cnt % 10 == 0:
            print cnt
        cnt += 1

        _d = G.degree(node)
        if _d == 1:
            continue

        _d = G.degree(node, weight="weight")
        _c = nx.clustering(G, node, weight="weight")

        degrees.append(_d)
        cc.append(_c)
        node_names.append(node)

    # plt.figure()
    # plt.xscale('log')
    # plt.yscale('log')
    #
    # plt.scatter(degrees, cc, c='g', label=r'$f(C(w))$')
    # plt.show()
    return degrees, cc, node_names
    # x = np.log10(x_lin)
    # y = np.log10(y_lin)
    #
    # [a, b] = np.polyfit(x, y, 1)
    # Y = 10**b * x_lin**a
    # plt.plot(x_lin, Y, c='r', linestyle='--', label=r'$\alpha = %.2f$' % a)
    #
    # plt.xlabel("C(k)")
    # plt.ylabel("Frequency")
    # plt.title("Clustering coefficients")
    # plt.legend()
    #
    # if work_dir:
    #     image_file = os.path.join(work_dir, 'clustering_coefficients.png')
    #     print "Saving to", image_file
    #     plt.savefig(image_file)
    # else:
    #     plt.show()


def clustering_coefficient_distribution(G, work_dir=None):
    c = nx.clustering(G, weight="weight")
    cx = c.values()

    plt.figure()
    plt.xscale('log')
    plt.yscale('log')

    x_lin, y_lin = linear_binning(cx)
    plt.scatter(x_lin, y_lin, c='g', label=r'$f(C(w))$')

    x = np.log10(x_lin)
    y = np.log10(y_lin)

    [a, b] = np.polyfit(x, y, 1)
    Y = 10 ** b * x_lin ** a
    plt.plot(x_lin, Y, c='r', linestyle='--', label=r'$\alpha = %.2f$' % a)

    plt.xlabel("C(k)")
    plt.ylabel("Frequency")
    plt.title("Clustering coefficients")
    plt.legend()

    if work_dir:
        image_file = os.path.join(work_dir, 'clustering_coefficients.png')
        print "Saving to", image_file
        plt.savefig(image_file)
    else:
        plt.show()


def subgraph(G, node_list, radius=0):

    _graph = nx.ego_graph(G, node_list[0], radius=radius)

    for node in node_list[1:]:

        if G.has_node(node):
            _graph = nx.compose(_graph, nx.ego_graph(G, node, radius=radius))

    return _graph


def subgraph_crosslink(G, node_list, remove_self_links=None):

    _graph = nx.ego_graph(G, node_list[0], radius=1)

    for node in node_list[1:]:
        if G.has_node(node):
            _graph = nx.compose(_graph, nx.ego_graph(G, node, radius=1))

    for node in _graph.nodes():
        if node not in node_list:
            _graph.remove_node(node)

    if remove_self_links:
        for edge in _graph.edges():
            (n1, n2) = edge
            _graph.remove_edge(n1, n2)

    return _graph


def write_edgelist(G, file_name):

    with open(file_name, "w") as outf:
        for edge in G.edges(data=True):
            (p1, p2, data) = edge
            outf.write("%s\t%s\t%s\n" % (p1, p2, data["weight"]))


def read_edgelist(file_name):

    G = nx.Graph()

    for l in open(file_name):

        if l.startswith("#"):
            continue
        p1, p2, weight = l.rstrip().split("\t")
        weight = float(weight)

        G.add_edge(p1, p2, weight=weight)

    return G