__author__ = 'Sanjarbek Hudaiberdiev'

import sys
import os

sys.path.append(os.path.join(os.path.expanduser('~'), 'Projects/lib/BioPy/'))

# import global_variables as gv
import dm_tools as dt
import tools as t

try:
    gnm2weight = t.map_genome2weight()
except:
    gnm2weight = {}

try:
    file2org = t.map_file2organism()
except:
    file2org = {}


class ProfileCount(object):
    def __init__(self, count, weight, gid):
        self.count = count
        self.weight = weight
        self.gids = [gid]


class Neighborhood(object):
    def __init__(self, source_file):
        self.source_file = source_file
        self.genes = dt.get_pty_file(source_file)
        for gene in self.genes:
            gene.tag = 'neighborhood'
        self.flank_extension = None

    def extend_flanks(self, flank_size, pty_path, cdd_map=None):

        first = self.genes[0].gid
        last = self.genes[-1].gid
        upstream, downstream = [], []

        pty_genes = dt.get_pty_file(pty_path)
        pty_genes.sort()

        for i in range(len(pty_genes)):
            if pty_genes[i].gid == first:
                upstream = pty_genes[i-flank_size: i]
            if pty_genes[i].gid == last:
                downstream = pty_genes[i+1:i+flank_size+1]
                break

        for gene in upstream:
            gene.tag = 'flank'
        for gene in downstream:
            gene.tag = 'flank'

        self.genes = upstream + self.genes + downstream

        if cdd_map:
            for gene in self.genes:
                if gene.cogid in ['-',''] and gene.gid in cdd_map:
                    gene.cogid = cdd_map[gene.gid]

        self.flank_extension = True


class Kplet(object):

    def __init__(self, id, profiles, count=None, weight=None, files=None):
        self.profiles = set(profiles)
        self.k = len(profiles)
        self.id = id
        self.count = count
        self.files = sorted(files)
        self.weight = weight

    def __cmp__(self, other):
        if self.weight > other.weight:
            return 1
        elif self.weight < other.weight:
            return -1
        else:
            return 0


