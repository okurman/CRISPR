__author__ = 'hudaiber'

import sys
import os

sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))

import dm_tools as dt
import scipy.spatial.distance as ssd
import numpy as np
import os


class LocusSkeleton(object):

    # merged_files = []
    # merged_base_files = []

    def overlaps(self, other):

        if self.organism != other.organism and self.source != other.source:

            return False

        # Assumes that the genes in the loci are sorted!
        dist1 = other.genes[0].pFrom - self.genes[-1].pTo
        dist2 = self.genes[0].pFrom - other.genes[-1].pTo

        return True if dist1*dist2>=0 else False


    def merge(self, other):

        self.merged_files.update(other.merged_files)
        self.merged_base_files.update(other.merged_base_files)

        pool = self.genes + other.genes
        seen = set()
        self.genes = [gene for gene in pool if gene.gid not in seen and not seen.add(gene.gid)]

    @staticmethod
    def calculate(first, second):

        score_intersection = sum([0.5 if '-' in term else 1 for term in first.intersection(second)])
        score_union        = sum([0.5 if '-' in term else 1 for term in first.union(second)])

        return score_intersection / score_union


    def score(self, other):

        ff = self.calculate(self.forward_set, other.forward_set)
        fr = self.calculate(self.forward_set, other.reverse_set)
        rf = self.calculate(self.reverse_set, other.forward_set)
        rr = self.calculate(self.reverse_set, other.reverse_set)

        return max(ff, fr, rf, rr)




class Locus(LocusSkeleton):

    def __init__(self, file_name, annotation_map=None, file_format='pty', profile2gene=None):

        self.file_name = file_name
        self.base_file_name = os.path.basename(file_name)
        if file_format=='pty':
            _genes = dt.get_pty_file(file_name, annotation_map=annotation_map)
        else:
            _genes = dt.get_pty_file_generic(file_name, profile2gene=profile2gene)

        locus_file = open(file_name)

        type_line = locus_file.readline()

        if type_line.startswith('#type:'):
            self.crispr_type = type_line.split(':')[1].strip()
        else:
            self.crispr_type = None

        self.summary_line = locus_file.readline()

        assert self.summary_line.startswith("===")

        locus_file.close()

        self.genes = _genes
        self.profiles = set(profile for g in _genes for profile in g.cogid.split() if profile != '')
        self.gene_names = set(gene_name for g in _genes for gene_name in g.gene_name.split(',') if gene_name != '')
        self.clusters = set(gene.cluster_id for gene in _genes)
        self.organism = _genes[0].organism
        self.source = _genes[0].src

        _forward = set()
        _reverse = set()

        for i in range(len(self.genes)):
            _gene = self.genes[i]
            for _cogid in _gene.cogid.split(','):
                _forward.update((_cogid,))
                if i == len(self.genes)-1:
                    continue
                _next_gene = self.genes[i+1]
                for _next_cogid in _next_gene.cogid.split(','):
                    _forward.update(("%s-%s" % (_cogid, _next_cogid),))

        self.genes.sort(reverse=True)

        for i in range(len(self.genes)):
            _gene = self.genes[i]
            for _cogid in _gene.cogid.split(','):
                _reverse.update((_cogid,))
                if i == len(self.genes)-1:
                    continue
                _next_gene = self.genes[i+1]
                for _next_cogid in _next_gene.cogid.split(','):
                    _reverse.update(("%s-%s" % (_cogid, _next_cogid),))

        self.forward_set = _forward
        self.reverse_set = _reverse

        self.feature_weights = None
        self.feature_labels  = None


    @staticmethod
    def calculate(first, second):

        score_intersection = sum([0.5 if '-' in term else 1 for term in first.intersection(second)])
        score_union        = sum([0.5 if '-' in term else 1 for term in first.union(second)])

        return score_intersection / score_union


    def score(self, other):

        ff = self.calculate(self.forward_set, other.forward_set)
        fr = self.calculate(self.forward_set, other.reverse_set)
        rf = self.calculate(self.reverse_set, other.forward_set)
        rr = self.calculate(self.reverse_set, other.reverse_set)

        return max(ff, fr, rf, rr)


    def set_features(self, feature_labels, feature_weights):

        self.feature_labels = feature_labels
        self.feature_weights = feature_weights[:]

        for i in range(len(feature_labels)):
            if feature_labels[i] not in self.clusters:
                self.feature_weights[i]=0

        self.feature_weights = np.asarray(self.feature_weights, dtype=float)


    def cosine_distance(self, other):

        if self.feature_labels==None or self.feature_weights==None:
            print("Set feature weights AND feature labels for loci before invoking cosine_distance")
            sys.exit()

        # if self.feature_weights or other.feature_weights is all zeros
        if not np.any(self.feature_weights) or not np.any(other.feature_weights):
            return 1.0
        distance = ssd.cosine(self.feature_weights, other.feature_weights)
        if distance < 1e-8:
            distance = 0
        # distance = np.dot(self.feature_weights, other.feature_weights)
        return distance


class BasicLocus(LocusSkeleton):

    def __init__(self, file_name, annotation_map=None, file_format='pty', profile2gene=None):

        self.file_name = file_name
        self.base_file_name = os.path.basename(file_name)

        self.merged_files = set([self.file_name])
        self.merged_base_files = set([self.base_file_name])

        if file_format=='pty':
            _genes = dt.get_pty_file(file_name, annotation_map=annotation_map, profile2gene=profile2gene)
        else:
            _genes = dt.get_pty_file_generic(file_name, profile2gene=profile2gene)

        # print "File name:", file_name, "Genes size:",len(_genes)

        self.genes = _genes
        self.profiles = set(profile for g in _genes for profile in g.cogid.split(',') if profile != '')

        self.gene_names = set(gene_name for g in _genes for gene_name in g.gene_name.split(',') if gene_name != '')
        self.organism = _genes[0].organism
        self.source = _genes[0].src

        _forward = set()
        _reverse = set()

        for i in range(len(self.genes)):
            _gene = self.genes[i]
            for _cogid in _gene.cogid.split(','):
                _forward.update((_cogid,))
                if i == len(self.genes)-1:
                    continue
                _next_gene = self.genes[i+1]
                for _next_cogid in _next_gene.cogid.split(','):
                    _forward.update(("%s-%s" % (_cogid, _next_cogid),))

        self.genes.sort(reverse=True)

        for i in range(len(self.genes)):
            _gene = self.genes[i]
            for _cogid in _gene.cogid.split(','):
                _reverse.update((_cogid,))
                if i == len(self.genes)-1:
                    continue
                _next_gene = self.genes[i+1]
                for _next_cogid in _next_gene.cogid.split(','):
                    _reverse.update(("%s-%s" % (_cogid, _next_cogid),))

        self.forward_set = _forward
        self.reverse_set = _reverse

        self.feature_weights = None
        self.feature_labels  = None


class CrisprLocus(LocusSkeleton):

    def __init__(self, organism, source, genes, status, type, complete=None, id=None):

        self.organism = organism
        self.source = source
        self.status = status
        self.genes = genes
        self.type = type
        self.complete = complete
        self.id = id
        self.cas4_genes = [g for g in genes if 'cas4' in g.gene_name.lower()]
        self.gis = set([g.gid for g in genes if 'cas4' in g.gene_name.lower()])

        # self.has_cas4 = False
        # for _gene in self.genes:
        #
        #     if 'cas4' in _gene.gene_name.lower():
        #         self.has_cas4=True
        #         break

        self.gis = [g.gid for g in genes]

        _seed_gene = [_gene for _gene in self.genes if _gene.is_seed]
        self.seed_gene = _seed_gene[0] if _seed_gene else None
