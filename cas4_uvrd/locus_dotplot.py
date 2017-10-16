#!/home/hudaiber/env3.4/bin/python

__description__ = """ Create megablast based dotplots between two loci """

import os
import sys
import configparser
import subprocess as sp
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import SeqIO
import numpy as np
import argparse

config_file = os.path.join(os.path.expanduser('~'),'paths.cfg')
cfg=configparser.ConfigParser()
cfg.read(config_file)
code_path = cfg['NewSystems']['code_path'][1:-1]
data_path = cfg['NewSystems']['data_path'][1:-1]
sys.path.append(code_path)

# threshold level identity to be considered
THR_ID = 70
PADDING_PLOTS = 800
PADDING_CARTOONS = 500
RECTANGLE_WIDTH = 500
SHORT_HIT_LENGTH = 500

GENE_NAME_FONT_SIZE=20
LABEL_FONT_SIZE=20

defense_profiles_file="/net/frosty/vol/export1/proteome/defense/defenseProfiles.tab"

crispr_related_profiles=set([l.split("\t")[0] for l in open(defense_profiles_file) if "CRISPR" in l])


class Locus(object):

    def __init__(self, genome, contig, locus_from, locus_to, genes):

        self.genome = genome
        self.contig = contig
        self.l_from = locus_from
        self.l_to = locus_to
        self.genes = genes
        self.length = abs(locus_from - locus_to)
        self.label = "%s %s %d..%d" % (genome, contig, locus_from, locus_to)

    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}".format(self.genome, self.contig, self.l_from, self.l_to, self.length)


class Gene(object):

    def __init__(self, gi, p_from, p_to, strand, name, profile):

        self.gi = gi
        self.p_from = p_from
        self.p_to = p_to
        self.strand = strand
        self.name = name
        self.profile = set(profile.split(","))
        self.length = abs(self.p_from - self.p_to)

    def shift_back(self, dist):

        self.p_from = self.p_from - dist
        self.p_to = self.p_to - dist

    def __str__(self):
        terms = [self.gi, self.p_from, self.p_to, self.length, self.strand, self.profile, self.name]
        return "\t".join(str(t) for t in terms)


class Hit(object):

    def __init__(self,hit_line):

        parts = hit_line.split("\t")
        self.hit_line = hit_line
        self.perc_id = float(parts[2])
        self.length = int(parts[3])
        self.q_start = int(parts[6])
        self.q_end = int(parts[7])
        self.q_len = abs(self.q_end - self.q_start)
        self.s_start = int(parts[8])
        self.s_end = int(parts[9])
        self.s_len = abs(self.s_end - self.s_start)

        self.slope = self.s_len / float(self.q_len)


    def overlaps(self, other):

        # of the hit is on the same direction
        if other.s_start < other.s_end:

            dist1 = other.q_start - self.q_end
            dist2 = self.q_start - other.q_end
            q_overlaps = True if dist1*dist2>=0 else False

            dist1 = other.s_start - self.s_end
            dist2 = self.s_start - other.s_end
            s_overlaps = True if dist1*dist2>=0 else False

            return True if q_overlaps and s_overlaps else False

        # inverse hits.
        else:

            dist1 = other.q_start - self.q_end
            dist2 = self.q_start - other.q_end
            q_overlaps = True if dist1*dist2>=0 else False

            other_s_start = other.s_end
            other_s_end = other.s_start

            self_s_start = self.s_end
            self_s_end = self.s_start
            
            dist1 = other_s_start - self_s_end
            dist2 = self_s_start - other_s_end

            s_overlaps = True if dist1*dist2>=0 else False

            return True if q_overlaps and s_overlaps else False


    def add_padding(self, dist):

        self.q_start += dist
        self.q_end += dist
        self.s_start += dist
        self.s_end += dist


def load_locus(fname):

    lines = open(fname).readlines()
    parts = lines[0].split("\t")
    genome = parts[3]
    contig = parts[4]

    genes = []

    for l in lines[1:]:
        parts = l.split("\t")
        gi=parts[0]
        coords = parts[1]
        p_from, p_to = coords.split("..")
        p_from, p_to = int(p_from), int(p_to)
        strand = parts[2]
        profile = parts[6]
        name = parts[8]

        _gene = Gene(gi=gi, p_from=p_from, p_to=p_to, strand=strand, profile=profile, name=name)
        genes.append(_gene)

    locus_from = genes[0].p_from
    locus_to = genes[-1].p_to

    [gene.shift_back(locus_from) for gene in genes]
    locus = Locus(genome=genome, contig=contig, locus_from=locus_from, locus_to=locus_to, genes=genes)

    return locus


def retreive_locus_nt_sequence(locus, fna_file):

    blast_db_path = "/panfs/pan1/prokdata/db/all1603.nt"
    cmd_format = "blastdbcmd -db {blast_db} -entry {accession} -range {l_from}-{l_to} -long_seqids -out {out_file}"
    cmd = cmd_format.format(blast_db=blast_db_path,
                            accession=locus.contig,
                            l_from=locus.l_from,
                            l_to=locus.l_to,
                            out_file=fna_file)

    sp.call(cmd, shell=True)


def blast_between_loci(locus1_file, locus2_file, blast_file):

    cmd_fmt = "blastn -task megablast -query {locus1} -subject {locus2} -outfmt 7 -out {blast_file}"
    cmd = cmd_fmt.format(locus1=locus1_file,
                         locus2=locus2_file,
                         blast_file=blast_file)
    print(cmd)
    sp.call(cmd, shell=True)


def read_blast_hits(blast_file, average_locus_length):

    all_hits = [Hit(l) for l in open(blast_file) if not l.startswith("#")]
    distinct_hits = []

    for hit in all_hits:

        if hit.perc_id < THR_ID:
            continue

        if hit.length < SHORT_HIT_LENGTH:
            continue

        overlap_found = False

        for i, distinct_hit in enumerate(distinct_hits):
            if distinct_hit.overlaps(hit):

                distinct_hits[i] = distinct_hit if distinct_hit.length >= hit.length else hit
                overlap_found = True

        if not overlap_found:
            distinct_hits.append(hit)

    [hit.add_padding(PADDING_PLOTS) for hit in distinct_hits]

    return distinct_hits


def draw_x_cartoons(ax, locus):

    for gene in locus.genes:

        gene_name = gene.name.split(",")[0]

        if gene_name == "Unknown":
            gene_name = ""

        _rotation=90
        _fill = False
        _color = None

        if "CRISPR" in gene.profile:

            patch = patches.Rectangle((gene.p_from + PADDING_PLOTS, PADDING_CARTOONS-RECTANGLE_WIDTH/2), gene.length, RECTANGLE_WIDTH, color="k")
            _rotation=0
        else:
            
            if gene_name.lower() == 'cas4':
                _fill = True
                _color = 'g'
            elif gene.profile.intersection(crispr_related_profiles):
                _fill = True
                _color = 'b'
            
            if gene.strand == "+":

                x_pos = gene.p_from + PADDING_PLOTS
                y_pos = PADDING_CARTOONS
                dx = gene.length
                dy = 0

            elif gene.strand == "-":

                x_pos = gene.p_to + PADDING_PLOTS
                y_pos = PADDING_CARTOONS
                dx = -gene.length
                dy = 0


            patch = patches.Arrow(
                x_pos,
                y_pos,
                dx,
                dy,
                width=2000,
                fill=_fill,
                color=_color)

        ax.add_patch(patch)
        ax.annotate(gene_name, (gene.p_from + gene.length/2 + PADDING_PLOTS, PADDING_CARTOONS + RECTANGLE_WIDTH), 
            rotation=_rotation, 
            va='bottom', 
            ha='center',
            # size=8,
            weight='bold',
            fontsize=GENE_NAME_FONT_SIZE)

    ax.set_xlabel(locus.label, weight='bold', size=LABEL_FONT_SIZE)


def draw_y_cartoons(ax, locus):

    for gene in locus.genes:

        gene_name = gene.name.split(",")[0]

        if gene_name == "Unknown":
            gene_name = ""

        # _rotation=0
        _fill = False
        _color = None

        if "CRISPR" in gene.profile:

            patch = patches.Rectangle((PADDING_CARTOONS-200, gene.p_from + PADDING_PLOTS), RECTANGLE_WIDTH, gene.length, color="k")
            _rotation=90
        else:
            
            if gene_name.lower() == 'cas4':
                _fill = True
                _color = 'g'
            elif gene.profile.intersection(crispr_related_profiles):
                _fill = True
                _color = 'b'
            

            if gene.strand == "+":

                x_pos = PADDING_CARTOONS
                y_pos = gene.p_from + PADDING_PLOTS
                dx = 0
                dy = gene.length

            elif gene.strand == "-":

                x_pos = PADDING_CARTOONS
                y_pos = gene.p_to + PADDING_PLOTS
                dx = 0
                dy = -gene.length

            patch = patches.Arrow(
                x_pos,
                y_pos,
                dx,
                dy,
                width=1000,
                fill=_fill,
                color=_color)

        ax.add_patch(patch)
        ax.annotate(gene_name, (PADDING_CARTOONS + 500, gene.p_from + gene.length/2 + PADDING_PLOTS), 
            # rotation=_rotation, 
            va='center', 
            ha='left',
            # size=8,
            weight='bold',
            fontsize=GENE_NAME_FONT_SIZE)

    ax.set_ylabel(locus.label, weight='bold', size=LABEL_FONT_SIZE)


def plot_loci(locus1, locus2, hits):

    fig, ax = plt.subplots(1, figsize=(14, 14))

    for hit in hits:
        ax.plot([hit.q_start, hit.q_end],[hit.s_start, hit.s_end])
        
        # annotate lines with perc_id
        text_x = hit.q_start + int((hit.q_end - hit.q_start)/2)
        text_y = hit.s_start + int((hit.s_end - hit.s_start)/2)
        
        text = "%%%d" % int(round(hit.perc_id))
        plt.annotate(text, xy=(text_x, text_y), rotation=np.rad2deg(np.arctan(hit.slope))-5, fontsize=20)
        
        # shade the projection on X axis
        x_polygon_vertices = [(hit.q_start, PADDING_PLOTS), (hit.q_start, hit.s_start), (hit.q_end, hit.s_end), (hit.q_end, PADDING_PLOTS)]
        x_polygon = patches.Polygon(x_polygon_vertices, closed=True, color='gray', alpha=0.3)
        ax.add_patch(x_polygon)
        
        # shade the projection on Y axis
        y_polygon_vertices = [(PADDING_PLOTS, hit.s_start), (hit.q_start, hit.s_start), (hit.q_end, hit.s_end), (PADDING_PLOTS, hit.s_end)]
        y_polygon = patches.Polygon(y_polygon_vertices, closed=True, color='gray', alpha=0.3)
        ax.add_patch(y_polygon)
        
    plt.axis([0, locus1.length + PADDING_PLOTS, 0, locus2.length + PADDING_PLOTS])
    
    draw_x_cartoons(ax, locus1)
    draw_y_cartoons(ax, locus2)
    
    # The tick labels are not populated until the plot is drawn.
    plt.draw()
    
    # Shift the ticks positions by PADDING_PLOTS
    labels = [label.get_text() for label in ax.get_yticklabels()]
    locations = [ tick+PADDING_PLOTS for tick in ax.get_yticks()]
    ax.set_yticks(locations)
    ax.set_yticklabels(labels)
    labels = [label.get_text() for label in ax.get_xticklabels()]
    locations = [ tick+PADDING_PLOTS for tick in ax.get_xticks()]
    ax.xaxis.set(ticks=locations, ticklabels=labels)
    ax.set_xticks(locations)
    ax.set_xticklabels(labels)
    plt.title("DNA sequence comparison", fontsize=LABEL_FONT_SIZE, weight="bold")
    plt.savefig(image_file)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(usage="%s -q query_locus -s subject_locus -o image_file [-d work_dir]"
                                           % os.path.basename(sys.argv[0]),
                                     description=__description__)
    parser.add_argument("-q", "--query_locus",
                        action="store", dest="query_locus",
                        help="Query locus file")
    parser.add_argument("-s", "--subject_locus",
                        action="store", dest="subject_locus",
                        help="Subject locus file")
    parser.add_argument("-o", "--output",
                        action="store", dest="image_file",
                        help="file to save image to")
    parser.add_argument("-d", "--work_dir",
                        action="store", dest="work_dir", default='/panfs/pan1/patternquest/Projects/NewSystems/data/cas4/in_situ/tmp/',
                        help="Work directory")
    args = parser.parse_args()

    if not args.query_locus or not args.subject_locus or not args.image_file:
        parser.print_help()
        sys.exit()

    locus_file_1 = args.query_locus
    locus_file_2 = args.subject_locus
    work_dir = args.work_dir
    image_file = args.image_file

    locus_dna_1 = os.path.join(work_dir, 'locus1.fna')
    locus_dna_2 = os.path.join(work_dir, 'locus2.fna')
    blast_file = os.path.join(work_dir, 'blast.out')

    locus1 = load_locus(locus_file_1)
    locus2 = load_locus(locus_file_2)

    print("Retrieving nt sequences for loci")
    retreive_locus_nt_sequence(locus1, locus_dna_1)
    retreive_locus_nt_sequence(locus2, locus_dna_2)

    print("Blasting the loci")
    blast_between_loci(locus_dna_1, locus_dna_2, blast_file)

    print("Reading hits and plotting")
    avegare_length = (locus1.length + locus2.length)/float(2)
    hits = read_blast_hits(blast_file, avegare_length)

    plot_loci(locus1, locus2, hits)

