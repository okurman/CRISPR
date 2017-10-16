#!/home/hudaiber/env2.7/bin/python

__description__ = """ Align loci segmenwise with annotations """

import os
import sys

import configparser
###############################################################
config_file = os.path.join(os.path.expanduser('~'), 'paths.cfg')
cfg = configparser.ConfigParser()
cfg.read(config_file)
code_path = cfg['NewSystems']['code_path']
data_path = cfg['NewSystems']['data_path']
sys.path.append(code_path)
###############################################################

import subprocess as sp
from lib.utils import blaster
from lib.utils import locus_helper as lh
import lib.utils.prok1603.parser as pk

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import argparse

# threshold level identity for blast hits to be considered
THR_ID = 70
SHORT_HIT_LENGTH = 500



def process_segment(loci, i):

    segment_seqs_file = os.path.join(work_dir, "segment_%s.fa" % i)
    segment_msa_file = os.path.join(work_dir, "segment_%s.msa.fa" % i)
    
    segment_temp_files = []

    for locus in loci:
        genes = sorted([gene for gene in locus.genes if gene.order==i])

        if not genes:
            continue

        _segment_file = pk.get_nt_sequence('all1603.nt',
                                           locus.contig,
                                           range=[genes[0].p_from, genes[-1].p_to],
                                           strand=locus.strand)
        segment_temp_files.append(_segment_file)

    if os.path.exists(segment_seqs_file):
        os.unlink(segment_seqs_file)

    for _file in segment_temp_files:
        sp.call("cat %s >> %s" % (_file, segment_seqs_file), shell=True)
        os.unlink(_file)

    cmd_muscle = "muscle -in %s -out %s -maxmb 10000" % (segment_seqs_file, segment_msa_file)

    sp.call(cmd_muscle, shell=True, stderr=open(os.devnull, 'wb'))

    return segment_msa_file


def process_inter_segment(loci, i, j):

    segment_seqs_file = os.path.join(work_dir, "inter_segment_%s-%s.fa" % (i, j))
    segment_msa_file = os.path.join(work_dir, "inter_segment_%s-%s.msa.fa" % (i, j))
    
    segment_temp_files = []

    for locus in loci:

        prior_segment = sorted([gene for gene in locus.genes if gene.order==i])
        if not prior_segment:
            continue

        post_segment = sorted([gene for gene in locus.genes if gene.order==j])
        if not post_segment:
            continue

        _segment_file = pk.get_nt_sequence('all1603.nt',
                                           locus.contig,
                                           range=[prior_segment[-1].p_to, post_segment[0].p_from],
                                           strand=locus.strand)
        segment_temp_files.append(_segment_file)

    if os.path.exists(segment_seqs_file):
        os.unlink(segment_seqs_file)

    for _file in segment_temp_files:
        sp.call("cat %s >> %s" % (_file, segment_seqs_file), shell=True)
        os.unlink(_file)

    cmd_muscle = "muscle -in %s -out %s -maxmb 10000" % (segment_seqs_file, segment_msa_file)
    sp.call(cmd_muscle, shell=True, stderr=open(os.devnull, 'wb'))

    return segment_msa_file


def align_segment_sequences(loci):

    orders = set()
    [orders.update(locus.orders) for locus in loci]
    orders = sorted(orders, key=lambda x: int(x))

    for ix, order in enumerate(orders):
        if ix > 0:
            prev_order = orders[ix-1]
            print("Aligning intersegment: %s-%s" % (prev_order, order))
            _msa_file = process_inter_segment(loci, prev_order, order)

        print("Aligning segment: %d/%d (order=%s)" % (ix,len(orders), order))
        _msa_file = process_segment(loci, order)


def join_to_master_alignment(big_alignment, subject_alignment, segment_type, segment_name=""):

    s_aln_len = len(list(subject_alignment.values())[0].seq)
    if segment_type == 'segment':
        middle_space=" "*(s_aln_len - 6 - 2*len(segment_name))
        desc_line = "==>{0}{1}{2}<==".format(segment_name,
                                             middle_space,
                                             segment_name)
    else:
        #  if it is an overlap
        if s_aln_len == 4 and list(subject_alignment.values())[0].seq == "____":
            desc_line = "ovrl"
        else:
            desc_line = " " * s_aln_len
    big_alignment[0].seq += desc_line

    for i in range(1, len(big_alignment)):
        _contig = big_alignment[i].id

        if _contig in subject_alignment:
            _seq = subject_alignment[_contig].seq
        else:
            _seq = "-" * s_aln_len

        big_alignment[i].seq += _seq


def merge_alignments(loci):

    orders = set()
    [orders.update(locus.orders) for locus in loci]
    orders = sorted(orders, key=lambda x: int(x))

    #crete an empty fasta record list
    master_alignment = [SeqRecord(seq=Seq(""), id="description", description="")]
    master_alignment += [SeqRecord(seq=Seq(""), id=locus.contig, description="") for locus in loci]

    prev_order = None
    for order in orders:

        if prev_order:
            aln_file = os.path.join(work_dir, "inter_segment_%s-%s.msa.fa" % (prev_order, order))
            if os.path.exists(aln_file):
                aln = SeqIO.to_dict(SeqIO.parse(aln_file,"fasta"))
            else:
                aln = {locus.contig:SeqRecord(seq=Seq("____"), id=locus.contig) for locus in loci}
            assert len(set([len(seq) for seq in aln.values()])) == 1

            join_to_master_alignment(master_alignment, aln, 'intersegment')
            # print(prev_order,order,len(master_alignment[0].seq), len(master_alignment[1].seq))

        aln_file = os.path.join(work_dir, "segment_%s.msa.fa" % order)
        aln = SeqIO.to_dict(SeqIO.parse(aln_file, "fasta"))
        assert len(set([len(seq) for seq in aln.values()])) == 1

        genes = []
        [genes.extend(locus.genes) for locus in loci]
        segment_name = "_".join(list(set(gene.name.split(",")[0] for gene in genes if gene.order==order)))

        join_to_master_alignment(master_alignment, aln, 'segment', segment_name)
        # print(order, len(master_alignment[0].seq), len(master_alignment[1].seq))

        prev_order = order

    acc2cas_clu = {locus.contig:locus.cas8_cluster for locus in loci}

    for seq in master_alignment[1:]:
        seq.id="%s_%s" % (seq.id, acc2cas_clu[seq.id])

    return master_alignment


if __name__ == '__main__':

    # parser = argparse.ArgumentParser(usage="%s -q query_locus -s subject_locus -o image_file [-d work_dir]"
    #                                        % os.path.basename(sys.argv[0]),
    #                                  description=__description__)
    # parser.add_argument("-q", "--query_locus",
    #                     action="store", dest="query_locus",
    #                     help="Query locus file")
    # parser.add_argument("-s", "--subject_locus",
    #                     action="store", dest="subject_locus",
    #                     help="Subject locus file")
    # parser.add_argument("-o", "--output",
    #                     action="store", dest="image_file",
    #                     help="file to save image to")
    # parser.add_argument("-d", "--work_dir",
    #                     action="store", dest="work_dir",
    #                     default='/panfs/pan1/patternquest/Projects/NewSystems/data/cas4/in_situ/tmp/',
    #                     help="Work directory")
    # args = parser.parse_args()
    #
    # if not args.query_locus or not args.subject_locus or not args.image_file:
    #     parser.print_help()
    #     sys.exit()

    global workdir
    work_dir = os.path.join(data_path, 'cas4/in_situ/align_loci/MethanoSarc/')
    islands_file = os.path.join(work_dir, 'loci.txt')
    result_aln_file = os.path.join(work_dir, 'master.msa.fa')
    result_clu_file = os.path.join(work_dir, 'master.msa.clu')

    loci = lh.load_islands(islands_file)

    align_segment_sequences(loci)
    # sys.exit()
    master_alignment = merge_alignments(loci)
    # for seq in master_alignment:
    #     print(len(seq))

    SeqIO.write(master_alignment,open(result_aln_file,"w"),"fasta")
    SeqIO.write(master_alignment, open(result_clu_file, "w"), "clustal")

