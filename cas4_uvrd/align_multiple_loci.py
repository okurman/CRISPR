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


def align_locus_sequences(loci):
    print("Retrieving sequences")
    out_files = []
    for locus in loci:
        seq_file=os.path.join(work_dir, "%s.fna" % locus.contig)

        pk.get_nt_sequence('all1603.nt',
                           locus.contig,
                           range=[locus.locus_from, locus.locus_to],
                           strand=locus.strand,
                           out_file=seq_file)
        out_files.append(seq_file)

    seq_file = os.path.join(work_dir, 'all.fna')
    aln_file = os.path.join(work_dir, 'all.msa.fna')

    print("Aligning sequences")
    command = "cat %s > %s" % (" ".join(out_files), seq_file)
    sp.call(command, shell=True)
    command = "muscle -in %s -out %s -maxmb 10000" % (seq_file, aln_file)
    sp.call(command, shell=True)





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

    global work_dir
    work_dir        = os.path.join(data_path, 'cas4/in_situ/align_loci/Porphyromonas_gulae/')
    islands_file    = os.path.join(work_dir, 'loci.txt')
    result_aln_file = os.path.join(work_dir, 'master.msa.fa')
    result_clu_file = os.path.join(work_dir, 'master.msa.clu')

    loci = lh.load_islands(islands_file)

    align_locus_sequences(loci)

