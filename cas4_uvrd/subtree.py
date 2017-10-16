#!/usr/bin/env python
__author__ = 'hudaiber'

import os
import sys
if sys.platform == 'darwin':
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))
elif sys.platform == 'linux2':
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))

import global_variables as gv
sys.path.append(gv.project_code_path)
import subprocess as sp

from Bio import SeqIO


def search_in_prok1603(leaves_file):

    # cmd_list = ['/home/wolf/bin/prok_extract_pty', '-f=7', leaves_file, '>', conf_leaves_file]
    cmd_list = ['/home/wolf/bin/prok_extract_pty', '-f=7', leaves_file]
    print "Total leaves requested:", len([l for l in open(leaves_file)])
    ps = sp.Popen(" ".join(cmd_list), shell=True, stderr=open(os.devnull, 'wb'), stdout=sp.PIPE)

    stdout = ps.communicate()[0]

    stdout = [l for l in stdout.split("\n") if l]

    gis = [l.split()[-1] for l in stdout]
    contigs = [l.split()[4] for l in stdout]

    contigs_file = os.path.join(work_dir, 'contigs.txt')
    with open(contigs_file, 'w') as outf:
        [outf.write("%s\n" % l) for l in contigs]

    print "Confirmed to be in Prok1603:", len(gis)

    cmd_list = ['/home/wolf/bin/prok_extract_pty', '-f=5', contigs_file]
    ps = sp.Popen(" ".join(cmd_list), shell=True, stderr=open(os.devnull, 'wb'), stdout=sp.PIPE)
    stdout = ps.communicate()[0]
    stdout = [l for l in stdout.split("\n") if l]

    rpob_genes = [l.split()[-1] for l in stdout if l.split()[-1] in all_rpob_gis]

    print "RpoB genes found:", len(rpob_genes)

    conf_leaves_file = os.path.join(work_dir, 'conf_leaves.txt')
    with open(conf_leaves_file, 'w') as outf:
        [outf.write("%s\n" % l) for l in gis]

    cas4_seq_file = os.path.join(work_dir, 'sub_cas4.fa')
    cmd_list = ['blastdbcmd',
                '-db=/panfs/pan1/prokdata/db/Prok1603/Prok1603',
                '-entry_batch',
                conf_leaves_file,
                '>',
                cas4_seq_file]
    sp.call(" ".join(cmd_list), shell=True, stderr=open(os.devnull, 'wb'), stdout=sp.PIPE)

    seqs = list(SeqIO.parse(cas4_seq_file, 'fasta'))

    unique_seqs = []
    orgs = []
    for i in range(len(seqs)):
        _seq = seqs[i]
        _org = "_".join(seqs[i].description.split("[")[-1][:-1].split())

        if _org in orgs:
            continue

        _seq.id = "_".join(_seq.description.split("[")[-1][:-1].split()) + "_" + _seq.id.split('|')[1]
        _seq.description = ""

        unique_seqs.append(_seq)
        orgs += [_org]

    SeqIO.write(unique_seqs, open(cas4_seq_file, 'w'), 'fasta')

    cas4_msa_file = os.path.join(work_dir, 'sub_cas4_msa.fa')
    print "Running MSA"
    rpob_msa_file = os.path.join(work_dir, 'rpob_msa.fa')
    cmd_list = ['mafft',
                '--auto',
                cas4_seq_file,
                '>',
                cas4_msa_file]

    sp.call(" ".join(cmd_list), shell=True, stderr=open(os.devnull, 'wb'), stdout=sp.PIPE)

    print "Running FastTree"
    cas4_tree_file = os.path.join(work_dir, 'sub_cas4.tre')
    cmd_list = ['FastTree',
                cas4_msa_file,
                '>',
                cas4_tree_file]
    sp.call(" ".join(cmd_list), shell=True, stderr=open(os.devnull, 'wb'), stdout=sp.PIPE)

    # cas4_tree_file = os.path.join(work_dir, 'confirmed_cas4.tre')
    # cmd_list = ['/home/wolf/bin/add/tree_subtree', global_cas4_tree, '-l=%s' % conf_leaves_file, '>', cas4_tree_file]
    # sp.Popen(" ".join(cmd_list), shell=True, stderr=open(os.devnull, 'wb'), stdout=sp.PIPE)

    rpob_gis_file = os.path.join(work_dir, 'rpob_gis.txt')
    with open(rpob_gis_file, 'w') as outf:
        [outf.write("%s\n" % l) for l in rpob_genes]

    rpob_seq_file = os.path.join(work_dir, 'rpob_orig.fa')
    cmd_list = ['blastdbcmd',
                '-db=/panfs/pan1/prokdata/db/Prok1603/Prok1603',
                '-entry_batch',
                rpob_gis_file,
                '>',
                rpob_seq_file]
    sp.call(" ".join(cmd_list), shell=True, stderr=open(os.devnull, 'wb'), stdout=sp.PIPE)

    seqs = list(SeqIO.parse(rpob_seq_file,'fasta'))

    unique_seqs = []
    orgs = []
    for i in range(len(seqs)):
        _seq = seqs[i]
        _org = "_".join(seqs[i].description.split("[")[-1][:-1].split())

        if _org in orgs:
            continue

        _seq.id = "_".join(_seq.description.split("[")[-1][:-1].split()) + "_" + _seq.id.split('|')[1]
        _seq.description = ""

        unique_seqs.append(_seq)
        orgs += [_org]

    rpob_seq_file = os.path.join(work_dir, 'rpob.fa')
    SeqIO.write(unique_seqs, open(rpob_seq_file, 'w'), 'fasta')

    print "Running MSA"
    rpob_msa_file = os.path.join(work_dir, 'rpob_msa.fa')
    cmd_list = ['mafft',
                '--auto',
                rpob_seq_file,
                '>',
                rpob_msa_file]

    sp.call(" ".join(cmd_list), shell=True, stderr=open(os.devnull, 'wb'), stdout=sp.PIPE)

    print "Running FastTree"
    rpob_tree_file = os.path.join(work_dir, 'rpob.tre')
    cmd_list = ['FastTree',
                rpob_msa_file,
                '>',
                rpob_tree_file]

    # sp.call(" ".join(cmd_list), shell=True, stderr=open(os.devnull, 'wb'), stdout=sp.PIPE)
    sp.call(" ".join(cmd_list), shell=True, stderr=open(os.devnull, 'wb'), stdout=sp.PIPE)

    # print os.system() " ".join(cmd_list)



if __name__=="__main__":

    work_dir = os.path.join(gv.project_data_path, 'cas4/unknowns/cas4_trees/tree_comp/')
    global_cas4_tree = os.path.join(gv.project_data_path, 'cas4/unknowns/cas4_trees/tree_comp/cas4.1.tre')
    leaves_file = os.path.join(work_dir, 'leaves.txt')

    all_rpob_gis = set([l.split(',')[0] for l in open(os.path.join(work_dir, 'COG0085.lst'))])

    search_in_prok1603(leaves_file)


