#!/home/hudaiber/env2.7/bin/python
__author__ = 'hudaiber'

import os
import sys
import glob
import configparser
import tempfile
import subprocess as sp
from collections import defaultdict, Counter


###############################################################
config_file = os.path.join(os.path.expanduser('~'),'paths.cfg')
cfg=configparser.ConfigParser()
cfg.read(config_file)
code_path = cfg['BioPy']['code_path']
sys.path.append(code_path)
code_path = cfg['NewSystems']['code_path']
sys.path.append(code_path)
###############################################################

prok1603_path = os.path.join(os.path.expanduser('~'),'data/Prok1603/')
pty_path = os.path.join(prok1603_path, 'pty')

prok_db_path = cfg['prok1603']['db_path']
import lib.utils.tools as t
from BioClasses import Gene

annotation_file="/panfs/pan1/prokdata/db_tmp/Prok1603/Annotation/Prok1603.ccp.csv"
crispr_annotation_file="/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/CRISPR_info_known_091216_polarised_strict.tsv"
cluster_annotation_file = "/panfs/pan1/patternquest/data/Prok1603/unannotated/cluster_assignments.txt"

devnull = open(os.devnull, "w")


def map_gi2profiles():

    gi2profiles=defaultdict(list)

    for l in open(annotation_file):
        parts = l.split(",")
        gi = parts[0]
        profile = parts[6]
        gi2profiles[gi].append(profile)

    return gi2profiles


def map_gi2cluster():

    gi2cluster = {}

    for l in open(cluster_annotation_file):
        parts = l.strip().split("\t")
        cluster = parts[0]

        for gi in parts[1].split(" "):
            gi2cluster[gi] = cluster

    return gi2cluster


def get_pty_organism(organism):

    source_files = glob.glob(os.path.join(pty_path,"%s/*.pty2"%organism))

    sources=defaultdict(dict)

    for _f in source_files:
        _source = os.path.splitext(os.path.basename(_f))[0]
        _genes = t.parse_pty_file(_f)
        sources[_source] = _genes

    return sources


def get_crispr_annotations():

    crispr_genes=defaultdict(list)

    for l in open(crispr_annotation_file):

        parts = l.split()

        pFrom = parts[2]
        pTo = parts[3]
        source = parts[1]
        organism = parts[10]

        _gene = Gene(source, pFrom, pTo, gi="crispr", organism=organism, strand="+",
                     profiles="crispr", gene_names="crispr")

        crispr_genes[source].append(_gene)

    return crispr_genes


def map_org2weight():

    _fname = os.path.join(prok1603_path, 'Prok1603_weights.tab')

    org2weight = defaultdict(float)

    for l in open(_fname):
        if l.startswith("#"):
            continue
        _org, _weight = l.split("\t")
        org2weight[_org] = float(_weight.strip())

    return org2weight


def map_profile2weight():

    _fname = os.path.join(prok1603_path, 'profile_weight.txt')

    profile2weight = defaultdict(float)

    for l in open(_fname):
        if l.startswith("#"):
            continue
        _pr, _weight = l.split("\t")
        profile2weight[_pr] = float(_weight.strip())

    return profile2weight


def get_db_path(db_name):
    """TODO: define path for other databases as well, as the need arises"""
    _db_path = ""
    if db_name == 'all1603':
        _db_path = os.path.join(prok_db_path,'all1603')
    elif db_name == 'all1603.nt':
        _db_path = os.path.join(prok_db_path, 'all1603.nt')

    if not _db_path:
        raise Exception("Database path could not be determined.")

    return _db_path


def get_nt_sequence(db, acc, range=None, strand=None, out_file=None):

    db_path = get_db_path(db)
    command = "blastdbcmd -db {path} -entry {acc}".format(path=db_path, acc=acc)

    if range:
        command += " -range {start}-{end}".format(start=range[0], end=range[1])

    if strand == "minus":
        command += " -strand {strand}".format(strand=strand)

    if not out_file:
        tempf = tempfile.NamedTemporaryFile(prefix="%s-%s" % (db, acc))
        tempf.close()
        out_file = tempf.name

    command += " > {file_name}".format(file_name=out_file)
    sp.call(command, shell=True, stderr=open(os.devnull, 'wb'))

    return out_file


def add_annotation_to_pty():

    pty_dir = os.path.join(prok1603_path, "pty")

    print "Loading dictionaries"
    gi2profiles = map_gi2profiles()
    gi2cluster = map_gi2cluster()
    profile2gene = t.map_cdd_defense2gene_name()

    print "Maps loaded. Starting annotations"
    cnt = 1
    for d in os.listdir(pty_dir):

        # if cnt < 1937:
        #     cnt += 1
        #     continue

        print cnt, d

        for f in glob.glob("%s/*.pty" % (os.path.join(pty_dir, d))):
            genes = t.parse_plain_pty_file(f)

            if not genes:
                continue

            for gene in genes:
                if gene.gid in gi2profiles:
                    gene.profiles = gi2profiles[gene.gid]
                    gene.gene_names = " ".join([profile2gene[p] for p in gene.profiles])
                elif gene.gid in gi2cluster:
                    gene.profiles = gi2cluster[gene.gid]

            new_f = f + "2"

            t.write_genes_to_pty(genes, new_f)

        cnt += 1

if __name__=="__main__":



    # print gi2profiles["1001253026"]

    # crisprs = get_crispr_annotations()
    # print crisprs
    add_annotation_to_pty()