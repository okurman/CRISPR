#!/home/hudaiber/env2.7/bin/python

import os
import sys
import configparser
from collections import defaultdict

###############################################################
config_file = os.path.join(os.path.expanduser('~'), 'paths.cfg')
cfg = configparser.ConfigParser()
cfg.read(config_file)

code_path = cfg.get('NewSystems', 'code_path')
data_path = cfg.get('NewSystems', 'data_path')
sys.path.append(code_path)
###############################################################

import lib.utils.tools as t
from lib.db.prok1402.duplets import insert_adjacent_duplets, insert_duplets, insert_source_duplets
import lib.db.prok1402.db_tools as dt
from itertools import combinations, product


def check_case_of_kplet():

    """
    Take this query:
        select *
        from prok1402_all_adj_duplet_source pds
        inner join sources s on pds.source_id=s.id
        where kplet_id=44061 and s.id=79;
    This will result in 32 rows. Let's check from the file, if it is correct.
    The source with id 79 is NC_014640.
    The kplet_id = 44061 is: pfam00126 pfam03401

    """

    print "Loading ccp"
    gi2profiles = t.map_gi2profiles()
    print "Loading pty"
    genes = t.parse_pty_file(
        "/panfs/pan1/patternquest/data/Pty/genomes/Achromobacter_xylosoxidans_A8_uid59899/NC_014640.pty")
    ccp_file = "/panfs/pan1/patternquest/data/Pty/genomes/Achromobacter_xylosoxidans_A8_uid59899/NC_014640.ccp"
    for gene in genes:
        gene.profiles = ",".join(gi2profiles[gene.gid])
    print "Writing pty"
    t.write_genes_to_pty(genes, ccp_file)

    """
    After this.
    "grep -B 1 -A 1 pfam00126 NC_014640.ccp | grep -c pfam03401"
    gives the expected 32. 
    """


def check_case_of_a_multidomain_kplet():

    """
    Take this query:
        select *
        from prok1402_all_adj_duplet_source pds
        inner join sources s on pds.source_id=s.id
        where kplet_id=68 and s.id=1;
    Gives 28 rows.
    The source with id 1 is NC_009925.
    The kplet_id = 68 is a multidomain: COG2205 cd00156
    """

    print "Loading ccp"
    gi2profiles = t.map_gi2profiles()
    print "Loading pty"
    genes = t.parse_pty_file(
        "/panfs/pan1/patternquest/data/Pty/genomes/Acaryochloris_marina_MBIC11017_uid58167/NC_009925.pty")
    ccp_file = "/panfs/pan1/patternquest/data/Pty/genomes/Acaryochloris_marina_MBIC11017_uid58167/NC_009925.ccp"
    for gene in genes:
        gene.profiles = ",".join(gi2profiles[gene.gid])
    print "Writing pty"
    t.write_genes_to_pty(genes, ccp_file)

    """
    Command:
     "grep -Ec "cd00156.*COG2205" NC_009925.ccp" gives 28
    """


if __name__ == "__main__":


    check_case_of_a_multidomain_kplet()
