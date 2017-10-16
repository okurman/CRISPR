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

import lib.utils.tools as t
from lib.db.prok1402.duplets import insert_adjacent_duplets, insert_duplets, insert_source_duplets
import lib.db.prok1402.db_tools as dt
from itertools import combinations, product


def extract_adjacent_duplets_from_file(f, genes=None):

    all_duplets = []

    if not genes:
        genes = t.parse_pty_file(f)
    
    previous_profiles = genes[0].profiles

    longest_multidomain = 0

    if previous_profiles and len(previous_profiles)>1:
        domain_duplets = list(combinations(previous_profiles,2))
        all_duplets.extend([(duplet, f, 1) for duplet in domain_duplets])
        
    for gene in genes[1:]:
        cur_profiles = gene.profiles

        if not cur_profiles:
            previous_profiles = cur_profiles
            continue

        if len(cur_profiles) > 1:
            domain_duplets = list(combinations(cur_profiles, 2))
            all_duplets.extend([(duplet, f, 1) for duplet in domain_duplets])

        if previous_profiles:
            adjacent_duplets = list(product(previous_profiles, cur_profiles))
            all_duplets.extend([(duplet, f, 0) for duplet in adjacent_duplets])

        previous_profiles = cur_profiles

        longest_multidomain = longest_multidomain if longest_multidomain >= len(cur_profiles) else len(cur_profiles)

    # In case the last gene is multidomain
    if previous_profiles and len(previous_profiles)>1:
        domain_duplets = list(combinations(previous_profiles,2))
        all_duplets.extend([(duplet, f, 1) for duplet in domain_duplets])

    return all_duplets


def extract_all_duplets_from_file(f, gi2profiles=None):

    genes = t.parse_pty_file(f)

    if not gi2profiles:
        profiles = [profile for gene in genes for profile in gene.profiles.split()]
    else:
        profiles = [profile for gene in genes for profile in gi2profiles[gene.gid]]

    duplets = list(combinations(profiles,2))
    all_duplets = [(duplet, f) for duplet in duplets]

    return all_duplets


def extract_all_duplets_from_prok1402():

    """
    Extraction adjacent duplets is done by means of recording them in the dictionary pair2weight

    The overall abundance of profiles is also needed. It's recorded in profile2weight
    """
    pty_path = "/panfs/pan1/patternquest/data/Pty/genomes/"
    work_dir = os.path.join(data_path, 'prok1402/graph/graph_files/')

    print("Loading dictionaries")
    gi2profiles = t.map_gi2profiles()
    genome2weight = t.map_genome2weight()
    pair2weight = defaultdict(float)
    pair2count = defaultdict(int)
    profile2weight=defaultdict(float)

    print("Reading Prok1402")
    for root, dirs, files in os.walk(pty_path):
        for f in files:

            if not f.endswith(".pty"):
                continue

            file_name = os.path.join(root, f)
            genome = os.path.basename(root)

            genes = t.parse_pty_file(file_name)
            for gene in genes:
                gene.profiles = gi2profiles[gene.gid]

                for profile in gene.profiles:
                    t.update_dictionary(profile2weight, profile, genome2weight[genome])

            previous_profiles = genes[0].profiles

            if len(previous_profiles) > 1:
                domain_duplets = list(combinations(previous_profiles,2))

                for duplet in domain_duplets:
                    [kplet_1, kplet_2] = sorted(duplet)
                    key = "%s-%s" % (kplet_1, kplet_2)
                    t.update_dictionary(pair2weight, key, genome2weight[genome])
                    t.update_dictionary(pair2count, key, 1)

            for gene in genes[1:]:
                    cur_profiles = gene.profiles

                    if not previous_profiles:
                        previous_profiles = cur_profiles
                        continue

                    if len(cur_profiles) > 1:
                        domain_duplets = list(combinations(previous_profiles, 2))

                        for duplet in domain_duplets:
                            [kplet_1, kplet_2] = sorted(duplet)
                            key = "%s-%s" % (kplet_1, kplet_2)
                            t.update_dictionary(pair2weight, key, genome2weight[genome])
                            t.update_dictionary(pair2count, key, 1)

                    adjacent_duplets = list(product(previous_profiles, cur_profiles))

                    for duplet in adjacent_duplets:
                        [kplet_1, kplet_2] = sorted(duplet)
                        key = "%s-%s" % (kplet_1, kplet_2)
                        t.update_dictionary(pair2weight, key, genome2weight[genome])
                        t.update_dictionary(pair2count, key, 1)

                    previous_profiles = cur_profiles

    print("Writing to files")
    with open(os.path.join(work_dir, "prok1402_adj_duplets_weights.txt"), "w") as outf:

        for (key,weight) in sorted(pair2weight.items(), key=lambda x: x[1], reverse=True):
            [kplet_1, kplet_2] = key.split("-")
            outf.write("%s\t%s\t%f\n" % (kplet_1, kplet_2, weight))

    with open(os.path.join(work_dir, "prok1402_profile_abundance.txt"), "w") as outf:
        for (profile,weight) in sorted(profile2weight.items(), key=lambda x: x[1], reverse=True):
            outf.write("%s\t%f\n" % (profile, weight))


def insert_all_adj_duplets_from_prok1402():

    pty_path = "/panfs/pan1/patternquest/data/Pty/genomes/"

    cnt = 0
    dirs = os.listdir(pty_path)

    for dir in dirs:

        print cnt+1, "/", len(dirs)
        print dir
        for f in os.listdir(os.path.join(pty_path, dir)):

            file_name = os.path.join(pty_path, dir, f)
            source = f.split(".")[0]
            source_id = source2id[source]

            genes = t.parse_pty_file(file_name)
            for gene in genes:
                gene.profiles = gi2profiles[gene.gid]

            all_duplets = extract_adjacent_duplets_from_file(source, genes)
            
            if all_duplets:
                insert_source_duplets(all_duplets, profile2id, source_id)
            
        cnt += 1


def insert_baited_adj_duplets(files):

    cnt = 1
    for f in files:
        print cnt, "/", len(files), f
        # f = "NC_008346_14.pty"
        file_path = os.path.join(files_dir, f)
        
        # Adjacent duplets
        duplets = extract_adjacent_duplets_from_file(file_path)
        
        if duplets:
            insert_adjacent_duplets(duplets, profile2id, file2id)
        cnt += 1


def insert_baited_duplets():

    cnt = 1
    for f in files:
        print cnt, "/", len(files), f

        file_path = os.path.join(files_dir, f)

        duplets = extract_all_duplets_from_file(file_path)
        if duplets:
            insert_duplets(duplets, profile2id, file2id)

        cnt += 1


def test_db_source():

    inf = "/panfs/pan1/patternquest/data/Pty/genomes/Acaryochloris_marina_MBIC11017_uid58167/NC_009925.pty"
    genes = t.parse_pty_file(inf)
    for gene in genes:
        gene.profiles = gi2profiles[gene.gid]
    outf = "db_test/NC_009925.pty"
    t.write_pty_file(genes, outf)


if __name__ == "__main__":

    extract_all_duplets_from_prok1402()
    # print "Loading ccp"
    # gi2profiles = t.map_gi2profiles()
    # print "Loading pty"
    # genes = t.parse_pty_file("/panfs/pan1/patternquest/data/Pty/genomes/Achromobacter_xylosoxidans_A8_uid59899/NC_014640.pty")
    # ccp_file = "/panfs/pan1/patternquest/data/Pty/genomes/Achromobacter_xylosoxidans_A8_uid59899/NC_014640.ccp"
    # for gene in genes:
    #     gene.profiles = ",".join(gi2profiles[gene.gid])
    # print "Writing pty"
    # t.write_genes_to_pty(genes,ccp_file)

    sys.exit()
    work_dir = os.path.join(data_path, 'prok1402/genes_and_flanks/win_10/')
    files_dir = os.path.join(work_dir, 'merged')

    profile2id = dt.map_profile2id("all")
    source2id = dt.map_source2id()
    file2id, _ = dt.map_baited_file2id()
    files = os.listdir(files_dir)
    gi2profiles = t.map_gi2profiles()

    # insert_baited_adj_duplets(files)
    # insert_baited_duplets()
    # insert_all_adj_duplets_from_prok1402()
    