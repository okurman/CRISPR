#!/usr/bin/env python
__author__ = 'Sanjarbek Hudaiberdiev'

import os
import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')
sys.path.append('../')

import global_variables as gv
# from lib import classes as cl
import lib.db.crispr
from lib.db.crispr import pentaplets as p
from lib.db.crispr import quadruplets as q
from lib.db.crispr import triplets as tr
from lib.db.crispr import duplets as d
from itertools import combinations
from itertools import product
import time


def flatten(the_list):
    is_nested = True

    while is_nested:
        keep_checking = False
        temp = []

        for element in the_list:
            if isinstance(element, list) or isinstance(element, tuple):
                temp.extend(element)
                keep_checking = True
            else:
                temp.append(element)

        is_nested = keep_checking
        the_list = temp[:]

    return the_list


def count_combinations(neighborhoods, profile_list, r, src2org, gnm2weight):

    key_fmt = "-".join(["%s" for i in range(r)])
    combination_count = {}
    for n in neighborhoods:
        cur_src = n.genes[0].src
        cur_org = src2org[cur_src]

        tmp_list = " ".join([g.cogid for g in n.genes])
        tmp_list = [p for p in tmp_list.split() if p in profile_list]
        combination_list = combinations(set(tmp_list), r)

        for comb in combination_list:
            k = key_fmt % comb
            if k in combination_count:
                combination_count[k].count += 1
                combination_count[k].weight += gnm2weight[cur_org]
            else:
                combination_count[k] = cl.ProfileCount(1, gnm2weight[cur_org])
    return combination_count


def extract_kplets(file, combination_size):

    annotations = list()
    lines = open(file).readlines()

    if len(lines) > 20:
        for i in range(len(lines)):
            if 'seed' in lines[i] or 'Seed' in lines[i]:
                lines = lines[i-10 : i+10]
                break

    for l in lines:
        terms = l.strip().split('\t')
        if terms[8] == "Unidentified":
            if len(terms) >= 17:
                _annotation = terms[16].strip()
            else:
                continue
        elif terms[8] == "CRISP":
            _annotation = "CRISPR"
        else:
            _annotation = terms[6]

        annotations.append(_annotation)

    annotations = [l.split(',') for l in annotations]

    singles = [a for l in annotations if len(l) == 1 for a in l]
    singles = list(set(singles))

    multiples = [l for l in annotations if len(l) > 1]

    if combination_size > len(multiples) > 0:
        tmp_comb_size = combination_size - len(multiples)
        single_combinations = combinations(singles, tmp_comb_size)

        tmp_comb = []
        for l in multiples:
            if not tmp_comb:
                tmp_comb = l
                continue
            tmp_comb = product(tmp_comb, l)
            tmp_comb = list(tmp_comb)

        flattened_tmp_comb = []

        for comb in tmp_comb:
            flattened_tmp_comb.append([comb] if isinstance(comb, str) else flatten(comb))
        joined_combinations = []
        for single_comb in single_combinations:
            for multi_comb in flattened_tmp_comb:
                joined_combinations.append(list(single_comb)+multi_comb)
        return joined_combinations

    if not multiples:
        return list(combinations(singles, combination_size))

    if len(multiples) >= combination_size:
        # TODO: re-implement this dirty hack properly.
        singles += flatten(multiples)
        singles = list(set(singles))
        return list(combinations(singles, combination_size))


def write_kmers_to_database_old(combination_size, neighborhoods_path):
    cnt = 0
    total = 0
    zeros = 0

    for f in os.listdir(neighborhoods_path):

        cnt += 1

        kplets = extract_kplets(os.path.join(neighborhoods_path, f), combination_size)

        if combination_size == 5:
            p.store_kplets(kplets, f)
        elif combination_size == 4:
            q.store_kplets(kplets, f)
        elif combination_size == 3:
            tr.store_kplets(kplets, f)
        elif combination_size == 2:
            d.store_kplets(kplets, f)

        total += len(kplets)
        if len(kplets) == 0:
            zeros += 1
        print cnt, f, len(kplets)
    print 'Nbrhds:', cnt, 'total combinations', total, 'zeros', zeros


def write_kmers_to_database(combination_size, dataset):
    cnt = 0
    total = 0
    zeros = 0

    local_cnt = 0

    kplet_pile = []
    extract_start = time.time()

    profile2id = lib.db.crispr.map_profile2id(dataset=dataset)
    file2id = lib.db.crispr.map_file_name2id(dataset=dataset)

    if combination_size == 2:
        from lib.db.crispr.duplets import store_kplets_pile
    elif combination_size == 3:
        from lib.db.crispr.triplets import store_kplets_pile
    elif combination_size == 4:
        from lib.db.crispr.quadruplets import store_kplets_pile
    elif combination_size == 5:
        from lib.db.crispr.pentaplets import store_kplets_pile

    # Temporary! Figure out how to deal with it!
    # to_read_files = []
    # [to_read_files.append(l.strip()) for l in open(os.path.join(gv.project_data_path, 'CRISPR/redundancy_elimination/cas1_1_singleton_loci.txt'))]
    # [to_read_files.append(l.strip().split()[-1]) for l in open(os.path.join(gv.project_data_path, 'CRISPR/redundancy_elimination/cas1_2_collapsed_loci_rep.txt'))]

    for f in os.listdir(neighborhoods_path):

        # if f not in to_read_files:
        #     continue

        cnt += 1

        if cnt < 5500:
            continue

        kplets = extract_kplets(os.path.join(neighborhoods_path, f), combination_size)

        if not kplets:
            zeros += 1
            continue

        local_cnt += len(kplets)
        kplet_pile += [(kplets, f)]
        
        if local_cnt > 10000:

            extract_time = time.time() - extract_start
            insert_start = time.time()
            store_kplets_pile(kplet_pile, profile2id, file2id)
            insert_time = time.time() - insert_start
            print 'File count:', cnt, 'Block size:', local_cnt, 'Extract time:', extract_time, extract_time/local_cnt, \
                  'Insert time:', insert_time, insert_time/local_cnt

            kplet_pile = []
            local_cnt = 0
            extract_start = time.time()

        total += len(kplets)
    print cnt
    print 'Nbrhds:', cnt, 'total combinations', total, 'zeros', zeros


if __name__=='__main__':

    neighborhoods_path = os.path.join(gv.project_data_path, 'CRISPR/datasets/crispr/wgs/')
    # count_profiles_in_neighborhoods(neighborhoods_path, './', 500, 10)
    dataset = 3
    combination_size = int(sys.argv[1])
    write_kmers_to_database(combination_size, dataset)
    # combination_size = 3
    # write_kmers_to_database(combination_size, neighborhoods_path)
