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
import lib.db.generic
from lib.db.generic import pentaplets as p
from lib.db.generic import quadruplets as q
from lib.db.generic import triplets as tr
from lib.db.generic import duplets as d
from itertools import combinations
from itertools import product
import time
import subprocess as sp


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

    for l in lines:
        terms = l.strip().split()
        _annotation = terms[8] if len(terms)==9 else terms[0]
        annotations.append(_annotation)

    annotations = [l.split(',') for l in annotations]

    annotations = flatten(annotations)
    annotations = set(annotations)

    return list(combinations(annotations, combination_size))


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


def write_kmers_to_database(combination_size, prefix):
    cnt = 0
    total = 0
    zeros = 0

    local_cnt = 0

    kplet_pile = []
    extract_start = time.time()

    profile2id = lib.db.generic.map_profile2id(prefix)
    file2id = lib.db.generic.map_file_name2id(prefix)

    if combination_size == 2:
        store_kplets_pile = d.store_kplets_pile
    elif combination_size == 3:
        store_kplets_pile = tr.store_kplets_pile
    elif combination_size == 4:
        store_kplets_pile = q.store_kplets_pile
    elif combination_size == 5:
        store_kplets_pile = p.store_kplets_pile

    for f in os.listdir(neighborhoods_path):

        cnt += 1

        # if cnt < 5500:
        #     continue

        kplets = extract_kplets(os.path.join(neighborhoods_path, f), combination_size)

        if not kplets:
            zeros += 1
            continue

        local_cnt += len(kplets)
        kplet_pile += [(kplets, f)]

        if local_cnt > 10000:

            extract_time = time.time() - extract_start
            insert_start = time.time()
            store_kplets_pile(prefix, kplet_pile, profile2id, file2id)
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

    neighborhoods_path = os.path.join(gv.project_data_path, 'cas1402/files/')

    # count_profiles_in_neighborhoods(neighborhoods_path, './', 500, 10)
    prefix = "cas"
    combination_size = int(sys.argv[1])

    write_kmers_to_database(combination_size, prefix)
