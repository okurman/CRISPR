#!/sw/bin/python2.7
__author__ = 'Sanjarbek Hudaiberdiev'

import os
import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')
sys.path.append('../../')

import global_variables as gv
from lib import classes as cl
from lib.utils import tools as t
import lib.db
from lib.db.bacteria import db_tools
from itertools import combinations
from itertools import product
import numpy as np
import pandas as pd
from operator import itemgetter
import time


def flatten(input_list):
    list_is_nested = True

    while list_is_nested:
        keep_checking = False
        temp = []

        for element in input_list:
            if isinstance(element,list) or isinstance(element,tuple):
                temp.extend(element)
                keep_checking = True
            else:
                temp.append(element)

        list_is_nested = keep_checking
        input_list = temp[:]

    return input_list


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


def count_profiles_in_neighborhoods(neighborhoods_path, save_path, limit_to, combination_size):
    target_profiles = [l.strip() for l in open(os.path.join(gv.project_data_path, 'Archea', 'arCOG/selected_arcogs.txt'))]
    src2org = t.map_src2org()
    gnm2weight = t.map_genome2weight()
    neighborhoods = [cl.Neighborhood(os.path.join(neighborhoods_path, f)) for f in os.listdir(neighborhoods_path)]

    # pickle.dump(neighborhoods, open('files/neighborhoods.p', 'w'))
    # neighborhoods = pickle.load(open('files/neighborhoods.p'))
    
    profile_stats = {}
    for nbr in neighborhoods:
        src_name = nbr.genes[0].src
        org_name = src2org[src_name]
        for g in nbr.genes:
            if g.cogid == "":
                continue
            for tmpCog in g.cogid.split():
                if tmpCog in target_profiles:
                    continue
                if tmpCog in profile_stats:
                    profile_stats[tmpCog].weight += gnm2weight[org_name]
                    profile_stats[tmpCog].count += 1
                else:
                    profile_stats[tmpCog] = cl.ProfileCount(1, gnm2weight[org_name])

    profile_weights = [(k, v.weight) for k, v in profile_stats.items()]
    profile_weights = sorted(profile_weights, key=itemgetter(1), reverse=True)

    # pickle.dump(profile_weights, open('files/profile_weights.p', 'w'))
    # profile_weights = pickle.load(open('files/profile_weights.p'))

    top_profiles = [k for (k, v) in profile_weights[:limit_to]]
    print 'started counting'
    counted_combinations = count_combinations(neighborhoods, top_profiles, combination_size, src2org, gnm2weight)
    print 'Done counting'
    weight_values = np.array([v.weight for v in counted_combinations.values()])
    weight_values.sort()
    weight_values = weight_values[::-1]
    pivot_ind = np.where(np.cumsum(weight_values)/np.sum(weight_values)>=0.9)[0][0]
    pivot_value = weight_values[pivot_ind]

    M = pd.DataFrame([], columns=['Comb', 'weight', 'count'])
    M['Comb'] = counted_combinations.keys()
    M['weight'] = [v.weight for v in counted_combinations.values()]
    M['count'] = [v.count for v in counted_combinations.values()]

    M = M[M['count'] > 1]
    M = M[M['weight'] > pivot_value]
    M = M.sort('weight',ascending=False)
    fname = '%d_%d.tab' % (limit_to, combination_size)
    fname = os.path.join(save_path, fname)
    print fname
    fout = open(fname, 'w')
    M.to_csv(fout, sep="\t", index=False)


def extract_kplets(file, combination_size):

    annotations = [l.strip().split('\t')[5] for l in open(file).readlines() if len(l.split('\t')) > 5]
    annotations = [l.split() for l in annotations]

    singles = [a for l in annotations if len(l)==1 for a in l]
    singles = list(set(singles))

    multiples = [l for l in annotations if len(l) > 1]

    if 0 < len(multiples) < combination_size:
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


def write_kmers_to_database(combination_size):
    cnt = 0
    total = 0
    zeros = 0

    local_cnt = 0

    kplet_pile = []
    extract_start = time.time()

    cdd2id = lib.db.map_cdd2id()
    file2id = db_tools.map_file_name2id()

    for f in os.listdir(neighborhoods_path):

        cnt += 1

        # if cnt < 5796:
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
            db.store_kplets_pile(kplet_pile, cdd2id, file2id)

            insert_time = time.time() - insert_start
            print 'File count:', cnt, 'Block size:', local_cnt, 'Extract time:', extract_time, extract_time/local_cnt, \
                  'Insert time:', insert_time, insert_time/local_cnt

            kplet_pile = []
            local_cnt = 0
            extract_start = time.time()

        total += len(kplets)

    print 'Nbrhds:', cnt, 'total combinations', total, 'zeros', zeros



if __name__ == '__main__':

    combination_size = int(sys.argv[1])
    print 'Combination size', combination_size
    neighborhoods_path = os.path.join(gv.project_data_path, 'Bacteria', 'genes_and_flanks', 'win_10', 'raw_nbr_files')

    if combination_size == 3:
        from lib.db.bacteria import triplets as db
        write_kmers_to_database(combination_size)

    elif combination_size == 2:
        from lib.db.bacteria import duplets as db
        write_kmers_to_database(combination_size)
