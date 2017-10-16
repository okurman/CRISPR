#!/sw/bin/python2.7
__author__ = 'Sanjarbek Hudaiberdiev'

import os
import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')

import global_variables as gv
import pickle
from lib_old import classes as cl
from lib_old import tools as t
# from lib import db
from itertools import combinations
from itertools import product
import numpy as np
import pandas as pd
from operator import itemgetter


def count_combinations(neighborhoods, profile_list, r, src2org, gnm2weight):

    key_fmt = "-".join(["%s" for i in range(r)])
    combination_count = {}
    for nbr in neighborhoods:
        cur_src = nbr.genes[0].src
        cur_org = src2org[cur_src]
        src_file = os.path.basename(nbr.source_file)

        tmp_list = " ".join([g.cogid for g in nbr.genes])
        tmp_list = [p for p in tmp_list.split() if p in profile_list]
        combination_list = combinations(set(tmp_list), r)

        for comb in combination_list:
            k = key_fmt % comb
            if k in combination_count:
                combination_count[k].count += 1
                combination_count[k].weight += gnm2weight[cur_org]
                combination_count[k].files.update([src_file])
            else:
                combination_count[k] = cl.ProfileCount(1, gnm2weight[cur_org], src_file)

    return combination_count


def count_profiles_in_neighborhoods(neighborhoods_path, save_path, limit_to, combination_size):
    target_profiles = [l.strip() for l in open('/Users/hudaiber/Projects/NewSystems/data/Bacteria/test/test_profiles.txt')]
    src2org = t.map_src2org()
    gnm2weight = t.map_genome2weight()
    neighborhoods = [cl.Neighborhood(os.path.join(neighborhoods_path, f)) for f in os.listdir(neighborhoods_path)]
    
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

    pickle.dump(profile_weights, open('profile_weights.p', 'w'))

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
    M['files'] = [" ".join(v.files) for v in counted_combinations.values()]

    M = M[M['count'] > 1]
    M = M[M['weight'] > pivot_value]
    M = M.sort('weight', ascending=False)
    fname = '%d_%d.tab' % (limit_to, combination_size)
    fname = os.path.join(save_path, fname)
    print fname
    fout = open(fname, 'w')
    M.to_csv(fout, sep="\t", index=False)


if __name__=='__main__':

    combination_size = 5

    neighborhoods_path = os.path.join('/Users/hudaiber/Projects/NewSystems/data/Bacteria/test/genes_and_flanks/pty')
    # neighborhoods_path = os.path.join('/Users/hudaiber/Projects/NewSystems/data/Bacteria/test/genes_and_flanks/merged')
    limit_to = 500
    save_path = '/Users/hudaiber/Projects/NewSystems/data/Bacteria/test/'
    count_profiles_in_neighborhoods(neighborhoods_path, save_path, limit_to, combination_size)

    # combination_size = 6
    # write_kmers_to_database(combination_size, neighborhoods_path)