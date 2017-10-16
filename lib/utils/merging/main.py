__author__ = 'Sanjarbek Hudaiberdiev'

import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')
import global_variables as gv

from lib.utils import tools as t
from lib.db.archea import neighborhoods_path
from lib.utils.classes import NeighborhoodFileSummary, Neighborhood, KpletList
import os

# Globals
_gid2arcog_cdd = t.map_gid2arcog_cdd()
# _neighborhoods_path = neighborhoods_path()


def get_file_summaries(profile_files, neighborhood_files_path, org2weight):

    file_summaries = list()

    for f in profile_files:

        _neighborhood = Neighborhood(os.path.join(neighborhood_files_path, f))
        _src = _neighborhood.genes[0].src
        _org = _neighborhood.genes[0].organism
        _weight = org2weight[_org]
        source_file = os.path.join(gv.pty_data_path, _org, "%s.pty" % _src)

        if os.path.exists(source_file):
            _neighborhood.extend_flanks(10, source_file , _gid2arcog_cdd)

        file_summaries.append(NeighborhoodFileSummary(f, [], _neighborhood, _org, _src, _weight))

    file_summaries.sort(key= lambda x: x.weight, reverse=True)
    return file_summaries


def merge_into_file_summaries(kplets, neighborhood_files_path, file2src_src2org_map, data_type='bacteria'):

    _org2weight = t.map_genome2weight()

    _file2kplets = dict()
    for kplet in kplets:
        for f in kplet.files:
            if f in _file2kplets:
                _file2kplets[f].append(kplet)
            else:
                _file2kplets[f] = [kplet]

    kplet_files = _file2kplets.keys()
    _file2src, _src2org = file2src_src2org_map(kplet_files)

    file_summaries = list()
    for f in kplet_files:
        _neighborhood = Neighborhood(os.path.join(neighborhood_files_path, f))
        _src = _file2src[f]
        _org = _src2org[_src]
        _weight = _org2weight[_org]
        kplets = _file2kplets[f]
        _neighborhood.extend_flanks(10, os.path.join(gv.pty_data_path, _org, "%s.pty" % _src), _gid2arcog_cdd)
        file_summaries.append(NeighborhoodFileSummary(f, kplets, _neighborhood, _org, _src, _weight))

    # file_summaries = trim_file_summary_list(file_summaries, data_type)
    # file_summaries = [fs for fs in file_summaries if fs]

    # Updating the map _file2src after trimming.
    # new_file_list = [ fs.file_name for fs in file_summaries]
    # for _file_name in _file2src.keys():
    #     if _file_name not in new_file_list:
    #         del _file2src[_file_name]

    # if len(file_summaries) < 2:
    #     return None, None, None, None, None, None

    file_summaries.sort(key= lambda x: x.weight, reverse=True)

    community_count_with_flanks = {}
    community_count = {}
    _org2weight = t.map_genome2weight()

    total_weight = 0

    for i in range(len(file_summaries)):
        cur_file_summary = file_summaries[i]
        _weight = _org2weight[cur_file_summary.org]
        total_weight += _weight
        for gene in cur_file_summary.neighborhood.genes:
            if gene.tag == 'flank':
                for k in gene.cogid.split():
                    t.update_dictionary(community_count_with_flanks, k, _weight)
            else:
                for k in gene.cogid.split():
                    t.update_dictionary(community_count_with_flanks, k, _weight)
                    t.update_dictionary(community_count, k, _weight)
    community = []
    return _src2org, file_summaries, community, community_count, community_count_with_flanks, total_weight


def merge_similar_files(files):

    neighborhoods = t.load_neighborhoods(_neighborhoods_path, files)

    merged_out = [0 for i in range(len(neighborhoods))]

    for i in range(len(neighborhoods)):
        if merged_out[i] == 1:
            continue

        pivot_gids = set([g.gid for g in neighborhoods[i].genes])

        for j in range(i+1, len(neighborhoods)):
            if merged_out[j] == 1:
                continue
            cur_gids = set([g.gid for g in neighborhoods[j].genes])

            if len(pivot_gids.intersection(cur_gids)) >= len(pivot_gids)/2:
                merged_out[j] = 1

    included_files, excluded_files = [], []
    for i in range(len(merged_out)):
        cur_file = neighborhoods[i].source_file.split('/')[-1]
        if merged_out[i] == 0:
            included_files.append(cur_file)
        else:
            excluded_files.append(cur_file)

    return included_files, excluded_files


def _similar_same_order(kplet_1, kplet_2):
    """Check if kplet_1 and kplet_2 are similar. Return boolean

    :param kplet_1: Kplet object of lib.utils.classes.Kplet
    :param kplet_2: Kplet object of lib.utils.classes.Kplet

    :return boolean: True(similar) or False(not similar)
    """

    k = kplet_1.k
    common_codes = len(kplet_1.codes.intersection(kplet_2.codes))

    files_1 = set(kplet_1.files)
    files_2 = set(kplet_2.files)

    common_files = len(files_1.intersection(files_2)) / float(min(len(kplet_1.files), len(kplet_2.files)))

    if k == 5:
        if common_codes == 4 and common_files > 0.5:
            return True
    elif k == 4:
        if common_codes == 3 and common_files > 0.5:
            return True
    elif k == 3:
        if common_codes == 2 and common_files > 0.5:
            return True
    elif k == 2:
        if common_codes == 1:
            return True

    return False


def kplet_in_list(kplet, kplet_list):

    for _kplet in kplet_list:
        if _similar_same_order(_kplet, kplet):
            return True
    return False


def basic_merge_within_orders(kplets):

    merged_kplets = []
    merge_order = []

    merged_out = [0 for _ in range(len(kplets))]

    for i in range(len(kplets)):
        if merged_out[i] == 1:
            continue
        outer_kplet = kplets[i]
        merge_order.append(i)

        if merged_out[i] == 1:
            continue

        to_move = []
        for j in range(i+1, len(kplets)):
            if merged_out[j] == 1:
                continue
            inner_kplet = kplets[j]

            if _similar_same_order(outer_kplet, inner_kplet):
                to_move.append(inner_kplet)
                merged_out[j] = 1
                merge_order.append(j)

        _merged_list = [outer_kplet] + to_move

        kplet_list = KpletList(_merged_list, merge_order)
        merged_kplets.append(kplet_list)

    merged_kplets.sort(key=lambda x: x.weight, reverse=True)
    return merged_kplets


def merge_kplets_within_orders_iterative_3(kplets):
    """ Merge the kplets of same size, if they carry similarity in composition.

    :param kplets: list of objects of lib.utils.classes.Kplet
    :return a: list of list of kplets of same order
    """

    # First round of merging

    new_kplets = list()

    # Second round of merging
    # Iterate the merging procedure until it converges
    # cnt = 1
    # while True:
    #     merged_out = [0 for _ in range(len(merged_kplets))]
    #
    #     gid_pool = []
    #     for merged_list in merged_kplets:
    #         gid_pool.append(set(gid for kplet in merged_list for gid in kplet.gids))
    #
    #     new_merged_kplets = []
    #     for i in range(len(merged_kplets)):
    #         if merged_out[i] == 1:
    #             continue
    #         to_move = []
    #
    #         outer_gids = gid_pool[i]
    #         outer_list = merged_kplets[i]
    #         for j in range(i+1, len(merged_kplets)):
    #             if merged_out[j] == 1:
    #                 continue
    #
    #             inner_gids = gid_pool[j]
    #             inner_list = merged_kplets[j]
    #
    #             common = len(inner_gids.intersection(outer_gids))
    #             smaller = min(len(inner_gids), len(outer_gids))
    #
    #             if not common:
    #                 continue
    #
    #             if float(common)/smaller > 0.8:
    #                 to_move += inner_list
    #                 merged_out[j] = 1
    #
    #         new_merged_kplets.append(outer_list + to_move)
    #     cnt += 1
    #     if len(merged_kplets) == len(new_merged_kplets):
    #         merged_kplets = new_merged_kplets
    #         break
    #
    #     merged_kplets = new_merged_kplets
    # return merged_kplets


def merge_kplets_within_orders_iterative_2(kplet_lists):

    # First round of merging

    # kplet_lists = basic_merge_within_orders(kplets)

    # Second round of merging
    # Iterate the merging procedure until it converges
    cnt = 1

    while True:
        merged_out = [0 for _ in range(len(kplet_lists))]
        new_kplet_lists = []

        for i in range(len(kplet_lists)):
            if merged_out[i] == 1:
                continue

            outer_kplet_list = kplet_lists[i]
            for j in range(i+1, len(kplet_lists)):
                if merged_out[j] == 1:
                    continue

                inner_kplet_list = kplet_lists[j]
                common = len(inner_kplet_list.files.intersection(outer_kplet_list.files))
                smaller = min(len(inner_kplet_list.files), len(outer_kplet_list.files))

                if not common:
                    continue

                if float(common)/smaller > 0.8:
                    outer_kplet_list.merge(inner_kplet_list)
                    merged_out[j] = 1

            new_kplet_lists.append(outer_kplet_list)
        cnt += 1

        if len(kplet_lists) == len(new_kplet_lists):
            kplet_lists = new_kplet_lists
            break

        kplet_lists = new_kplet_lists

    return kplet_lists


def merge_kplets_within_orders_iterative(kplets):
    """ Merge the kplets of same size, if they carry similarity in composition.

    :param kplets: list of objects of lib.utils.classes.Kplet
    :return a: list of list of kplets of same order
    """

    # First round of merging

    target_profiles = t.target_profiles()
    new_kplets = list()

    for kplet in kplets:
        if kplet.codes.intersection(target_profiles):
            new_kplets.append(kplet)

    kplets = new_kplets

    merged_kplets = basic_merge_within_orders(kplets)
    return merged_kplets
    print 'Over return'

    # Second round of merging
    # Iterate the merging procedure until it converges
    # cnt = 1
    # while True:
    #     merged_out = [0 for _ in range(len(merged_kplets))]
    #
    #     gid_pool = []
    #     for merged_list in merged_kplets:
    #         gid_pool.append(set(gid for kplet in merged_list for gid in kplet.gids))
    #
    #     new_merged_kplets = []
    #     for i in range(len(merged_kplets)):
    #         if merged_out[i] == 1:
    #             continue
    #         to_move = []
    #
    #         outer_gids = gid_pool[i]
    #         outer_list = merged_kplets[i]
    #         for j in range(i+1, len(merged_kplets)):
    #             if merged_out[j] == 1:
    #                 continue
    #
    #             inner_gids = gid_pool[j]
    #             inner_list = merged_kplets[j]
    #
    #             common = len(inner_gids.intersection(outer_gids))
    #             smaller = min(len(inner_gids), len(outer_gids))
    #
    #             if not common:
    #                 continue
    #
    #             if float(common)/smaller > 0.8:
    #                 to_move += inner_list
    #                 merged_out[j] = 1
    #
    #         new_merged_kplets.append(outer_list + to_move)
    #     cnt += 1
    #     if len(merged_kplets) == len(new_merged_kplets):
    #         merged_kplets = new_merged_kplets
    #         break
    #
    #     merged_kplets = new_merged_kplets
    # return merged_kplets


def merge_kplets_across_orders(superplets_pool, subplets_pool):
    """ Merge kplets of different size.

    Input arguments:
    superplets -- higher level kplets
    subplets -- lower level kplets

    In case of pentaplets and quadruplets, merge out the quadruplets into the pentaplets, if it happens to be a
    subset of any pentaplet.
    """

    superplet_gi_pool = []
    for superplet_list in superplets_pool:
        tmp_gi_list = set([])
        for kplet in superplet_list:
            tmp_gi_list.update(kplet.gids)
        superplet_gi_pool.append(tmp_gi_list)

    assert len(superplet_gi_pool) == len(superplets_pool)

    merged_out = [[0]*len(subplets_list) for subplets_list in subplets_pool]

    for subplet_outer_ind in range(len(subplets_pool)):
        subplet_list = subplets_pool[subplet_outer_ind]

        for subplet_inner_ind in range(len(subplet_list)):
            cur_subplet = subplet_list[subplet_inner_ind]

            for superplet_ind in range(len(superplets_pool)):

                if cur_subplet.gids.issubset(superplet_gi_pool[superplet_ind]):
                    superplets_pool[superplet_ind].append(cur_subplet)
                    merged_out[subplet_outer_ind][subplet_inner_ind] = 1
                    break

    subplets_pool = [[subplet_list[i] for i in range(len(subplet_list)) if not merged_out[cnt][i]] for cnt, subplet_list in enumerate(subplets_pool)]
    subplets_pool = [l for l in subplets_pool if l]

    return superplets_pool, subplets_pool


def arcog_profiles_pool_into_classes_pool(profile_community):

    _arcog2class = t.map_arcog2class()
    class_community = list()

    for profiles in profile_community:
        classes = dict()
        for k in profiles:
            if k in _arcog2class:
                t.update_dictionary(classes, _arcog2class[k], profiles[k])
        class_community.append(classes)
    return class_community


def arcog_profile_count_into_class_count(community_count):

    _arcog2class = t.map_arcog2class()
    class2count = dict()
    class2profiles = dict()

    for k in community_count:
        if k in _arcog2class:
            _classes = _arcog2class[k]
            for _class in _classes:
                t.update_dictionary(class2count, _class, community_count[k])
                # t.update_dictionary_list_value(class2profiles, _class, k)
                t.update_dictionary(class2profiles, _class, [k])
        else:
            _class = 'Unclassified'
            t.update_dictionary(class2count, _class , community_count[k])
            # t.update_dictionary_list_value(class2profiles, _class, k)
            t.update_dictionary(class2profiles, _class, [k])

    return class2count, class2profiles


def cdd_profile_count_into_class_count(community_count):

    _cdd2class = t.map_cdd2class()
    class_count = dict()
    class2profiles = dict()

    for k in community_count:
        if k in _cdd2class:
            _classes = _cdd2class[k]
            for _class in _classes:
                t.update_dictionary(class_count,    _class, community_count[k])
                t.update_dictionary(class2profiles, _class, [k])
        else:
            _class = 'Unclassified'
            t.update_dictionary(class_count,    _class, community_count[k])
            t.update_dictionary(class2profiles, _class, [k])

    return class_count, class2profiles


def trim_file_summary_list(file_summaries, data_type='bacteria'):

    if data_type == 'bacteria':
        return [c for c in file_summaries if len(c.kplets) > 20]
    elif data_type == 'archaea':
        return [c for c in file_summaries if len(c.kplets) > 10]
    else:
        raise NotImplementedError


