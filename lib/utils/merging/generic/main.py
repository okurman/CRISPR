import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')
import global_variables as gv
import dm_tools as dt
from lib.utils import tools as t
from lib.utils.classes import NeighborhoodFileSummary, KpletList

from lib.utils.arguments import GenericMergingKplets2FsOutput


import os

# Globals
_gid2arcog_cdd = t.map_gid2arcog_cdd()
_org2weight = t.map_genome2weight()

cas_type_file = os.path.join(gv.project_data_path, 'cas1402/cas1402.type.tab')
_gi2castype =  {l.strip().split()[0]: l.strip().split()[1].split(';') for l in open(cas_type_file).readlines()}


def _similar_same_order(kplet_1, kplet_2, loci_threshold=0.5):
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
        if common_codes == 4 and common_files > loci_threshold:
            return True
    elif k == 4:
        if common_codes == 3 and common_files > loci_threshold:
            return True
    elif k == 3:
        if common_codes == 2 and common_files > loci_threshold:
            return True
    elif k == 2:
        if common_codes == 1 and common_files > loci_threshold:
            return True

    return False


def kplet_in_list(kplet, kplet_list):

    for _kplet in kplet_list:
        if _similar_same_order(_kplet, kplet):
            return True
    return False


def basic_merge_within_orders(kplets, loci_threshold):

    merged_kplets = []

    merged_out = [0 for _ in range(len(kplets))]

    for i in range(len(kplets)):
        if merged_out[i] == 1:
            continue
        outer_kplet = kplets[i]

        if merged_out[i] == 1:
            continue

        to_move = []
        for j in range(i+1, len(kplets)):
            if merged_out[j] == 1:
                continue
            inner_kplet = kplets[j]

            if _similar_same_order(outer_kplet, inner_kplet, loci_threshold):
                to_move.append(inner_kplet)
                merged_out[j] = 1

        _merged_list = [outer_kplet] + to_move

        kplet_list = KpletList(_merged_list)
        merged_kplets.append(kplet_list)

    return merged_kplets


def merge_kplets_within_orders_iterative(kplet_lists, loci_threshold=0.5):

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
                avg = (len(inner_kplet_list.files) + len(outer_kplet_list.files)) / 2.0

                if not common:
                    continue

                if float(common)/avg > loci_threshold:
                    outer_kplet_list.merge(inner_kplet_list)
                    merged_out[j] = 1

            new_kplet_lists.append(outer_kplet_list)
        cnt += 1
        kplet_lists = new_kplet_lists
        if len(kplet_lists) == len(new_kplet_lists):
            break

    return kplet_lists


def kplet_list_to_file_summaries(kplets, neighborhood_files_path, filter_weak_hits=True):

    file_summaries = list()
    organisms = set()

    _file2kplets = dict()
    _kplet2count_af = dict() # kplet2count after filtration
    _kplet2count_bf = dict() # kplet2count before filtration

    _profile2count_bf = dict()
    _profile2count_af = dict()
    _cas_type2count = dict()
    # filter_size = 5

    for kplet in kplets:
        for f in kplet.files:
            t.update_dictionary(_file2kplets, f, [kplet])

    initial_length = len(_file2kplets)

    for f in _file2kplets.keys():
        [t.update_dictionary(_kplet2count_bf, kplet.id, 1) for kplet in _file2kplets[f]]
    del f

    # if filter_weak_hits:
    #     _file2kplets = {k: v for (k,v) in _file2kplets.items() if len(v) > filter_size}

    if len(_file2kplets) < 2: return None

    _file2genes = {f: dt.get_pty_file_generic(os.path.join(neighborhood_files_path, f)) for f in _file2kplets.keys()}
    _files = set(_file2kplets.keys())

    for _f in _files:

        _genes = _file2genes[_f]
        _src = _genes[0].src
        _org = _genes[0].organism

        organisms.update([_org])
        _nfs = NeighborhoodFileSummary(_f, _file2kplets[_f], _genes, _org, _src, _org2weight[_org])

        for _gene in _genes:
            if _gene.gid in _gi2castype:
                _nfs.cas_type = ";".join(_gi2castype[_gene.gid])
                for _cas_type in _gi2castype[_gene.gid]:
                    t.update_dictionary(_cas_type2count,_cas_type,1)
                break

        [t.update_dictionary(_kplet2count_af, kplet.id, 1) for kplet in _file2kplets[_f]]

        for _gene in _genes:
            for _c in _gene.cogid.split(','):
                t.update_dictionary(_profile2count_af, _c, 1)


        file_summaries.append(_nfs)
    # file_summaries = [fs for fs in file_summaries if len(fs.kplets)>1]
    # _files = [fs.file_name for fs in file_summaries]

    # for _f in _files:
        # [t.update_dictionary(_kplet2count_af, kplet.id, 1) for kplet in _file2kplets[_f]]
        #
        # _gene_list = _file2genes[_f]
        # for _gene in _gene_list:
        #     for _c in _gene.cogid.split(','):
        #         t.update_dictionary(_profile2count_af, _c, 1)

    file_summaries.sort(key=lambda x: x.org)
    retval = GenericMergingKplets2FsOutput()
    retval.file_summaries = file_summaries
    retval.organisms = organisms
    retval.profile2count_bf = _profile2count_bf
    retval.profile2count_af = _profile2count_af
    retval.kplet2count_af = _kplet2count_af
    retval.kplet2count_bf = _kplet2count_bf
    retval.weight = sum(fs.weight for fs in file_summaries)
    retval.cas_type2count = _cas_type2count

    return retval
