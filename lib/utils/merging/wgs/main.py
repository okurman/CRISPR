__author__ = 'hudaiber'

import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/lib/BioPy/')
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/lib/BioPy/')
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')
import global_variables as gv
import dm_tools as dt
from lib.utils import tools as t
from lib.utils.classes import WGSNeighborhoodFileSummary, WgsKpletList
from lib.utils.arguments import CrisprMergingKplets2FsOutput
from lib.utils.merging.wgs import get_clustered_loci, get_singleton_loci
import os


def kplet_to_file_summaries(kplet, neighborhood_files_path):

    file_summaries = list()
    organisms = set()

    _crispr_type2files = dict()

    for f in kplet.files:
        _genes = dt.get_wgs_file(os.path.join(neighborhood_files_path, f))

        _src = _genes[0].src
        _org = _genes[0].organism

        _crispr_type = _genes[0].crispr_type
        t.update_dictionary_set(_crispr_type2files, _crispr_type, f)

        file_summaries.append(WGSNeighborhoodFileSummary(f,[kplet],_genes,_org,_src))
        organisms.update(_org)

    file_summaries.sort(key=lambda x: x.org)
    return file_summaries, organisms, _crispr_type2files


def kplet_list_to_file_summaries(kplets, neighborhood_files_path, filter_weak_hits=True, dataset=None):

    file_summaries = list()
    organisms = set()
    _crispr_type2files = dict()
    _file2kplets = dict()
    _kplet2count_af = dict() # kplet2count after filtration
    _kplet2count_bf = dict() # kplet2count before filtration

    _profile2count_bf = dict()
    _profile2count_af = dict()

    filter_size = 5

    singletons = get_singleton_loci(dataset)
    clusters = get_clustered_loci(dataset)

    for kplet in kplets:
        for f in kplet.files:
            t.update_dictionary(_file2kplets, f, [kplet])

    initial_length = len(_file2kplets)

    for f in _file2kplets.keys():
        [t.update_dictionary(_kplet2count_bf, kplet.id, 1) for kplet in _file2kplets[f]]
    del f

    kplet_ids = [k.id for k in kplets]

    if filter_weak_hits:
        _file2kplets = {k: v for (k,v) in _file2kplets.items() if len(v) > filter_size}

    if len(_file2kplets) < 2: return None

    _file2genes = {f: dt.get_wgs_file(os.path.join(neighborhood_files_path, f)) for f in _file2kplets.keys()}
    _files = set(_file2kplets.keys())

    for _gene_list in _file2genes.values():
        for _gene in _gene_list:
            for _c in _gene.cogid.split(','):
                t.update_dictionary(_profile2count_bf, _c, 1)
    del _gene_list, _gene, _c

    while _files:
        _f = _files.pop()
        if _f in singletons:
            _genes = _file2genes[_f]
            _src = _genes[0].src
            _org = _genes[0].organism
            _crispr_type = _genes[0].crispr_type
            t.update_dictionary_set(_crispr_type2files, _crispr_type, _f)

            file_summaries.append(WGSNeighborhoodFileSummary(_f, _file2kplets[_f], _genes, _org, _src, 'singleton'))
            organisms.update(set([_org]))

        else:
            _cluster = None
            for cl in clusters:
                if _f in cl.files:
                    _cluster = cl
                    break
            if not _cluster:
                continue
            del cl

            _cl_files = _cluster.files.intersection(_files)
            _representative = _f
            del _f

            for _cl_file in _cl_files:
                if len(_file2genes[_cl_file]) > len(_file2genes[_representative]):
                    _representative = _cl_file

            _genes = _file2genes[_representative]
            _src = _genes[0].src
            _org = _genes[0].organism
            _crispr_type = _genes[0].crispr_type
            t.update_dictionary_set(_crispr_type2files, _crispr_type, _representative)

            _file_summary = WGSNeighborhoodFileSummary(_representative, _file2kplets[_representative], _genes, _org,
                                                       _src, _cluster)
            _file_summary.cluster_local_count = len(_cl_files)+1

            file_summaries.append(_file_summary)
            organisms.update(set([_org]))

            _files = _files.difference(_cl_files)

    file_summaries = [fs for fs in file_summaries if len(fs.kplets)>1]

    _files = [fs.file_name for fs in file_summaries]
    for _f in _files:
        [t.update_dictionary(_kplet2count_af, kplet.id, 1) for kplet in _file2kplets[_f]]

        _gene_list = _file2genes[_f]
        for _gene in _gene_list:
            for _c in _gene.cogid.split(','):
                t.update_dictionary(_profile2count_af, _c, 1)

    file_summaries.sort(key=lambda x: x.org)
    retval = CrisprMergingKplets2FsOutput()
    retval.file_summaries = file_summaries
    retval.organisms = organisms
    retval.crispr_type2files = _crispr_type2files
    retval.kplet2count_af = _kplet2count_af
    retval.kplet2count_bf = _kplet2count_bf
    retval.initial_length = initial_length
    retval.kplets = kplets
    retval.profile2count_bf = _profile2count_bf
    retval.profile2count_af = _profile2count_af

    return retval


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
    # print len(files_1), len(files_2), len(files_1.intersection(files_2)) ,common_files

    if k == 5:
        if (common_codes == 4 and common_files > 0.9) or common_files == 1.0:
            return True
    elif k == 4:
        if (common_codes == 3 and common_files > 0.9) or common_files == 1.0:
            return True
    elif k == 3:
        if (common_codes == 2 and common_files > 0.9) or common_files == 1.0:
            return True
    elif k == 2:
        if (common_codes == 1 and common_files > 0.9) or common_files == 1.0:
            return True

    return False


def basic_merge_within_orders(kplets):

    merged_kplets = []

    merged_out = [0] * len(kplets)

    for i in range(len(kplets)):
        if merged_out[i] == 1:
            continue
        outer_kplet = kplets[i]

        to_move = []
        for j in range(i+1, len(kplets)):
            if merged_out[j] == 1:
                continue
            inner_kplet = kplets[j]
            if _similar_same_order(outer_kplet, inner_kplet):
                to_move.append(inner_kplet)
                merged_out[j] = 1

        _merged_list = [outer_kplet] + to_move

        kplet_list = WgsKpletList(_merged_list)
        merged_kplets.append(kplet_list)

    merged_kplets.sort(key=lambda x: x.count, reverse=True)
    return merged_kplets


def merge_kplets_within_orders_iterative(kplet_lists):

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

                if float(common)/avg > 0.5:
                    outer_kplet_list.merge(inner_kplet_list)
                    merged_out[j] = 1

            new_kplet_lists.append(outer_kplet_list)
        cnt += 1
        kplet_lists = new_kplet_lists
        if len(kplet_lists) == len(new_kplet_lists):
            break


    return kplet_lists
