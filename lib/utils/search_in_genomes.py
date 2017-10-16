__author__ = 'hudaiber'

import os
import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
    sys.path.append('/Users/hudaiber/Projects/BioPy/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')
    sys.path.append('/hudaiber/hudaiber/Projects/BioPy/')

import global_variables as gv
import lib.utils.tools as t
import dm_tools as dt

gid2cdd = t.map_gid2arcog_cdd()

def search_kplet_in_genomes(kplet_codes, target_profiles, max_dist=20, block_size=4):

    org2src = dict()
    src2blocks = dict()

    for _org in os.listdir(gv.pty_data_path):

        _org_path = os.path.join(gv.pty_data_path, _org)
        for _src in os.listdir(_org_path):
            _genes = dt.get_pty(os.path.join(_org_path, _src))

            blocks = list()
            cur_block = list()

            last_ind = None

            for (ind, _gene) in enumerate(_genes):
                _cogids = set(gid2cdd[_gene.gid].split() if _gene.gid in gid2cdd else set([]))
                _gene.cogid = _cogids
                if _cogids.intersection(kplet_codes):
                    if not last_ind:
                        cur_block.append(_gene)
                        last_ind = ind+1
                        continue

                    if ind - last_ind < max_dist:
                        cur_block += _genes[last_ind: ind+1]
                    else:
                        blocks.append(cur_block)
                        cur_block = [_gene]
                    last_ind = ind+1
            blocks.append(cur_block)

            filtered_blocks = list()

            for block in blocks:
                block_codes = set([])
                block_all_codes = set([])

                for _gene in block:
                    block_codes.update(_gene.cogid.intersection(kplet_codes))
                    block_all_codes.update(_gene.cogid)

                if len(block_codes) >= block_size and not block_all_codes.intersection(target_profiles):
                    filtered_blocks.append(block)
            del blocks

            if filtered_blocks:
                t.update_dictionary_set(org2src, _org, _src)
                src2blocks[_src] = filtered_blocks

    return org2src, src2blocks
