__author__ = 'hudaiber'

import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
    sys.path.append('/Users/hudaiber/Projects/BioPy/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')
    sys.path.append('/hudaiber/hudaiber/Projects/BioPy/')
import global_variables as gv

sys.path.append('../')
import os
import dm_tools as dt
import lib.utils.tools as t
import lib.utils.reporting as r
import lib.utils.search_in_genomes as sig
from lib.db import map_id2cdd_clusters
from lib.db.bacteria import pentaplets as p


if __name__=='__main__':

    work_path = os.path.join(gv.project_data_path, 'Bacteria/cases')
    pty_path = gv.pty_data_path

    kplet_id = 306123

    id2cdd = map_id2cdd_clusters()
    kplet = p.get_report_kplet(kplet_id, id2cdd, load_locations=True)

    target_profiles = set(t.bacteria_target_profiles())
    dump_file = os.path.join(work_path, 'kplet.p.bz2')
    t.dump_compressed_pickle(dump_file, kplet)

    kplet = t.load_compressed_pickle(dump_file)
    kplet_codes = kplet.codes.difference(target_profiles)

    org2src, src2blocks = sig.search_kplet_in_genomes(kplet_codes, target_profiles, max_dist=4)

    # dump_file = os.path.join(work_path, 'org2src_global.p.bz2')
    # t.dump_compressed_pickle(dump_file, org2src)
    # dump_file = os.path.join(work_path, 'src2blocks_global.p.bz2')
    # t.dump_compressed_pickle(dump_file, src2blocks)
    #
    # dump_file = os.path.join(work_path, 'org2src_global.p.bz2')
    # org2src = t.load_compressed_pickle(dump_file)
    # dump_file = os.path.join(work_path, 'src2blocks_global.p.bz2')
    # src2blocks = t.load_compressed_pickle(dump_file)

    xls_file = os.path.join(work_path,'kplet_%d_tight.xlw'%kplet_id)
    r.write_to_xls_extra_search(xls_file,org2src,src2blocks,kplet)