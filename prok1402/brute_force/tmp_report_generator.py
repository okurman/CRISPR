#!/sw/bin/python2.7
__author__ = 'Sanjarbek Hudaiberdiev'

import os
import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')

import global_variables as gv
sys.path.append(gv.project_code_path)


from lib.db.bacteria import neighborhoods_path
from lib.db.bacteria.db_tools import file2src_src2org_map
import lib.utils.merging as merging
import lib.utils.tools as t
import lib.utils.reporting as r

import cPickle
import bz2


if __name__ == '__main__':

    print 'Pre-Loading dictionaries'
    target_profiles = t.bacteria_target_profiles()
    profile2def = t.map_cdd_profile2def()
    gid2arcog_cdd = t.map_gid2arcog_cdd()
    neighborhood_files_path = neighborhoods_path()
    # profile_id2code = map_id2cdd()
    # pickle.dump(profile_id2code, open('profile_id2code.p','w'))
    profile_id2code = cPickle.load(open('/Users/hudaiber/Projects/NewSystems/code/Bacteria/profile_id2code.p'))

    # fname = '/Users/hudaiber/Projects/NewSystems/data/Bacteria/pickle/100000/pentaplets_merged_across.p.bz2'
    fname = '/Users/hudaiber/Projects/NewSystems/data/Bacteria/pickle/100000/merged_kplets.p.bz2'
    f = bz2.BZ2File(fname, 'rb')
    merged_kplets = cPickle.load(f)

    dest_dir = '/Users/hudaiber/Projects/NewSystems/data/Bacteria/reports/tmp/'

    print 'starting reports'
    for cnt, kplet_sublist in enumerate(merged_kplets):
        print cnt+1
        xls_file_name = os.path.join(dest_dir,  "%d.xls" % (cnt+1))

        src2org, file_summaries, community = merging.merge_into_file_summaries(kplet_sublist, neighborhood_files_path, file2src_src2org_map, 'bacteria')
        params = {}
        params['xls_file_name'] = xls_file_name
        params['src2org'] = src2org
        params['file_summaries'] = file_summaries
        params['community'] = community
        params['target_profiles'] = target_profiles
        params['profile2def'] = profile2def
        params['gid2arcog_cdd'] = gid2arcog_cdd

        r.write_to_xls(params)
        break
        if cnt == 200:
            break
