__author__ = 'Sanjarbek Hudaiberdiev'

import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')

import os
import global_variables as gv
sys.path.append('../')
from lib.utils import tools
from lib.db.bacteria import db_tools
from lib.db.bacteria import duplets as db
import time


def kplet_gids(neighborhood, kplet_codes):

    out_gids = set()

    for term in kplet_codes:
        for g in neighborhood.genes:
            if term in g.cogid:
                out_gids.add(g.gid)
                break
    return out_gids


if __name__=='__main__':
    tic = time.time()
    print 'reading neighborhoods'
    files_path = os.path.join(gv.project_data_path, 'Bacteria/genes_and_flanks/win_10/raw_nbr_files/')
    neighborhoods = tools.load_neighborhoods(files_path)

    print 'loading kplets'
    kplets = db.get_multiple_kplets()
    # p.dump(kplets, open('multiple_kplets.p', 'w'))
    print 'loading file maps'
    file_id2name = db_tools.map_file_id2name()
    name2file_id = db_tools.map_file_name2id()
    # kplets = p.load(open('multiple_kplets.p'))

    cnt = 0
    cnt_2 = 0
    duplicate_ids = []
    print 'Total kplets:', len(kplets)
    for kplet in kplets:

        kplet_id = kplet[0]
        kplet_file_ids = kplet[2].split(',')

        kplet_file_names = [file_id2name[id] for id in kplet_file_ids]
        kplet_file_names.sort()

        kplet_codes = db.get_code_kplet(kplet_id)
        neighborhoods = tools.load_neighborhoods(files_path, kplet_file_names)
        n_files = [n.source_file.split('/')[-1] for n in neighborhoods]
        files_to_remove = []

        cur_gids = set([])
        for nbr in neighborhoods:
            if not cur_gids:
                cur_gids = kplet_gids(nbr, kplet_codes)
                continue

            tmp_gids = kplet_gids(nbr, kplet_codes)

            if tmp_gids == cur_gids:
                files_to_remove.append(nbr.source_file)
            else:
                cur_gids = tmp_gids

        if files_to_remove:
            files_to_remove = [l.split('/')[-1] for l in files_to_remove]
            file_ids_to_remove = [name2file_id[name] for name in files_to_remove]
            duplicate_ids += [(kplet_id, file_id) for file_id in file_ids_to_remove]

            cnt_2 += 1
            if cnt_2 < 10:
                print kplet_id, files_to_remove
                print kplet_id, file_ids_to_remove
                print

        cnt += 1
        if cnt % 1000 == 0:
            print cnt

    fout = open('bacteria_remove_duplicate_2plets_script.sql', 'w')

    fout.write("delete from bacteria_2plets_win10 where (kplet_id, file_id) in (\n")
    for comb in duplicate_ids:
        fout.write("(%s, %s),\n" % comb)
    fout.write(");")
    fout.close()
