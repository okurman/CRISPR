__author__ = 'Sanjarbek Hudaiberdiev'

from lib.db import DbClass
from lib.utils.classes import Kplet



def store_kplets_pile(kplets_pile, profile2id, file2id):

    kplet_table_name = "crispr_5plets"
    kplet_table_file_name = "crispr_5plets_wgs"

    _sql_kplet = """insert ignore into %s (kplet_1, kplet_2, kplet_3, kplet_4, kplet_5) values \n""" % kplet_table_name
    _sql_kplet_file = """insert ignore into %s (kplet_id, file_id) values \n""" % kplet_table_file_name

    for (kplets, fname) in kplets_pile:
        for kplet in kplets:
            kplet = list(kplet)
            kplet.sort()

            kplet = tuple([int(profile2id[k]) for k in kplet])

            _sql_kplet += """(%d, %d, %d, %d, %d),\n""" % kplet

            _sql_kplet_file += ("""((select id from crispr_5plets where """ +
                                """kplet_1=%d and kplet_2=%d and kplet_3=%d and kplet_4=%d and kplet_5=%d),""" +
                                """%d),\n""") % (kplet + (int(file2id[fname]),))

    _sql_kplet = _sql_kplet[:-2]
    _sql_kplet += ';'

    _sql_kplet_file = _sql_kplet_file[:-2]
    _sql_kplet_file += ';'

    _db = DbClass()

    _db.cmd = _sql_kplet
    _db.execute()
    _db.commit()

    _db.cmd = _sql_kplet_file
    _db.execute()
    _db.commit()


def get_report_kplets(id2cdd, limit_to=500):

    _db = DbClass()
    _db.cmd = """SET group_concat_max_len=15000000"""
    _db.execute()

    _db.cmd = """select cp.*, count(distinct cwf.id) as cnt, group_concat(cwf.name) as files
                from crispr_5plets cp
                inner join crispr_5plets_wgs cpw on cp.id = cpw.kplet_id
                inner join crispr_wgs_files cwf on cpw.file_id = cwf.id
                group by cp.id
                having count(distinct cwf.id)>1
                order by cnt desc
                limit 0, %d""" % limit_to

    out_list = []

    for row in _db.retrieve():
        id = row[0]
        kplet_codes = ([id2cdd[int(_id)] for _id in row[1:6]])
        if len(set(kplet_codes)) != 5:
            continue
        count = row[6]
        files = row[7].split(',')
        tmp_kplet = Kplet(id=id, codes=kplet_codes, count=count, files=files)
        out_list.append(tmp_kplet)

    return out_list