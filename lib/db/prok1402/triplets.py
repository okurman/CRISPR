__author__ = 'Sanjarbek Hudaiberdiev'

from lib.db import DbClass
from lib.utils.classes import Kplet
from lib.db.bacteria import neighborhoods_path


def store_kplets_pile(kplets_pile, cdd2id, file2id):

    _sql_kplet = """insert ignore into bacteria_3plets (kplet_1, kplet_2, kplet_3) values \n"""

    _sql_kplet_file = """insert ignore into bacteria_3plets_win10 (kplet_id, file_id) values \n"""

    for (kplets, fname) in kplets_pile:

        for kplet in kplets:
            kplet = list(kplet)
            kplet.sort()
            kplet = tuple([int(cdd2id[k]) for k in kplet])

            _sql_kplet += """(%d, %d, %d),\n""" % kplet

            _sql_kplet_file += ("""((select id from bacteria_3plets where """ +
                                """kplet_1=%d and kplet_2=%d and kplet_3=%d),""" +
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

def get_multiple_kplets():

    _db = DbClass()
    _db.cmd = "SET group_concat_max_len = 100000;"
    _db.execute()
    _db.cmd = """ select  ap.id, count(*) cnt, group_concat(convert(apw.file_id, char(15))) as file_ids
                  from bacteria_3plets ap
                  inner join bacteria_3plets_win10 apw on ap.id = apw.kplet_id
                  group by ap.id
                  having count(*)>1
                  order by cnt desc"""

    return _db.retrieve()


def get_code_kplet(kplet_id, id2cdd=None):

    _db = DbClass()

    if not id2cdd:
        _db.cmd = """select cp1.code, cp2.code, cp3.code
                from bacteria_3plets bp
                inner join cdd_profiles cp1 on cp1.id = bp.kplet_1
                inner join cdd_profiles cp2 on cp2.id = bp.kplet_2
                inner join cdd_profiles cp3 on cp3.id = bp.kplet_3
                where bp.id = %d""" % kplet_id
        retval = _db.retrieve()[0]

    else:

        _db.cmd = """select kplet_1, kplet_2, kplet_3
                     from bacteria_3plets where id = %d""" % kplet_id

        retval = _db.retrieve()[0]
        retval = set([id2cdd[id] for id in retval])

    return retval


def get_report_kplets(id2cdd, limit_to=500, load_locations=None):

    _db = DbClass()
    _db.cmd = """SET group_concat_max_len=1500000"""
    _db.execute()

    _db.cmd = """select ap.* ,count(*) as cnt, sum(w.weight) as wgt, group_concat(awf.name) as an
                 from bacteria_3plets ap
                 inner join bacteria_3plets_win10 apw on ap.id = apw.kplet_id
                 inner join bacteria_win10_files awf on apw.file_id = awf.id
                 inner join sources s on awf.source_id=s.id
                 inner join weights w on w.genome_id=s.genome_id
                 group by ap.id
                 having count(distinct s.genome_id)>1
                 order by wgt desc
                 limit 0, %d""" % limit_to

    out_list = []

    for row in _db.retrieve():
        id = row[0]
        kplet_codes = ([id2cdd[int(_id)] for _id in row[1:4]])
        if len(set(kplet_codes)) != 3:
            continue
        count = row[4]
        weight = row[5]
        files = row[6].split(',')
        tmp_kplet = Kplet(id=id, codes=kplet_codes, count=count, files=files)
        out_list.append(tmp_kplet)

    _path = neighborhoods_path()
    if load_locations:
        [kplet.load_locations(_path) for kplet in out_list]

    return out_list
