__author__ = 'Sanjarbek Hudaiberdiev'

from lib.db import DbClass
from lib.utils.classes import Kplet
from lib.db.archea import neighborhoods_path


def get_multiple_kplets():

    _sql_cmd = """select  ap.id, count(*) cnt, group_concat(convert(apw.file_id, char(15))) as file_ids
                  from archea_5plets ap
                  inner join archea_5plets_win10 apw on ap.id = apw.kplet_id
                  group by ap.id
                  having count(*)>1
                  order by cnt desc"""
    _db = DbClass()
    _db.cmd = _sql_cmd

    return _db.retrieve()


def get_archaea_kplets():

    sql_cmd = """select apc.*, s1.cnt, s1.wgt
                 from (
                        select ap.id ,count(*) as cnt, sum(w.weight) as wgt
                        from archea_5plets ap
                        inner join archea_5plets_win10 apw on ap.id = apw.kplet_id
                        inner join archea_win10_files awf on apw.file_id = awf.id
                        inner join sources s on awf.source_id=s.id
                        inner join weights w on w.genome_id=s.genome_id
                        group by ap.id ) s1
                 inner join archea_5plets_codes apc on s1.id=apc.id
                 order by s1.wgt desc"""

    _db = DbClass()
    _db.cmd = sql_cmd

    return _db.retrieve()


def archea_kplet_ids2files(id_list):

    _sql_cmd = """select distinct awf.name
                  from archea_5plets ap
                  inner join archea_5plets_win10 apw on ap.id = apw.kplet_id
                  inner join archea_win10_files awf on apw.file_id = awf.id
                  where ap.id in (%s)""" % " , ".join(id_list)

    _db = DbClass()
    _db.cmd = _sql_cmd

    return [l[0] for l in _db.retrieve()]


def get_code_kplet(kplet_id):
    _sql_cmd = """select kplet_1, kplet_2, kplet_3, kplet_4, kplet_5 from archea_5plets_codes where id = %s""" % kplet_id
    _db = DbClass()
    _db.cmd = _sql_cmd
    return _db.retrieve()


def get_report_kplets(limit_to=300, load_locations=None):

    sql_cmd = """select apc.*, s1.cnt, s1.wgt, s1.an
                from (
                        select ap.id ,count(*) as cnt, sum(w.weight) as wgt, group_concat(awf.name) as an
                        from archea_5plets ap
                        inner join archea_5plets_win10 apw on ap.id = apw.kplet_id
                        inner join archea_win10_files awf on apw.file_id = awf.id
                        inner join sources s on awf.source_id=s.id
                        inner join weights w on w.genome_id=s.genome_id
                        group by ap.id
                        having count(distinct s.genome_id)>1 ) s1
                inner join archea_5plets_codes apc on s1.id=apc.id
                order by s1.wgt desc
                limit 0,%d""" % limit_to

    _db = DbClass()
    _db.cmd = sql_cmd
    out_list = []

    for row in _db.retrieve():
        id = row[0]
        kplet_codes = (row[1:6])
        if len(set(kplet_codes)) != 5:
            continue
        count = row[6]
        files = row[8].split(',')
        tmp_kplet = Kplet(id=id, codes=kplet_codes, count=count, files=files)
        out_list.append(tmp_kplet)

    _path = neighborhoods_path()
    if load_locations:
        [kplet.load_locations(_path) for kplet in out_list]

    return out_list
