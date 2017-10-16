__author__ = 'Sanjarbek Hudaiberdiev'

from lib.db import DbClass
import db_tools as t
from lib.utils.classes import Kplet
from lib.db.archea import neighborhoods_path


def store_kplets(kplets, fname):

    for kplet in kplets:
        l = list(kplet)
        l.sort()
        kplet_id = get_kplet_id(l)

        if not kplet_id:
            insert_kplet(l)
            kplet_id = get_kplet_id(l)
        connection.commit()

        file_id = t.get_file_id(fname)
        insert_kplet_file(kplet_id, file_id)


def insert_kplet_file(kplet_id, file_id):
    sql_cmd = """insert ignore into archea_2plets_win10 values (%d, %d)"""
    sql_cmd = sql_cmd % (kplet_id, file_id)

    cursor = setup_cursor()
    cursor.execute(sql_cmd)
    connection.commit()


def get_kplet_id(profiles):
    profile2id = t.map_profile2id(profiles)

    sql_cmd = """select id from archea_2plets where kplet_1=%d and kplet_2=%d"""
    sql_cmd = sql_cmd % (profile2id[profiles[0]], profile2id[profiles[1]])

    cursor = setup_cursor()
    cursor.execute(sql_cmd)
    rows = cursor.fetchall()

    if rows:
        return rows[0][0]
    return None


def insert_kplet(profiles):
    profile2id = t.map_profile2id(profiles)

    sql_cmd = """insert into archea_2plets (kplet_1, kplet_2) values (%d, %d)"""
    sql_cmd = sql_cmd % (profile2id[profiles[0]], profile2id[profiles[1]])

    cursor = setup_cursor()
    cursor.execute(sql_cmd)


def get_code_kplet(kplet_id):
    _sql_cmd = """SELECT cp1.code, cp2.code
                    FROM PatternQuest.archea_2plets ap
                    inner join cdd_profiles cp1 on cp1.id = ap.kplet_1
                    inner join cdd_profiles cp2 on cp2.id = ap.kplet_2
                    where ap.id = %d""" % kplet_id
    _cursor = setup_cursor()
    _cursor.execute(_sql_cmd)
    return _cursor.fetchall()[0]


def get_multiple_kplets():

    _sql_cmd = """select  ap.id, count(*) cnt, group_concat(convert(apw.file_id, char(15))) as file_ids
                  from archea_2plets ap
                  inner join archea_2plets_win10 apw on ap.id = apw.kplet_id
                  group by ap.id
                  having count(*)>1
                  order by cnt desc"""

    _cursor = setup_cursor()
    _cursor.execute(_sql_cmd)

    return _cursor.fetchall()


def get_report_kplets(limit_to=300, count_filter=1, load_locations=None):

    _db = DbClass()

    _db.cmd = """SET group_concat_max_len=15000"""
    _db.execute()

    _db.cmd = """select apc.*, s1.cnt, s1.wgt, s1.an
                from (
                        select ap.id ,count(*) as cnt, sum(w.weight) as wgt, group_concat(awf.name) as an
                        from archea_2plets ap
                        inner join archea_2plets_win10 apw on ap.id = apw.kplet_id
                        inner join archea_win10_files awf on apw.file_id = awf.id
                        inner join sources s on awf.source_id=s.id
                        inner join weights w on w.genome_id=s.genome_id
                        group by ap.id
                        having count(distinct s.genome_id) > %d ) s1
                inner join archea_2plets_codes apc on s1.id=apc.id
                order by s1.wgt desc
                limit 0,%d""" % (count_filter, limit_to)

    out_list = []

    for row in _db.retrieve():
        id = row[0]
        kplet_codes = (row[1:3])
        if len(set(kplet_codes)) != 2:
            continue
        count = row[3]
        files = row[5].split(',')
        tmp_kplet = Kplet(id=id, codes=kplet_codes, count=count, files=files)
        out_list.append(tmp_kplet)


    if load_locations:
        _path = neighborhoods_path()
        [kplet.load_locations(_path) for kplet in out_list]

    return out_list