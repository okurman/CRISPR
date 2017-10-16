__author__ = 'Sanjarbek Hudaiberdiev'

from lib.db import DbClass
# from lib.utils.classes import Kplet
from lib.db.prok1402 import neighborhoods_path
import os
import sys


def store_kplets_pile(kplets_pile, cdd2id, file2id):

    _sql_kplet = """insert ignore into bacteria_2plets (kplet_1, kplet_2) values \n"""

    _sql_kplet_file = """insert ignore into bacteria_2plets_win10 (kplet_id, file_id) values \n"""

    for (kplets, fname) in kplets_pile:

        for kplet in kplets:
            kplet = list(kplet)
            kplet.sort()
            kplet = tuple([int(cdd2id[k]) for k in kplet])

            _sql_kplet += """(%d, %d),\n""" % kplet

            _sql_kplet_file += ("""((select id from bacteria_2plets where """ +
                                """kplet_1=%d and kplet_2=%d),""" +
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
    _db.cmd = "SET group_concat_max_len = 10000000;"
    _db.execute()
    _db.cmd = """ select  ap.id, count(*) cnt, group_concat(convert(apw.file_id, char(15))) as file_ids
                  from bacteria_2plets ap
                  inner join bacteria_2plets_win10 apw on ap.id = apw.kplet_id
                  group by ap.id
                  having count(*)>1
                  order by cnt desc"""

    return _db.retrieve()


def get_code_kplet(kplet_id, id2cdd=None):

    _db = DbClass()

    if not id2cdd:
        _db.cmd = """select cp1.code, cp2.code
                from bacteria_2plets bp
                inner join cdd_profiles cp1 on cp1.id = bp.kplet_1
                inner join cdd_profiles cp2 on cp2.id = bp.kplet_2
                where bp.id = %d""" % kplet_id
        retval = _db.retrieve()[0]

    else:

        _db.cmd = """select kplet_1, kplet_2
                     from bacteria_2plets where id = %d""" % kplet_id

        retval = _db.retrieve()[0]
        retval = set([id2cdd[id] for id in retval])

    return retval


def get_report_kplets(id2cdd, limit_to=500, load_locations=None):

    _db = DbClass()
    _db.cmd = """SET group_concat_max_len=1500000"""
    _db.execute()

    _db.cmd = """select ap.* ,count(*) as cnt, sum(w.weight) as wgt, group_concat(awf.name) as an
                 from bacteria_2plets ap
                 inner join bacteria_2plets_win10 apw on ap.id = apw.kplet_id
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
        kplet_codes = ([id2cdd[int(_id)] for _id in row[1:3]])
        if len(set(kplet_codes)) != 2:
            continue
        count = row[3]
        weight = row[4]
        files = row[5].split(',')
        tmp_kplet = Kplet(id=id, codes=kplet_codes, count=count, files=files)
        out_list.append(tmp_kplet)

    _path = neighborhoods_path()
    if load_locations:
        [kplet.load_locations(_path) for kplet in out_list]

    return out_list


#######################################################################################################################
################           Graph based approach starts here          ##################################################
#######################################################################################################################


def insert_adjacent_duplets(kplets_pile, profile2id, file2id):

    sql_insert_kplet = """insert ignore into prok1402_baited_adj_duplet (kplet_1, kplet_2, multidomain) values \n"""
    sql_insert_kplet_file = """insert ignore into prok1402_baited_adj_duplet_file (kplet_id, file_id) values \n"""

    for (duplet, fname, multidomain) in kplets_pile:

        terms = sorted(list(duplet))
        terms = [profile2id[_] for _ in terms]

        fname = os.path.basename(fname)
        file_id = file2id[fname]

        sql_insert_kplet += """(%s, %s, %d),\n""" % (terms[0], terms[1], multidomain)
        sql_insert_kplet_file += """((select id from prok1402_baited_adj_duplet where """ + \
                                  """kplet_1=%s and kplet_2=%s and multidomain=%d), %d),\n""" \
                                    % (terms[0], terms[1], multidomain ,file_id)

    sql_insert_kplet_file = sql_insert_kplet_file[:-2]
    sql_insert_kplet_file += ";"

    sql_insert_kplet = sql_insert_kplet[:-2]
    sql_insert_kplet += ";"

    _db = DbClass()

    _db.cmd = sql_insert_kplet
    _db.execute()

    _db.cmd = sql_insert_kplet_file
    _db.execute()

    _db.commit()


def insert_duplets(kplets_pile, profile2id, file2id):

    sql_insert_kplet = """insert ignore into prok1402_baited_duplet (kplet_1, kplet_2) values \n"""
    sql_insert_kplet_file = """insert ignore into prok1402_baited_duplet_file (kplet_id, file_id) values \n"""

    for (duplet, fname) in kplets_pile:

        terms = sorted(list(duplet))
        terms = [profile2id[_] for _ in terms]

        fname = os.path.basename(fname)
        file_id = file2id[fname]

        sql_insert_kplet += """(%s, %s),\n""" % (terms[0], terms[1])

        sql_insert_kplet_file += """((select id from prok1402_baited_duplet where """ + \
                                  """kplet_1=%s and kplet_2=%s), %d),\n""" % (terms[0], terms[1] ,file_id)

    sql_insert_kplet_file = sql_insert_kplet_file[:-2]
    sql_insert_kplet_file += ";"

    sql_insert_kplet = sql_insert_kplet[:-2]
    sql_insert_kplet += ";"

    _db = DbClass()

    _db.cmd = sql_insert_kplet
    _db.execute()

    _db.cmd = sql_insert_kplet_file
    _db.execute()

    _db.commit()


def insert_source_duplets(kplets_pile, profile2id, source_id):

    sql_insert_kplet = """insert ignore into prok1402_all_adj_duplet (kplet_1, kplet_2, multidomain) values \n"""
    sql_insert_kplet_file = """insert into prok1402_all_adj_duplet_source (kplet_id, source_id) values \n"""

    for (duplet, fname, multidomain) in kplets_pile:

        terms = sorted(list(duplet))
        terms = [profile2id[_] for _ in terms]

        sql_insert_kplet += """(%s, %s, %d),\n""" % (terms[0], terms[1], multidomain)

        sql_insert_kplet_file += """((select id from prok1402_all_adj_duplet where """ + \
                                  """kplet_1=%s and kplet_2=%s and multidomain=%d), %d),\n""" % \
                                  (terms[0], terms[1], multidomain, source_id)

    sql_insert_kplet_file = sql_insert_kplet_file[:-2]
    sql_insert_kplet_file += ";"

    sql_insert_kplet = sql_insert_kplet[:-2]
    sql_insert_kplet += ";"

    _db = DbClass()

    _db.cmd = sql_insert_kplet
    _db.execute()

    _db.cmd = sql_insert_kplet_file
    _db.execute()

    _db.commit()



def extract_baited_duplet_aggregates(duplet_type="adj"):

    if duplet_type == "adj":
        sql_cmd = """select d.id, 
                        p1.code as profile_1, 
                        p2.code as profile_2, 
                        count(*) as loci, 
                        sum(w.weight) as loci_weight, 
                        group_concat(CONVERT(df.file_id, CHAR(20)) separator ',') as files  
                    from prok1402_baited_adj_duplet d 
                    inner join prok1402_baited_adj_duplet_file df on d.id = df.kplet_id 
                    inner join cdd_profiles p1 on d.kplet_1 = p1.id 
                    inner join cdd_profiles p2 on d.kplet_2 = p2.id 
                    inner join prok1402_baited_files bf on bf.id = df.file_id 
                    inner join sources s on bf.source_id = s.id 
                    inner join genomes g on g.id = s.genome_id 
                    inner join weights w on w.genome_id = g.id 
                    group by d.id 
                    order by loci desc;"""
    else:
        sql_cmd = """select d.id, 
                        p1.code as profile_1, 
                        p2.code as profile_2, 
                        count(*) as loci, 
                        sum(w.weight) as loci_weight, 
                        group_concat(CONVERT(df.file_id, CHAR(20)) separator ',') as files 
                    from prok1402_baited_duplet d 
                    inner join prok1402_baited_duplet_file df on d.id = df.kplet_id
                    inner join cdd_profiles p1 on d.kplet_1 = p1.id
                    inner join cdd_profiles p2 on d.kplet_2 = p2.id
                    inner join prok1402_baited_files bf on bf.id = df.file_id
                    inner join sources s on bf.source_id = s.id
                    inner join genomes g on g.id = s.genome_id
                    inner join weights w on w.genome_id = g.id
                    group by d.id
                    order by loci desc;"""

    _db = DbClass()
    _db.cmd = sql_cmd

    rows = _db.retrieve()

    return rows


