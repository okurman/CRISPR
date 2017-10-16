__author__ = 'Sanjarbek Hudaiberdiev'

from lib.db import DbClass


def map_profile2id(profile_list):

    _db = DbClass()
    if profile_list == 'all':
        sql_cmd = """select code, id from cdd_profiles"""
    else:
        sql_cmd = """select code, id from cdd_profiles where code in ('%s')"""
        sql_cmd = sql_cmd % "','".join(profile_list)

    _db.cmd = sql_cmd

    return {row[0]: int(row[1]) for row in _db.retrieve()}


def get_file_id(fname):
    sql_cmd = """select id from archea_win10_files where name='%s'"""
    sql_cmd = sql_cmd % fname
    cursor = setup_cursor()
    cursor.execute(sql_cmd)
    rows = cursor.fetchall()
    assert len(rows) <= 1

    if not rows:
        raise Exception("File not found in database:%s" % fname)

    return rows[0][0]


def org2src_src2files_map(files):

    _db = DbClass()

    _db.cmd = """select awf.name, s.name, g.name
                  from bacteria_win10_files awf
                  inner join sources s on awf.source_id=s.id
                  inner join genomes g on s.genome_id = g.id
                  where awf.name in ('%s')""" % "','".join(files)

    _org2src = {}
    _src2files = {}

    for row in _db.retrieve():
        parts = row
        [_file, _src, _org] = parts
        if _org in _org2src:
            _org2src[_org].update([_src])
        else:
            _org2src[_org] = set([_src])

        if _src in _src2files:
            _src2files[_src].update([_file])
        else:
            _src2files[_src] = set([_file])

    return _org2src, _src2files



def file2src_src2org_map(files):

    _db = DbClass()
    _db.cmd = """select awf.name, s.name, g.name
                  from bacteria_win10_files awf
                  inner join sources s on awf.source_id=s.id
                  inner join genomes g on s.genome_id = g.id
                  where awf.name in ('%s')""" % "','".join(files)

    _src2org = {}
    _file2src = {}
    for row in _db.retrieve():
        parts = row
        [_file, _src, _org] = parts
        _src2org[_src] = _org
        _file2src[_file] = _src

    return _file2src, _src2org



def map_file_id2name():

    _db = DbClass()
    _db.cmd = """select id, name from bacteria_win10_files"""
    return {str(l[0]): l[1] for l in _db.retrieve()}


def map_file_name2id():

    _db = DbClass()
    _db.cmd = """ select name, id from  bacteria_win10_files"""
    rows = _db.retrieve()
    return {row[0]: row[1] for row in rows}


def map_baited_file2id():

    _db = DbClass()
    _db.cmd = """ select name, id from prok1402_baited_files"""
    rows = _db.retrieve()

    file2id = {}
    id2file = {}

    for row in rows:
        file2id[row[0]] = row[1]
        id2file[row[1]] = row[0]

    return file2id, id2file


def map_source2id():

    _db = DbClass()
    _db.cmd = """ SELECT id, name FROM PatternQuest.sources; """
    rows = _db.retrieve()
    return {row[1]: row[0] for row in rows}



if __name__ == '__main__':
    pass
    # kplets_ranked = get_heavy_archaea_kplets()

    # for kplet in kplets_ranked[:100]:
    #     print kplet
