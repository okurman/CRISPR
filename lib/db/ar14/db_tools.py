__author__ = 'Sanjarbek Hudaiberdiev'

import MySQLdb as mdb
connection = mdb.connect(host='mysql-dev', user='hudaiber', db='PatternQuest', passwd='buP!est9')

def setup_cursor():
    try:
        cursor = connection.cursor()
        return cursor
    except ConnectionDoesNotExist:
        print "Database not configured"
        return None


def file_name2id():
    _cursor = setup_cursor()
    _sqlcmd = """ select id, name
                  from  archea_win10_files"""

    _cursor.execute(_sqlcmd)
    print { v:int(k) for k,v in _cursor.fetchall()}


def map_profile2id(profile_list):

    sql_cmd = """select code, id from cdd_profiles where code in ('%s')"""
    sql_cmd = sql_cmd%"','".join(profile_list)
    cursor = setup_cursor()
    cursor.execute(sql_cmd)
    return {row[0]: int(row[1]) for row in cursor.fetchall()}


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

    _sql_cmd = """select awf.name, s.name, g.name
                  from archea_win10_files awf
                  inner join sources s on awf.source_id=s.id
                  inner join genomes g on s.genome_id = g.id
                  where awf.name in ('%s')"""

    _sql_cmd %= "','".join(list(files))

    _cursor = setup_cursor()
    _cursor.execute(_sql_cmd)
    _org2src = {}
    _src2files = {}

    for row in _cursor.fetchall():
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

    _sql_cmd = """select awf.name, s.name, g.name
                  from archea_win10_files awf
                  inner join sources s on awf.source_id=s.id
                  inner join genomes g on s.genome_id = g.id
                  where awf.name in ('%s')"""

    _sql_cmd = _sql_cmd % "','".join(files)

    _cursor = setup_cursor()
    _cursor.execute(_sql_cmd)
    _src2org = {}
    _file2src = {}

    for row in _cursor.fetchall():
        parts = row
        [_file, _src, _org] = parts
        _src2org[_src] = _org
        _file2src[_file] = _src

    return _file2src, _src2org






if __name__ == '__main__':

    kplets_ranked = get_heavy_archaea_kplets()

    # for kplet in kplets_ranked[:100]:
    #     print kplet