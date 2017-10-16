__author__ = 'hudaiber'

import sys
sys.path.append('../..')
from lib.db import DbClass

_db = DbClass()

# for table_name in ['crispr_4plets_wgs', 'crispr_3plets_wgs', 'crispr_2plets_wgs']:
for table_name in ['crispr_5plets_wgs']:
    print "Starting for table:", table_name
    print
    print "Loading the duplicates"
    _db.cmd = """SELECT count(*) as cnt, kplet_id, file_id
                FROM PatternQuest.%s
                group by kplet_id, file_id
                having count(*)>1
                order by cnt desc""" % table_name

    table_rows = [row for row in _db.retrieve()]

    print("Number of duplicate entries:", len(table_rows))
    fname = "duplicates_%s.txt" % table_name
    with open(fname, 'w') as f:
        [f.write("%d\t%d\t%d\n"%(row)) for row in table_rows]
    print "Duplicates saved in file:", fname
    continue
    sys.exit()
    print("Starting to remove duplicates")

    total_cnt = 0
    chunk_size = 10000
    chunks = [table_rows[i: i+chunk_size] for i in range(0, len(table_rows), chunk_size)]

    for row in table_rows:

        (_cnt, _kplet_id, _file_id) = row

        assert _cnt > 1

        _db.cmd = """delete from %s where kplet_id=%d and file_id=%d""" % (table_name, _kplet_id, _file_id)
        _db.execute()
        _db.commit()

        _db.cmd = """insert into %s(kplet_id, file_id) values (%d, %d)""" % (table_name, _kplet_id, _file_id)
        _db.execute()
        _db.commit()

        total_cnt += 1
        if total_cnt % 10000 ==0:
            print total_cnt

    print total_cnt
    print "Finished:", table_name
    print
    print