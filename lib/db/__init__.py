__author__ = 'hudaiber'

global connection
import MySQLdb as mdb
from warnings import filterwarnings
filterwarnings('ignore', category = mdb.Warning)
connection = mdb.connect(host='mysqldev21', user='hudaiber', db='PatternQuest', passwd='buP!est9')

# connection = mdb.connect(host='localhost', user='root', db='PatternQuest', passwd='')

# import sqlite3


# def create_mysql(db_name):
#     global connection
#     connection = sqlite3.connect(db_name)
#
# def create_mysql_in_memory():
#     global connection
#     connection = sqlite3.connect(":memory:")

def setup_cursor():
    try:
        cursor = connection.cursor()
        return cursor
    except ConnectionDoesNotExist:
        print "Database not configured"
        return None


class DbClass(object):
    def __init__(self, cmd=None):
        self.cursor = setup_cursor()
        self.cmd = cmd

        _cmd = """SET group_concat_max_len=100000000"""
        self.cursor.execute(_cmd)

    def execute(self):
        self.cursor.execute(self.cmd)

    @staticmethod
    def commit():
        connection.commit()

    def retrieve(self):

        if not self.cmd:
            raise ValueError("Provide valid SQL command")

        self.execute()
        return self.cursor.fetchall()

    def get_last_id(self):
        return self.cursor


def map_cdd2id():

    _db = DbClass()
    _db.cmd = """select code, id from cdd_profiles"""

    rows = _db.retrieve()

    return {row[0]: int(row[1]) for row in rows}


def map_id2cdd():
    _db = DbClass()
    _db.cmd = """select id, code from cdd_profiles"""

    rows = _db.retrieve()

    return {int(row[0]): row[1] for row in rows}


def map_id2cdd_clusters():
    _db = DbClass()
    _db.cmd = """select id, code from cdd_profiles"""

    rows = _db.retrieve()

    return {int(row[0]): row[1] for row in rows}


def map_id2cdd_cdd2id_cdd2def():
    _db = DbClass()
    _db.cmd = """select code, id, description from cdd_profiles"""

    rows = _db.retrieve()

    _cdd2id, _id2cdd, _cdd2def = {}, {}, {}

    for row in rows:
        (code, id, desctiption) = row

        _cdd2id[code] = id
        _id2cdd[id] = code
        _cdd2def[code] = desctiption

    return _id2cdd, _cdd2id, _cdd2def

