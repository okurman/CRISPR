__author__ = 'hudaiber'

import sys
import os
from .. import DbClass

if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')
import global_variables as gv


def neighborhoods_path():

    _path = os.path.join(gv.project_data_path, 'CRISPR/Cas1/genes_and_flanks/')
    if os.path.exists(_path):
        return _path
    else:
        print _path
        raise IOError("The neighborhoods path doesn't exist. Check the project data paths.")

def map_file2org():

    _fname = os.path.join(gv.project_data_path, 'cas1402/file2org.txt')
    out_map = {}
    for l in open(_fname).readlines():
        terms = l.strip().split()
        out_map[terms[0]] = terms[1]

    return out_map


def map_profile2id(prefix):
    _db = DbClass()
    _db.cmd = """select code, id from %s_profiles""" % prefix
    rows = _db.retrieve()
    return {row[0]: int(row[1]) for row in rows}


def map_file_id2name(prefix):
    _db = DbClass()
    _db.cmd = """select id, name from %s_files""" % prefix
    return {str(l[0]): l[1] for l in _db.retrieve()}


def map_file_name2id(prefix):
    _db = DbClass()
    _db.cmd = """ select name, id from  %s_files""" % prefix
    rows = _db.retrieve()
    return {row[0]: row[1] for row in rows}


def map_profiles_id2code_code2def(prefix):
    _db = DbClass()
    _db.cmd = """select id, code, description from cdd_profiles"""
    rows = _db.retrieve()

    _id2code = dict()
    _code2def = dict()

    for row in rows:
        (_id, _code, _def) = row
        _id2code[_id] = _code
        _code2def[_code] = _def

    return _id2code, _code2def