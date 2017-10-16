__author__ = 'hudaiber'

import sys
import os
sys.path.append('../')
from .. import DbClass

if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')
import global_variables as gv

def neighborhoods_path():

    _path = os.path.join(gv.project_data_path, 'Archea/genes_and_flanks/win_10/pty/')
    if os.path.exists(_path):
        return _path

    _path = os.path.join('/panfs/pan1.be-md.ncbi.nlm.nih.gov/patternquest/Projects/NewSystems/data/Archea/genes_and_flanks/win_10/pty/')

    if os.path.exists(_path):
        return _path

    raise IOError("The neighborhoods path doesn't exist. Check the project data paths.")

def merged_neighborhoods_path():

    _path = os.path.join(gv.project_data_path, 'Archea/genes_and_flanks/win_10/merged/')
    if os.path.exists(_path):
        return _path

    _path = os.path.join('/panfs/pan1.be-md.ncbi.nlm.nih.gov/patternquest/Projects/NewSystems/data/Archea/genes_and_flanks/win_10/merged/')

    if os.path.exists(_path):
        return _path

    raise IOError("The neighborhoods path doesn't exist. Check the project data paths.")


def map_file_id2name():
    _db = DbClass()
    _db.cmd = """select id, name from archea_win10_files"""

    return {str(l[0]): l[1] for l in _db.retrieve()}


def map_name2file_id():
    _db = DbClass()
    _db.cmd = """select name, id from archea_win10_files"""

    return {str(l[0]): l[1] for l in _db.retrieve()}

