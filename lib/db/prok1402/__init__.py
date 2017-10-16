__author__ = 'hudaiber'

import sys
import os
sys.path.append('../')
from .. import DbClass


def neighborhoods_path():

    _path = os.path.join(gv.project_data_path, 'Bacteria/genes_and_flanks/win_10/raw_nbr_files/')
    if os.path.exists(_path):
        return _path

    _path = os.path.join('/panfs/pan1.be-md.ncbi.nlm.nih.gov/patternquest/Projects/NewSystems/data/Bacteria/genes_and_flanks/win_10/raw_nbr_files/')

    if os.path.exists(_path):
        return _path

    raise IOError("The neighborhoods path doesn't exist. Check the project data paths.")


def merged_neighborhoods_path():

    _path = os.path.join(gv.project_data_path, 'Bacteria/genes_and_flanks/win_10/merged/')
    if os.path.exists(_path):
        return _path

    _path = os.path.join('/panfs/pan1.be-md.ncbi.nlm.nih.gov/patternquest/Projects/NewSystems/data/Bacteria/genes_and_flanks/win_10/merged/')

    if os.path.exists(_path):
        return _path

    raise IOError("The neighborhoods path doesn't exist. Check the project data paths.")


