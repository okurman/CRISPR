__author__ = 'Sanjarbek Hudaiberdiev'

import os
import sys
import configparser
from collections import defaultdict
import re

###############################################################
config_file = os.path.join(os.path.expanduser('~'), 'paths.cfg')
cfg = configparser.ConfigParser()
cfg.read(config_file)

code_path = cfg.get('NewSystems', 'code_path')
data_path = cfg.get('NewSystems', 'data_path')
sys.path.append(code_path)

###############################################################


HHPRED_PROB_THR = 80
HHPRED_EVALUE_THR = 1e-3


class Hit(object):
    """docstring for Hit"""

    def __init__(self, line):

        self.line = line
        parts = line.strip().split()

        self.profile = parts[1]
        self.prob = float(parts[2])
        self.evalue = float(parts[3])
        self.score = float(parts[5])

        q_coords = parts[8]
        self.q_from = int(q_coords.split('-')[0])
        self.q_to = int(q_coords.split('-')[1])
        t_coords = parts[9]
        if "(" in t_coords:
            t_coords = t_coords[:t_coords.index("(")]

        self.t_from = int(t_coords.split('-')[0])
        self.t_to = int(t_coords.split('-')[1])

        self.t_len = int(re.search(r"\(([0-9]+)\)",line).groups()[0])


    def __str__(self):
        return "%s\t%f\t%s\t%d\t%d" % (self.profile, self.prob, self.evalue, self.q_from, self.q_to)

    def overlaps(self, other):
        dist1 = other.q_from - self.q_to
        dist2 = self.q_from - other.q_to
        return True if dist1 * dist2 >= 0 else False


def hhsearch_parse(hhr_file, prob_thr=HHPRED_PROB_THR, eval_thr=HHPRED_EVALUE_THR, query_cov_thr=0, target_cov_thr=0):

    lines = open(hhr_file).readlines()
    query_length = int(lines[1].split()[1])

    lines = lines[9:]

    ind = lines.index("\n")
    lines = lines[:ind]

    hits = []

    for line in lines:
        _hit = Hit(line)

        query_cov = (_hit.q_to - _hit.q_from) / float(query_length)
        target_cov = (_hit.t_to - _hit.t_from) / float(_hit.t_len)

        assert query_cov >= 0 and target_cov >= 0

        if _hit.prob < prob_thr or _hit.evalue > eval_thr or query_cov < query_cov_thr or target_cov < target_cov_thr:
            break

        hits.append(_hit)

    return hits
