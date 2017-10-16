__author__ = 'hudaiber'


class CrisprMergingKplets2FsOutput(object):

    def __init__(self):
        self.file_summaries = None
        self.organisms = None
        self.crispr_type2files = None
        self.kplet2count_af = None
        self.kplet2count_bf = None
        self.initial_length = None

class CrisprReportingWgsInput(object):

    def __init__(self):
        self.xls_file_name = None
        self.file_summaries = None
        self.organisms = None
        self.profile_code2def = None
        self.crispr_type2files = None
        self.local_bf_kplet2count = None
        self.local_af_kplet2count = None
        self.initial_length = None
        self.wgs_global_profile_count = None

class CrisprReportingWgsSummaryInput(object):

    def __init__(self):
        self.worksheet = None
        self.cur_worksheet = None
        self.kplet_list = None
        self.ind = None
        self.local_bf_kplet2count = None
        self.local_af_kplet2count = None
        self.wgs_global_profile_count = None
        self.crispr_type_summary = None


class GenericMergingKplets2FsOutput(object):

    def __init__(self):
        self.file_summaries = None
        self.organisms = None


class GenericReportingInput(object):

    def __init__(self):
        self.xls_file_name = None
        self.file_summaries = None
        self.organisms = None
        self.weight = None
        self.profile_code2def = None
        self.local_bf_kplet2count = None
        self.local_af_kplet2count = None
