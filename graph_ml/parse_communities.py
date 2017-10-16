
import os
import sys
import configparser
from collections import defaultdict

###############################################################
config_file = os.path.join(os.path.expanduser('~'),'paths.cfg')
cfg=configparser.ConfigParser()
cfg.read(config_file)

code_path = cfg.get('NewSystems','code_path')
data_path = cfg.get('NewSystems','data_path')
sys.path.append(code_path)

###############################################################

import networkx as nx




if __name__ == "__main__":
    pass