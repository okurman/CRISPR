#! /usr/bin/env python
__author__ = 'hudaiber'

import sys
sys.path.append('../')
import os
import platform

if sys.platform == 'darwin':
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))
elif sys.platform == 'linux2':
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))

if platform.node()=='okurman':
	mount_point='/mnt/ncbi_home/'
	sys.path.append(os.path.join(mount_point,'Projects/lib/BioPy/'))
	sys.path.append(os.path.join(mount_point,'Projects/SystemFiles/'))
	work_dir = os.path.join(mount_point,'/mnt/panfs/Projects/NewSystems/data/cas4/in_situ/')


import global_variables as gv

