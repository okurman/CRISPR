#!/home/hudaiber/env2.7/bin/python

import sys
import numpy as np

input_file = sys.argv[1]

values = np.loadtxt(input_file)
values.sort()

_size=np.size(values)

means = np.array([np.mean(values),
                  np.mean(values[:_size/4]),
                  np.mean(values[_size/4:_size*3/4]),
                  np.mean(values[_size*3/4:])])

print [ float(np.round(a, 2)) for a in means]