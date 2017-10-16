__author__ = 'hudaiber'
import os
import sys
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')
import global_variables as gv
# sys.path.append(gv.project_code_path)
import dm_tools as dt

cas_type_file = os.path.join(gv.project_data_path, 'cas1402/cas1402.type.tab')
_gi2castype =  {l.strip().split('\t')[0]: l.strip().split('\t')[1].split(';') for l in open(cas_type_file).readlines()}

files_path = os.path.join(gv.project_data_path, 'cas1402/files')

file2type = {}

for f in os.listdir(files_path):
    _genes = dt.get_pty_file_generic(os.path.join(files_path, f))
    for _gene in _genes:
        if _gene.gid in _gi2castype:
            if f in file2type:
                file2type[f] += _gi2castype[_gene.gid]
            else:
                file2type[f] = _gi2castype[_gene.gid]

print len(file2type)
print

outf = os.path.join(gv.project_data_path, 'cas1402/file2type.tab')

with open(outf,'w') as of:
    for file, cas_type in file2type.items():
        of.write("%s\t%s\n" % (file, ";".join(cas_type)))
