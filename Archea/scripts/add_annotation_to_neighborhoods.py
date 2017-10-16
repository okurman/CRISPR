#!/usr/bin/env python

import os
import sys
import time


def load_map():
	annot_file = '/Users/hudaiber/data/CDD/all_Prok1402.ccp.csv'
	out_map={}
	with open(annot_file) as f:
		for l in f:
			parts = l.split(',')
			if out_map.has_key(parts[0]):
				out_map[parts[0]].append(parts[6])
			else:
				out_map[parts[0]] = [parts[6]]
	return out_map


def load_arcog_map():
	annot_file = '/Users/hudaiber/data/Archea/arCOG/ar14.arCOG.csv'
	out_map={}
	with open(annot_file) as f:
		for l in f:
			if 'arCOG' not in l:
				continue
			parts = l.split(',')
			if out_map.has_key(parts[0]):
				if parts[6] not in out_map[parts[0]]:
					out_map[parts[0]].append(parts[6])
			else:
				out_map[parts[0]] = [parts[6]]
	return out_map


if __name__=='__main__':

	gid_list_file = sys.argv[1]
	out_dest = sys.argv[2]

	print 'Loading map'
	annot_map = load_map()
	print 'Dictionary loaded'
	arcog_map = load_arcog_map()
	print 'ArCOG Dictionary loaded'
	print 'Starting annotation adding.'

	with open(gid_list_file) as f:
		for gid_file in f:
			with open(os.path.join(out_dest, '%s_annot%s'%(os.path.splitext(gid_file.strip()))),'w') as outf:
				for l in open(os.path.join(out_dest, gid_file.strip())):
					cur_gid = l.split()[0]
					
					if cur_gid in arcog_map:
						annotation = ' '.join(arcog_map[cur_gid])
					elif cur_gid in annot_map:
						annotation = ' '.join(annot_map[cur_gid])
					else:
						annotation = None
					if annotation:
						outf.write(l.strip() +'\t'+annotation+'\n')
					else:
						outf.write(l)
	print 'Done'