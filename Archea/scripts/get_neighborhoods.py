import sys
import subprocess as sp
import os

pty_dir='/Users/hudaiber/data/Pty/genomes/'
annotation_file = '/Users/hudaiber/Projects/NewSystems/data/Archea/arCOG/selected_annotations.csv'

dest_dir = '/Users/hudaiber/Projects/NewSystems/data/Archea/genes_and_flanks/win_40/pty/'
flank = 20
for l in open(annotation_file):
	parts = l.split(',')
	gi = parts[0]
	gnm = parts[1]
	
	for pty_f in os.listdir(os.path.join(pty_dir, gnm)):
		src_file = os.path.join(os.path.join(pty_dir, gnm, pty_f))
		if gi in [l.split()[0] for l in open(src_file)]:
			dest_file = dest_dir+'%s.pty'%gi
			cmd = ("grep", "-B", str(flank), "-A", str(flank), gi, src_file, ">", dest_file)
			sp.call("%s %s %s %s %s %s %s %s %s"%cmd, shell=True)