#!/usr/bin/env python

import os
import sys
sys.path.append('/Users/hudaiber/Projects/lib/BioPy/')
import dm_tools
sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
import global_variables as gv
import lib.utils.tools as t
import time
import pickle

nbrhds_path = '/Users/hudaiber/Projects/NewSystems/data/Archea/genes_and_flanks/win_40/pty/'
inter_dist = 21

def build_src2gids_map():

    src_to_gids = {}
    cnt=0
    t=time.time()
    for f in [l for l in os.listdir(nbrhds_path) if l.endswith('_annot.pty')]:
        tmp_l = open(os.path.join(nbrhds_path, f)).readline()
        cur_gid = f.split('_')[0]
        cur_src = tmp_l.split()[4]
        if src_to_gids.has_key(cur_src):
            src_to_gids[cur_src].append(cur_gid)
        else:
            src_to_gids[cur_src] = [cur_gid]
        cnt+=1
        if cnt%10000==0:
            print 'Processed: %d, Time: %f' % (cnt, time.time()-t)
            t=time.time()
    pickle.dump(src_to_gids, open('src_to_gids.p','w'))
    src_to_gids = pickle.load(open('src_to_gids.p'))
    outf=open('src_to_gids.txt','w')

    for k,v in src_to_gids.items():
        outf.write('%d\t%s\t%s\n' % (len(v), k, ' '.join(v)))


def write_merged_block(block, background_genes, fname, source_path):

    subject_files = ['%s_annot.pty'%background_genes[i].gid for i in block]
    subject_files = [ os.path.join(source_path, f) for f in subject_files ]
    lines = []
    for f in subject_files:
        lines += open(f).readlines()
    lines = set(lines)
    lines = sorted(lines)

    outf = open(fname, 'w')
    [outf.write(l) for l in lines]


def merge_neighborhoods(src_to_gid_file, dest_path, source_path):

    """If the neighborhood center (subject gene) of two neighborhoods are located closer than 20 genes,
    then merge them into one"""

    src_to_org = t.map_src2org()

    for l in open(src_to_gid_file):
        src, gids = l.split('\t')[1], l.split('\t')[2].split()
        org = src_to_org[src]
        subject_genes = dm_tools.gid_to_gene(os.path.join(gv.pty_genomes_path, org, "%s.pty"%src), gids)
        subject_genes.sort()

        background_genes = dm_tools.get_pty_file(os.path.join(gv.pty_genomes_path, org, '%s.pty'%src))
        background_genes.sort()
        background_genes_gids = [g.gid for g in background_genes]

        blocks = []
        indices = [background_genes_gids.index(g.gid) for g in subject_genes]
        indices.sort()

        cur_block=[]
        for ind in indices:
            if not cur_block:
                cur_block.append(ind)
                continue
            diff = ind - cur_block[-1]
            assert(diff>0)
            if diff < inter_dist:
                cur_block.append(ind)
            else:
                blocks.append(cur_block)
                cur_block=[ind]
        blocks.append(cur_block)

        cnt=1
        for block in blocks:
            write_merged_block(block, background_genes, os.path.join(dest_path, '%s_%d.pty'%(src, cnt)), source_path)
            cnt+=1

        assert (sum([len(b) for b in blocks]) == len(indices))




if __name__=='__main__':

    # gid_list_file = sys.argv[1]
    # out_dest = sys.argv[2]

    # build_src2gids_map()

    merge_neighborhoods('src_to_gids.txt', '../../../data/Archea/genes_and_flanks/win_40/merged/', '../../../data/Archea/genes_and_flanks/win_40/pty/')

