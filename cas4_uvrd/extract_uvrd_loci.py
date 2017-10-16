#!/home/hudaiber/env2.7/bin/python
__author__ = 'hudaiber'

import os
import sys
sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))

import global_variables as gv
sys.path.append(gv.project_code_path)

from lib.utils import BasicLocus
import dm_tools as dt
from lib.utils import tools as t

def extract_UvrD_baited_loci():

    print "Loading annotations"
    gi2cdd = {}

    with open(annotation_file) as inf:
        for l in inf:
            terms = l.split(',')
            gi = terms[0]
            profile = terms[6]

            if gi not in gi2cdd:
                gi2cdd[terms[0]] = [terms[6]]
            else:
                gi2cdd[terms[0]].append(terms[6])

    print "Extracting loci"

    cnt =0
    for l in open(os.path.join(work_dir, 'selected_annotations_pty.txt')):

        terms = l.split()
        org = terms[3]
        source = terms[4]
        gid = terms[6]

        pty_file = os.path.join(pty_dir, org, "%s.pty"%source)

        _genes = dt.get_prok1603_pty_file(pty_file, annotation_map=gi2cdd)
        _genes.sort()

        # if len(_genes) == 0:
        #     print "Len:", len(_genes), pty_file
        #     sys.exit()

        locus_file = os.path.join(work_dir, 'files/%s.pty' % gid)

        for i in range(len(_genes)):
            if gid == _genes[i].gid:
                dt.write_pty_file(_genes[i-10 if i>10 else i: i+10], locus_file)

        cnt += 1
        if not cnt % 1000:
            print cnt


def merge_loci(source_path, dest_path):

    profile2gene = t.map_cdd_profile2gene_name()

    print "Loading loci"
    loci = [BasicLocus(os.path.join(source_path, f), profile2gene=profile2gene) for f in os.listdir(source_path)]

    print "Merging loci"

    cnt = 1
    while True:

        merged_out = [0]*len(loci)
        new_loci = []

        for i in range(len(loci)-1):

            if merged_out[i]:
                continue

            for j in range(i+1, len(loci)):

                if merged_out[i]:
                    continue

                if loci[i].overlaps(loci[j]):
                    loci[i].merge(loci[j])
                    merged_out[j] = 1

            new_loci.append(loci[i])

        if not merged_out[len(loci)-1]:
            new_loci.append(loci[-1])

        print "Iteration %d results. Old list: %d, New list: %d" % (cnt, len(loci), len(new_loci))
        cnt += 1

        if len(loci) == len(new_loci):
            loci = new_loci
            break

        loci = new_loci

    # with open(work_dir+'merged_loci_name.txt', 'wb') as outf:
    #
    #     [outf.write("%s\n"%_file) for locus in loci for _file in locus.merged_files]

    print "Writing merged loci to files"

    for locus in loci:

        fname = os.path.join(dest_path, locus.merged_base_files.pop())

        while os.path.exists(fname):
            print "Duplicate file name:", fname
            fname = os.path.join(dest_path, locus.merged_base_files.pop())

        header_lines=[]
        if locus.merged_base_files:
            header_line = "Merged files: %s" % ",".join(f for f in locus.merged_base_files)
            header_lines = [header_line]

        t.write_genes_to_pty(locus.genes, fname, header_lines=header_lines)





if __name__=='__main__':

    work_dir = os.path.join(gv.project_data_path, 'UvrD/prok1603/')
    pty_dir = os.path.join(os.path.expanduser('~'),'data/Prok1603/pty/')
    annotation_file = '/panfs/pan1.be-md.ncbi.nlm.nih.gov/prokdata/db_tmp/Prok1603/Annotation/Prok1603.ccp.csv'

    merge_loci(os.path.join(work_dir, 'files'), os.path.join(work_dir, 'merged_files'))
