#!/home/hudaiber/env2.7/bin/python

import os
import sys
if sys.platform == 'darwin':
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))
elif sys.platform == 'linux2':
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))

import global_variables as gv
sys.path.append(gv.project_code_path)
from lib.utils import tools as t


if __name__=='__main__':

    loci_file =             '/home/hudaiber/Projects/NewSystems/data/cas4/Islands_ID.ann_clust'
    tree_file =             '/home/hudaiber/Projects/NewSystems/data/cas4/trees/cas4_muscle.tre'
    tree_leaves_file =      '/home/hudaiber/Projects/NewSystems/data/cas4/trees/cas4_muscle_leaves.tre'

    _dir = tree_annotation_file = '/home/hudaiber/Projects/NewSystems/data/cas4/trees/'

    loci = t.parse_crispr_loci(loci_file)
    tree = open(tree_file).readline()

    gene2annot = {}
    gene2locus = {}

    for locus in loci:

        UvrDs = len([g for g in locus.genes if "uvrd" in g.gene_name.lower()])

        for i in range(len(locus.genes)):

            gene = locus.genes[i]

            _gene_name = gene.gene_name.lower()

            if "cas4" in _gene_name:

                gene2locus[gene.gid] = locus

                if "uvrd" in _gene_name:
                    gene2annot[gene.gid] = 'fusion'
                    continue

                _before = "" if i == 0 else locus.genes[i - 1].gene_name
                _after = "" if i == len(locus.genes) - 1 else locus.genes[i + 1].gene_name

                if "uvrd" in _before.lower() or "uvrd" in _after.lower():
                    gene2annot[gene.gid] = 'tandem'
                    continue

                if UvrDs:
                    gene2annot[gene.gid] = 'in_vicinity'
                    continue

                gene2annot[gene.gid] = 'else'

    # fname = os.path.join(_dir, 'cas4_tree_uvrd_annotations.txt')
    # with open(fname, 'w') as outf:
    #     outf.write("Taxa\tUvrD\n")
    #     for leaf in open(tree_leaves_file).readlines():
    #         leaf = leaf.strip()
    #         _locus = gene2locus[leaf]
    #         outf.write("%s\t%s\n" % (leaf, gene2annot[leaf]))
    #
    # fname = os.path.join(_dir, 'cas4_tree_status_annotations.txt')
    # with open(fname, 'w') as outf:
    #     outf.write("Taxa\tstatus\n")
    #     for leaf in open(tree_leaves_file).readlines():
    #         leaf = leaf.strip()
    #         _locus = gene2locus[leaf]
    #         outf.write("%s\t%s\n" % (leaf, _locus.status))
    #
    # fname = os.path.join(_dir, 'cas4_tree_complete_annotations.txt')
    # with open(fname, 'w') as outf:
    #     outf.write("Taxa\tcomplete\n")
    #     for leaf in open(tree_leaves_file).readlines():
    #         leaf = leaf.strip()
    #         _locus = gene2locus[leaf]
    #         outf.write("%s\t%d\n" % (leaf, int(_locus.complete)))

    fname = os.path.join(_dir, 'cas4_tree_type_annotations.txt')
    with open(fname, 'w') as outf:
        outf.write("Taxa\tType\n")
        for leaf in open(tree_leaves_file).readlines():
            leaf = leaf.strip()
            _locus = gene2locus[leaf]

            if "part" in _locus.type or "unk" in _locus.type:
                _type = "A_incomplete"
            else:
                _type = _locus.type.split(";")[0].replace("-", "_")
            outf.write("%s\t%s\n" % (leaf, _type))
