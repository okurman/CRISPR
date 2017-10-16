#!/home/hudaiber/env2.7/bin/python
__author__ = 'hudaiber'

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

# from lib.utils import BasicLocus, Locus
# import dm_tools as dt

# from lib.utils.prok1603 import parser

data_path = os.path.join(gv.project_data_path, 'cas4')
target_path = os.path.join(gv.project_data_path, 'cas4/figures/phylo_tree')

taxonomy_file = '/panfs/pan1.be-md.ncbi.nlm.nih.gov/prokdata/db/all1603/all1603.tax.tab'

genome2domain = {l.split('\t')[0]: l.split('\t')[4] for l in open(taxonomy_file)}
genome2phylum = {l.split('\t')[0]: l.split('\t')[5] for l in open(taxonomy_file)}
genome2class = {l.split('\t')[0]: l.split('\t')[6] for l in open(taxonomy_file)}
genome2Order = {l.split('\t')[0]: l.split('\t')[7] for l in open(taxonomy_file)}
genome2Family = {l.split('\t')[0]: l.split('\t')[8] for l in open(taxonomy_file)}
genome2Genus = {l.split('\t')[0]: l.split('\t')[9] for l in open(taxonomy_file)}
genome2Species = {l.split('\t')[0]:  l.split('\t')[10] for l in open(taxonomy_file)}

genome2tax_lineage = {l.split('\t')[0]:l.split('\t')[3] for l in open(taxonomy_file)}


def census_table():

    cas42locus = {}
    genome2cas4 = {}
    loci_file = os.path.join(data_path, 'Islands_ID.ann_clust')
    loci = t.parse_crispr_loci(loci_file)

    census_outf = open(os.path.join(data_path, 'cas4_census.txt'), "w")
    # columns = ["GI", "profiles", "Gene name(s)", "Status", "Complete", "Type", "Organism", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    columns = ["GI", "profiles", "Gene name(s)", "Status", "Complete", "Type", "Organism", "Total Cas4 in genome", "Kingdom"]
    census_outf.write("#"+"\t".join(columns)+"\n")

    # annotation_outf = open(os.path.join(target_path, 'cas4_annotation.txt'), "w")
    # annotation_outf.write("\t".join(["GI","annotation"]) + "\n")

    for locus in loci:
        cas42locus.update({gene:locus for gene in locus.cas4_genes})
        [t.update_dictionary_set(genome2cas4, locus.organism, _cas4.gid) for _cas4 in locus.cas4_genes]

    for cas4_gene, locus in sorted(cas42locus.items(), key=lambda x: x[0].gid):

        _gi = cas4_gene.gid
        _profiles = cas4_gene.cogid
        _gene_names = ",".join(set(cas4_gene.gene_name.split(',')))
        _status = locus.status
        _complete = "complete" if locus.complete else "partial"
        _type = locus.type if ";" not in locus.type else locus.type.split(";")[0]
        _org = locus.organism

        _domain = genome2domain[_org]
        # _phylum = genome2phylum[_org]
        # _class = genome2class[_org]
        # _order = genome2Order[_org]
        # _family = genome2Family[_org]
        # _genus = genome2Genus[_org]
        # _species = genome2Species[_org].strip()
        # _lineage = genome2tax_lineage[_org]
        # terms = [_gi, _profiles, _gene_names, _status, _complete, _type, _org, _domain, _phylum, _class, _order, _family, _genus, _species]
        terms = [_gi, _profiles, _gene_names, _status, _complete, _type, _org, str(len(genome2cas4[_org])), _domain]
        census_outf.write("\t".join(terms)+"\n")

        #------------------------------------------------
        # For annotation: Known, Unknow, fusion
        # if "uvrd" in _gene_names.lower():
        #     _annotation = "fusion"
        #
        # elif locus.complete and locus.status == "Known":
        #     _annotation = _type
        # else:
        #     _annotation = "unknown"
        # ------------------------------------------------


        # ------------------------------------------------
        # For annotation: Known, Unknow, fusion
        # if "uvrd" in _gene_names.lower():
        #     _annotation = "A_fusion"
        # elif locus.status == "Known":
        #     _annotation = "C_"+_type.replace("part ", "")
        # else:
        #     _annotation = "B_unknown"
        # ------------------------------------------------


        # annotation_outf.write(_gi+"\t"+_annotation+"\n")


if __name__ == "__main__":

    census_table()