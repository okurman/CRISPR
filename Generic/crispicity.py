#! /usr/bin/env python
__author__ = 'hudaiber'

import sys
sys.path.append('../')
import os
if sys.platform == 'darwin':
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))
elif sys.platform == 'linux2':
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
    sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))
import global_variables as gv

import dm_tools as dt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import lib.utils.tools as t


class ProfileInOrganismCount:
    def __init__(self, org, profile):
        self.org = org
        self.profile = profile
        self.in_crispr = 0
        self.outside = 0


def load_maps(fname, cas1402_organisms):

    gid2profiles, profile2orgs2gis = {}, {}

    with open(fname) as inf:
        for l in inf:
            terms = l.split(',')
            gid = terms[0]
            org = terms[1]
            profile = terms[6]

            # if gid in cas1402_gis or org not in cas1402_organisms:
            if org not in cas1402_organisms:
                continue

            if gid not in gid2profiles:
                gid2profiles[gid] = set([profile])
            else:
                gid2profiles[gid].update([profile])

            if profile not in profile2orgs2gis:
                profile2orgs2gis[profile] = {}

            if org not in profile2orgs2gis[profile]:
                profile2orgs2gis[profile][org] = set([gid])
            else:
                profile2orgs2gis[profile][org].update([gid])

    return gid2profiles, profile2orgs2gis


def load_maps_simple(fname, cas1402_gis):

    profile2orgs2gis = {}
    print "fname:", fname
    with open(fname) as inf:
        for l in inf:
            terms = l.split(',')
            gid = terms[0]

            if gid in cas1402_gis:
                continue

            org = terms[1]
            profile = terms[6]

            if profile not in profile2orgs2gis:
                profile2orgs2gis[profile] = {}

            if org not in profile2orgs2gis[profile]:
                profile2orgs2gis[profile][org] = set([gid])
            else:
                profile2orgs2gis[profile][org].update([gid])

    return profile2orgs2gis


def load_gi2cas_dom():

    fname = os.path.join(gv.project_data_path, 'cas1402/cas1402.gi2dom.tab')

    return {l.split()[0]:l.split()[1] for l in open(fname).readlines()}


def calculate_profile_based_crispricity_old(cas1402_loci, cas1402_gis, cas1402_organisms, prok1402_path_file):

    print "Loding global maps"
    global_gid2profiles, global_profile2orgs2gis = load_maps(prok1402_path_file, cas1402_gis, cas1402_organisms)
    print "Loading weights"
    gnm2weight = t.map_genome2weight()

    print "Counting in CRISPR loci"

    profile2orgs2count = {}

    for locus in cas1402_loci:
        for gene in locus:

            for _cogid in gene.cogid.split(','):

                if _cogid not in profile2orgs2count:
                    profile2orgs2count[_cogid] = {}

                if gene.organism not in profile2orgs2count[_cogid]:

                    _orgObj = ProfileInOrganismCount(gene.organism, _cogid)

                    outside_count = 0

                    if gene.gid in global_gid2profiles:
                        for _profile in global_gid2profiles[gene.gid]:
                            for _org in global_profile2orgs2gis[_profile]:
                                outside_count += 1 if global_profile2orgs2gis[_profile][_org].difference(cas1402_gis) else 0

                    _orgObj.outside = outside_count
                    profile2orgs2count[_cogid][gene.organism] = _orgObj

                profile2orgs2count[_cogid][gene.organism].in_crispr = 1

    in_crispr_all   = []
    crispricity_all = []
    profiles_all    = []

    print "Writing to files"
    with open('crispricity_profiles.tab', 'w') as outf_profiles:
        with open('crispricity_gis.tab', 'w') as outf_gis:

            for profile in profile2orgs2count:
                in_crispr = 0
                everywhere = 0

                for org in profile2orgs2count[profile]:
                    _org = profile2orgs2count[profile][org]
                    in_crispr +=  _org.in_crispr * gnm2weight[org]
                    everywhere += (_org.in_crispr + _org.outside) * gnm2weight[org]

                crispricity = in_crispr / everywhere

                in_crispr_all.append(in_crispr)
                crispricity_all.append(crispricity)
                profiles_all.append(profile)
                if profile.isdigit():
                    outf_gis.write("%s\t%f\t%f\n"%(profile, in_crispr, crispricity))
                else:
                    outf_profiles.write("%s\t%f\t%f\n"%(profile, in_crispr, crispricity))


    in_crispr_all   = np.asarray(in_crispr_all)
    in_crispr_all   = np.log10(in_crispr_all)
    crispricity_all = np.asarray(crispricity_all)
    # crispricity_all = np.log(crispricity_all)

    plt.ioff()
    fig, ax = plt.subplots()
    ax.scatter(in_crispr_all, crispricity_all,s=1)

    plt.xlabel("Effective orcurrence in CRISPR loci (log10)")
    plt.ylabel("X-axis / Effective occurrences")

    plt.savefig('crispricity.png')




def calculate_profile_based_crispricity(cas1402_loci, cas1402_gis, cas1402_organisms, prok1402_path_file):

    print "Loding global maps"
    global_profile2orgs2gis = load_maps_simple(prok1402_path_file, cas1402_gis)
    print "Loading weights"
    gnm2weight = t.map_genome2weight()
    print "Loading CDD definitions"
    profile2def = t.map_cdd_profile2def()

    print "Counting in CRISPR loci"

    profile2orgs2obj = {}

    for locus in cas1402_loci:
        for gene in locus:

            for _cogid in gene.cogid.split():

                if _cogid not in profile2orgs2obj:
                    profile2orgs2obj[_cogid] = {}
                    for _org in cas1402_organisms:
                        _orgObj = ProfileInOrganismCount(_org, _cogid)

                        if _cogid in global_profile2orgs2gis:
                            _orgObj.outside = len(global_profile2orgs2gis[_cogid][_org]) \
                                              if _org in global_profile2orgs2gis[_cogid] \
                                              else 0
                        else:
                            _orgObj.outside = 0

                        profile2orgs2obj[_cogid][_org] = _orgObj

                profile2orgs2obj[_cogid][gene.organism].in_crispr += 1

    out_file = os.path.join(gv.project_data_path, 'cas1402/crispricity_count.tab')

    in_crispr_all   = []
    crispricity_all = []
    profiles_all    = []

    print "Writing to file:", out_file
    with open(out_file, 'w') as outf:

        outf.write("Profile\tOccurrence in CRISPR loci\tCrispricity\tDefinition\n")

        for profile in profile2orgs2obj:
            in_crispr = 0
            everywhere = 0

            # for org in profile2orgs2obj[profile]:
            #     _org = profile2orgs2obj[profile][org]
            #     in_crispr +=  _org.in_crispr * gnm2weight[org]
            #     everywhere += (_org.in_crispr + _org.outside) * gnm2weight[org]

            for org in profile2orgs2obj[profile]:
                _org = profile2orgs2obj[profile][org]
                in_crispr += _org.in_crispr
                everywhere += (_org.in_crispr + _org.outside)

            crispricity = in_crispr / everywhere

            in_crispr_all.append(in_crispr)
            crispricity_all.append(crispricity)
            profiles_all.append(profile)
            outf.write("%s\t%f\t%f\t%s\n"%(profile, in_crispr, crispricity, profile2def[profile]))

    in_crispr_all   = np.asarray(in_crispr_all)
    in_crispr_all   = np.log10(in_crispr_all)
    crispricity_all = np.asarray(crispricity_all)
    # crispricity_all = np.log(crispricity_all)

    plt.ioff()
    fig, ax = plt.subplots()
    ax.scatter(in_crispr_all, crispricity_all, s=1)

    plt.xlabel("Effective orcurrence in CRISPR loci (log10)")
    plt.ylabel("X-axis / Effective occurrences")

    # fig.savefig('first.png')
    plt.savefig('first_count.png')

    # for i, profile in enumerate(profiles_all):
    #     ax.annotate(profile, (in_crispr_all[i], crispricity_all[i]))
    # fig.savefig('second.png')
    # plt.savefig('second.png')


def dull_gene_name():

    cas_gene_names = [l.strip() for l in open(os.path.join(gv.project_data_path,'cas1402/all_gene_names.txt'))]

    gene_name2gids = { gene:set() for gene in cas_gene_names }

    cnt = 0
    with open(os.path.join(gv.project_data_path,'cas1402/cas1402.arrisl.lst')) as inf:

        for in_line in inf:

            if in_line.startswith("==="):
                continue

            parts = in_line.strip().split('\t')
            if len(parts) < 9:
                continue

            _gene = parts[8]

            if _gene in cas_gene_names:

                gene_name2gids[_gene].update([parts[0]])

    cdd_gid2profiles = t.map_gid2cdd()

    cas_gene2profile = { gene:{} for gene in cas_gene_names }

    for _cas_gene in cas_gene_names:
        for _gid in gene_name2gids[_cas_gene]:
            if not _gid in cdd_gid2profiles:
                # t.update_dictionary(cas_gene2profile[_cas_gene], "NA", 1)
                continue
            for _profile in cdd_gid2profiles[_gid].split():

                t.update_dictionary(cas_gene2profile[_cas_gene], _profile, 1)

    work_dir = os.path.join(gv.project_data_path,'cas1402/crispricity/')

    with open(os.path.join(work_dir, 'gene_name2profiles.txt'), 'w') as outf:
        for _gene_name in cas_gene_names:
            for _profile in cas_gene2profile[_gene_name]:
                outf.write("%s\t%s\t%d\n" % (_gene_name, _profile, cas_gene2profile[_gene_name][_profile]))

    cas_related_profiles = set([_profile for _gene in cas_gene_names for _profile in cas_gene2profile[_gene].keys()])

    cr_occurrence = []
    cr_crispricity = []
    ncr_occurrence = []
    ncr_crispricity = []

    for l in open(os.path.join(work_dir, 'crispricity.tab')).readlines()[1:]:
        if not l:
            continue
        parts = l.split('\t')

        if parts[0] in cas_related_profiles:
            cr_occurrence.append(parts[1])
            cr_crispricity.append(parts[2])
        else:
            ncr_occurrence.append(parts[1])
            ncr_crispricity.append(parts[2])

    cr_occurrence = np.asarray(cr_occurrence, dtype=np.float)
    cr_occurrence = np.log(cr_occurrence)
    cr_crispricity = np.asarray(cr_crispricity)

    ncr_occurrence = np.asarray(ncr_occurrence, dtype=np.float)
    ncr_occurrence = np.log(ncr_occurrence)
    ncr_crispricity = np.asarray(ncr_crispricity)

    plt.ioff()
    fig, ax = plt.subplots()
    ax.scatter(cr_occurrence, cr_crispricity, color='r', s=1, label="Cas related")
    ax.scatter(ncr_occurrence, ncr_crispricity, color='b', s=1, label="Not Cas related")

    ax.axvline(1.6 ,color='g', linewidth=0.5)
    ax.axhline(0.5 ,color='g', linewidth=0.5)

    plt.xlabel("Effective orcurrence in CRISPR loci (log)")
    plt.ylabel("Crispricity")

    plt.legend(loc="upper left", fontsize=7)

    plt.savefig(os.path.join(work_dir, 'crispricity_log.png'))


if __name__=='__main__':

    # prok1402_path_file = '/mnt/storage/data/CDD/all_Prok1402.ccp.csv'
    prok1402_path_file = os.path.join(os.path.expanduser('~'),'data/CDD/all_Prok1402.ccp.csv')

    print "Loading map_gid2cdd"
    gi2annotation = t.map_gid2cdd()

    print "Loading CRISPR loci"
    # cas1402_loci_path = os.path.join(gv.project_data_path,'cas1402/files/')
    cas1402_loci_path = os.path.join(gv.project_data_path,'cas1402/files/')
    cas1402_loci = [dt.get_pty_file_generic(os.path.join(cas1402_loci_path, f), annotation_map = gi2annotation) for f in os.listdir(cas1402_loci_path)]

    cas1402_gis = set([gene.gid for locus in cas1402_loci for gene in locus])
    cas1402_organisms = set([locus[0].organism for locus in cas1402_loci])

    calculate_profile_based_crispricity(cas1402_loci, cas1402_gis, cas1402_organisms, prok1402_path_file)

    # all_orgs = set()
    # for profile, org2gis in global_profile2orgs2gis.items():
    #     if 'Yersinia_pseudotuberculosis_PB1__uid59153' in org2gis:
    #         print "Yersinia_pseudotuberculosis_PB1__uid59153"
    #     for _org in org2gis.keys():
    #         all_orgs.update([_org])