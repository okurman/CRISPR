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
        self.in_locus = 0
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



def calculate_profile_based_baiticity(bacteria_loci, loci_gis,
                                      loci_organisms,
                                      prok1402_path_file,
                                      bait_profiles,
                                      filter_threshold,
                                      save_path):

    print "Loding global maps"
    global_profile2orgs2gis = load_maps_simple(prok1402_path_file, loci_gis)
    print "Loading weights"
    gnm2weight = t.map_genome2weight()
    print "Loading CDD definitions"
    profile2def = t.map_cdd_profile2def()

    print "Counting in loci"

    profile2orgs2obj = {}

    gi_checklist = set()

    for locus in bacteria_loci:
        for gene in locus:

            if gene.gid in gi_checklist:
                continue

            for _cogid in gene.cogid.split():

                if _cogid in bait_profiles:
                    continue

                if _cogid not in profile2orgs2obj:
                    profile2orgs2obj[_cogid] = {}
                    for _org in loci_organisms:
                        _orgObj = ProfileInOrganismCount(_org, _cogid)

                        if _cogid in global_profile2orgs2gis:
                            _orgObj.outside = len(global_profile2orgs2gis[_cogid][_org]) \
                                              if _org in global_profile2orgs2gis[_cogid] \
                                              else 0
                        else:
                            _orgObj.outside = 0

                        profile2orgs2obj[_cogid][_org] = _orgObj

                profile2orgs2obj[_cogid][gene.organism].in_locus += 1

            gi_checklist.update([gene.gid])

    out_file = os.path.join(save_path, 'baiticity.tab')

    profiles = []

    in_loci_count = []
    baiticity_count = []

    in_loci_weight = []
    baiticity_weight = []

    rare_profiles_file = open(os.path.join(save_path, 'rare_profiles.tab'), 'w')
    rare_profiles_file.write("Profile\tOccurence everywhere\tOccurrence in loci\tBaiticity\tDefinition\n")

    print "Writing to file:", out_file
    with open(out_file, 'w') as outf:

        outf.write("Profile\tOccurrence in loci(count)\tBaiticity(count)\tOccurrence in loci(weight)\tBaiticity(weight)\tDefinition\n")

        for profile in profile2orgs2obj:

            in_locus_count = 0
            everywhere_count = 0
            in_locus_weight = 0
            everywhere_weight = 0

            for org in profile2orgs2obj[profile]:
                _org = profile2orgs2obj[profile][org]

                in_locus_count   += _org.in_locus
                everywhere_count += (_org.in_locus + _org.outside)

                in_locus_weight   += _org.in_locus * gnm2weight[org]
                everywhere_weight += (_org.in_locus + _org.outside) * gnm2weight[org]

            _baiticity_count  = 1.0 * in_locus_count / everywhere_count
            _baiticity_weight = in_locus_weight / everywhere_weight

            if everywhere_weight < filter_threshold:

                rare_profiles_file.write("%s\t%f\t%f\t%f\t%s\n"%(profile, everywhere_count, in_locus_count, _baiticity_count, profile2def[profile]))
                continue

            in_loci_count.append(in_locus_count)
            baiticity_count.append(_baiticity_count)
            in_loci_weight.append(in_locus_weight)
            baiticity_weight.append(_baiticity_weight)

            profiles.append(profile)
            outf.write("%s\t%f\t%f\t%f\t%f\t%s\n"%(profile,
                                                    in_locus_count,
                                                    _baiticity_count,
                                                    in_locus_weight,
                                                    _baiticity_weight,
                                                    profile2def[profile]))

    in_loci_weight = np.asarray(in_loci_weight)
    in_loci_weight = np.log10(in_loci_weight)
    baiticity_weight = np.asarray(baiticity_weight)

    plt.ioff()
    fig, ax = plt.subplots()
    ax.scatter(in_loci_weight, baiticity_weight, s=1)

    plt.xlabel("Effective orcurrence in loci (log10)")
    plt.ylabel("Baiticity")

    image_file = os.path.join(save_path, 'baiticity.png')
    plt.savefig(image_file)

    # for i, profile in enumerate(profiles_all):
    #     ax.annotate(profile, (in_loci_all[i], crispricity_all[i]))
    # fig.savefig('second.png')
    # plt.savefig('second.png')

    rare_profiles_file.close()


if __name__=='__main__':

    # prok1402_path_file = '/mnt/storage/data/CDD/all_Prok1402.ccp.csv'
    prok1402_path_file = os.path.join(os.path.expanduser('~'),'data/CDD/all_Prok1402.ccp.csv')

    pan_data_path = '/panfs/pan1.be-md.ncbi.nlm.nih.gov/patternquest/Projects/NewSystems/data/Bacteria/'

    fpath = '/panfs/pan1/patternquest/Projects/NewSystems/data/Bacteria/CDD/cdd_integrase.txt'
    integrases = [l.split()[1] for l in open(fpath).readlines()]
    fpath = '/panfs/pan1/patternquest/Projects/NewSystems/data/Bacteria/CDD/cdd_recombinase.txt'
    recombinases = [l.split()[1] for l in open(fpath).readlines()]
    fpath = '/panfs/pan1/patternquest/Projects/NewSystems/data/Bacteria/CDD/cdd_resolvase.txt'
    resolvases = [l.split()[1] for l in open(fpath).readlines()]
    fpath = '/panfs/pan1/patternquest/Projects/NewSystems/data/Bacteria/CDD/cdd_transpos.txt'
    transposases = [l.split()[1] for l in open(fpath).readlines()]

    pty_files_path = pan_data_path + 'genes_and_flanks/win_10/raw_nbr_files/'

    print "Loading loci"

    pickle_file = os.path.join(pan_data_path, 'pickle/10000/loci.p.bz2')

    bacteria_loci = t.load_compressed_pickle(pickle_file)
    # all_bait_profiles = integrases + recombinases + resolvases + transposases

    # bacteria_loci = []
    #
    # cnt = 0
    # for f in os.listdir(pty_files_path):
    #     bacteria_loci.append(dt.get_pty_file(pty_files_path + f))
    #     if cnt % 1000 == 0:
    #         print cnt
    #     cnt += 1
    #
    # print "Dumping loci"
    # t.dump_compressed_pickle(pickle_file, bacteria_loci)
    # print "Finished dumping"

    loci_organisms = set([locus[0].organism for locus in bacteria_loci])
    loci_gis       = set([gene.gid for locus in bacteria_loci for gene in locus])

    filter_threshold = 3

    arg = int(sys.argv[1])

    if arg == 1:

        save_path = os.path.join(gv.project_data_path, 'baiticity/bacteria/transposase/')
        bait_profiles = transposases

    elif arg == 2:

        save_path = os.path.join(gv.project_data_path, 'baiticity/bacteria/resolvase/')
        bait_profiles = resolvases

    elif arg == 3:

        save_path = os.path.join(gv.project_data_path, 'baiticity/bacteria/recombinase/')
        bait_profiles = recombinases

    elif arg == 4:

        save_path = os.path.join(gv.project_data_path, 'baiticity/bacteria/integrase/')
        bait_profiles = integrases

    elif arg == 5:

        save_path = os.path.join(gv.project_data_path, 'baiticity/bacteria/all/')
        bait_profiles = integrases + recombinases + resolvases + transposases

    print arg, save_path

    calculate_profile_based_baiticity(bacteria_loci,
                                      loci_gis,
                                      loci_organisms,
                                      prok1402_path_file,
                                      bait_profiles,
                                      filter_threshold,
                                      save_path)

    # all_orgs = set()
    # for profile, org2gis in global_profile2orgs2gis.items():
    #     if 'Yersinia_pseudotuberculosis_PB1__uid59153' in org2gis:
    #         print "Yersinia_pseudotuberculosis_PB1__uid59153"
    #     for _org in org2gis.keys():
    #         all_orgs.update([_org])