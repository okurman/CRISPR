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

            if len(terms) < 7:
                continue

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
                                      arcog_path_file,
                                      bait_profiles,
                                      filter_threshold,
                                      save_path):

    print "Loding global maps"
    global_profile2orgs2gis = load_maps_simple(arcog_path_file, loci_gis)
    print "Loading weights"
    gnm2weight = t.map_genome2weight()
    print "Loading CDD definitions"
    profile2def = t.map_cdd_profile2def()

    profile2def.update(t.map_profile2def())

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

    print len(profile2orgs2obj['arCOG08578'])
    # print profile2orgs2obj['arCOG08578'].keys()

    for org, obj in profile2orgs2obj['arCOG08578'].items():
        if obj.in_locus + obj.outside > 0:
            print org, obj.in_locus, obj.outside

    sys.exit()
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

            if profile=='arCOG14077':
                continue

            in_locus_count = 0
            everywhere_count = 0
            in_locus_weight = 0
            everywhere_weight = 0

            for org in profile2orgs2obj[profile]:

                if org in ['Nitrosoarchaeum_koreensis_MY1_MY1','Nitrosoarchaeum_limnia_SFB1']:
                    continue

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

    arcog_path_file = os.path.join(os.path.expanduser('~'),'data/Archea/arCOG/ar14.arCOG.csv')

    pan_data_path = '/panfs/pan1/patternquest/Projects/NewSystems/data/Archea/'

    fpath = '/panfs/pan1/patternquest/Projects/NewSystems/data/Archea/arCOG/arcogs_integrase.txt'
    integrases = [l.strip() for l in open(fpath).readlines()]
    fpath = '/panfs/pan1/patternquest/Projects/NewSystems/data/Archea/arCOG/arcogs_recombinase.txt'
    recombinases = [l.strip() for l in open(fpath).readlines()]
    fpath = '/panfs/pan1/patternquest/Projects/NewSystems/data/Archea/arCOG/arcogs_transposase.txt'
    transposases = [l.strip() for l in open(fpath).readlines()]

    pty_files_path = pan_data_path + 'genes_and_flanks/win_10/pty/'

    print "Loading loci"

    pickle_file = os.path.join(pan_data_path, 'pickle/10000/loci.p.bz2')

    loci = t.load_compressed_pickle(pickle_file)

    # all_bait_profiles = integrases + recombinases + resolvases + transposases

    # loci = []
    #
    # cnt = 0
    # for f in os.listdir(pty_files_path):
    #     loci.append(dt.get_pty_file(pty_files_path + f))
    #     if cnt % 1000 == 0:
    #         print cnt
    #     cnt += 1
    #
    # print "Dumping loci"
    # t.dump_compressed_pickle(pickle_file, loci)
    # print "Finished dumping"


    loci_organisms = set([locus[0].organism for locus in loci])
    loci_gis       = set([gene.gid for locus in loci for gene in locus])

    filter_threshold = 3

    save_path = os.path.join(gv.project_data_path, 'baiticity/archaea/all/')
    bait_profiles = integrases + recombinases + transposases

    calculate_profile_based_baiticity(loci,
                                      loci_gis,
                                      loci_organisms,
                                      arcog_path_file,
                                      bait_profiles,
                                      filter_threshold,
                                      save_path)

    # all_orgs = set()
    # for profile, org2gis in global_profile2orgs2gis.items():
    #     if 'Yersinia_pseudotuberculosis_PB1__uid59153' in org2gis:
    #         print "Yersinia_pseudotuberculosis_PB1__uid59153"
    #     for _org in org2gis.keys():
    #         all_orgs.update([_org])