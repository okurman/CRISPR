# __author__ = 'Sanjarbek Hudaiberdiev'

import sys
import os
# sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/lib/BioPy/'))
# sys.path.append(os.path.join(os.path.expanduser('~'),'Projects/SystemFiles/'))
import configparser

###############################################################
config_file = os.path.join(os.path.expanduser('~'), 'paths.cfg')
cfg = configparser.ConfigParser()
cfg.read(config_file)

code_path = cfg.get('NewSystems', 'code_path')
data_path = cfg.get('NewSystems', 'data_path')
prok1402_path = cfg.get('prok1402', 'db_path')
sys.path.append(code_path)

code_path = cfg.get('BioPy', 'code_path')
sys.path.append(code_path)
###############################################################

import os
from operator import itemgetter
import cPickle
import bz2
import csv
# from . import CrisprLocus
# from BioClasses import CrisprGene
from BioClasses import Gene
from collections import defaultdict
import hhpred
import networkx as nx
import shutil

cdd_def_file = '/net/frosty/vol/export1/cdd/current.version/entrez/cdfind_pub.dat'


def target_profiles():
    profiles_file = os.path.join(gv.project_data_path, 'Archea/arCOG/selected_arcogs.txt')

    if not os.path.exists(profiles_file):
        profiles_file = '/panfs/pan1.be-md.ncbi.nlm.nih.gov/patternquest/Projects/NewSystems/data/Archea/arCOG/selected_arcogs.txt'

    if not os.path.exists(profiles_file):
        raise IOError("file selected_arcogs wasn't found")

    return [l.strip() for l in open(profiles_file).readlines()]


def bacteria_target_profiles():
    profiles_file = os.path.join(gv.project_data_path, 'Bacteria/CDD/profile_ids.txt')

    if not os.path.exists(profiles_file):
        profiles_file = '/panfs/pan1.be-md.ncbi.nlm.nih.gov/patternquest/Projects/NewSystems/data/Bacteria/CDD/profile_ids.txt'

    return [l.strip() for l in open(profiles_file).readlines()]


def map_src2org():
    return {l.split()[1]: l.split()[0] for l in
            open(os.path.join(gv.data_path, 'info', 'map_gnm_src.txt')).readlines()}


def map_genome2weight():
    return {l.split()[0]: float(l.split()[1]) for l in
            open(os.path.join(gv.data_path, 'CDD', 'Prok1402_ad.weight.tab')).readlines()}


def map_profile2def():
    def_map = {l.split('\t')[0]: l.split('\t')[3] for l in
               open(os.path.join(gv.data_path, 'Archea/arCOG/ar14.arCOGdef.tab')).readlines()}
    def_map["-"] = " "
    return def_map


def map_cdd_profile2def():
    # def_map = {l.split('\t')[1]: l.split('\t')[3] for l in
    #            open(os.path.join(gv.data_path, 'CDD/cdfind_pub_ad.dat')).readlines()}
    def_map = {l.split('\t')[1]: l.split('\t')[3] for l in open(cdd_def_file).readlines()}
    def_map["-"] = " "

    return def_map


def map_defense_profile2gene_name(defense_map_file=os.path.join(data_path, 'cas4/profiles/defenseProfiles.tab')):
    def_map = {}

    for l in open(defense_map_file):
        terms = l.split('\t')
        profile = terms[0]
        gene_names = terms[3].split(',')
        if len(gene_names) > 1:
            def_map[profile] = gene_names[1]
        else:
            def_map[profile] = gene_names[0]

    return def_map


def map_cdd_defense2gene_name():
    cdd_map = map_cdd_profile2gene_name()

    cdd_map.update(map_defense_profile2gene_name())

    return cdd_map


def map_cdd_profile2gene_name():
    # def_map = {l.split('\t')[1]: l.split('\t')[2] for l in
    #            open(os.path.join(gv.data_path, 'CDD/cdfind_pub_ad.dat')).readlines()}
    def_map = defaultdict(str)
    for l in open(cdd_def_file).readlines():
        def_map[l.split('\t')[1]] = l.split('\t')[2].strip()

    def_map["-"] = " "
    def_map[""] = " "

    def_map['cd09637'] = 'Cas4'
    def_map['cd09659'] = 'Cas4'
    def_map['cl00641'] = 'Cas4'

    return def_map


def map_cdd_ar14_profile2def():
    def_map = map_cdd_profile2def()
    def_map.update(map_profile2def())

    return def_map


def update_dictionary(map, key, value):
    if key in map:
        map[key] += value
    else:
        map[key] = value


def update_dictionary_set(map, key, value):
    if key in map:
        if isinstance(value, set):
            map[key].update(value)
        else:
            map[key].update({value})
    else:
        if isinstance(value, set):
            map[key] = value
        else:
            map[key] = {value}


def update_dictionary_list(map, key, value):
    if key not in map:
        map[key] = []

    if isinstance(value, list):
        map[key] += value
    else:
        map[key].append(value)


def merge(d1, d2, merge_fn):
    """
    Merges two dictionaries, non-destructively, combining
    values on duplicate keys as defined by the optional merge
    function.  The default behavior replaces the values in d1
    with corresponding values in d2.

    Examples:

    >>> d1
    {'a': 1, 'c': 3, 'b': 2}
    >>> merge(d1, d1, lambda x,y: y)
    {'a': 1, 'c': 3, 'b': 2}
    >>> merge(d1, d1)
    {'a': 2, 'c': 6, 'b': 4}
    """
    result = dict(d1)
    for k, v in d2.items():
        if k in result:
            result[k] = merge_fn(result[k], v)
        else:
            result[k] = v
    return result


def merge_set_val(d1, d2):
    result = dict(d1)
    for k, v in d2.items():
        if k in result:
            result[k].update(v)
        else:
            result[k] = v
    return result


def merge_dict_list(dict_list, merge_fn=lambda x, y: x + y):
    """
    Merge list of dictionaries into one dictionary.
    It uses the merge(d1, d2, merge_fn=lambda x,y: x+y)
    function from above. For the purpose of this library,
    it uses the default function for merging, which is addition.
    """
    result = dict()
    for _dict in dict_list:
        result = merge(result, _dict, merge_fn)
    return result


def merge_dict_set_list(dict_list, gid2weights):
    result = dict()
    for _dict in dict_list:
        result = merge_set_val(result, _dict)

    for k, v in result.items():
        result[k] = sum(gid2weights[gid] for gid in v)

    return result


# def merge_cog2gids_list(cog2gids_list):
#
#     result = dict()
#
#     for _dict in cog2gids_list:
#
#


def map_arcog2class():
    _profile2class_code = {l.split('\t')[0]: l.split('\t')[1] for l in
                           open(os.path.join(gv.data_path, 'Archea/arCOG/ar14.arCOGdef.tab')).readlines()}
    _class2def = {l.split('\t')[0]: l.split('\t')[2].strip() for l in
                  open(os.path.join(gv.data_path, 'Archea/arCOG/ar14/funclass.tab'))}

    _profile2class_def = dict()
    for _profile, _class_code in _profile2class_code.items():
        for _sub_class_code in _class_code:
            _class_def = _class2def[_sub_class_code]
            # update_dictionary_list_value(_profile2class_def, _profile, _class_def)
            update_dictionary(_profile2class_def, _profile, [_class_def])

    return _profile2class_def


def map_cdd2class():
    base_path = os.path.join(gv.project_code_path, 'scripts/cdd_clusters')
    cdd_class_file = os.path.join(base_path, 'cdd_to_class.tab')
    cdd_class_def_file = os.path.join(base_path, 'fun2003-2014.tab')

    _profile2class_code = {l.split('\t')[0]: l.split('\t')[1].strip() for l in open(cdd_class_file).readlines()}
    _class2def = {l.split('\t')[0]: l.split('\t')[1].strip() for l in open(cdd_class_def_file).readlines()}

    _profile2class_def = dict()
    for _profile, _class_code in _profile2class_code.items():
        for _sub_class_code in _class_code:
            _class_def = _class2def[_sub_class_code]
            update_dictionary(_profile2class_def, _profile, [_class_def])

    return _profile2class_def


def map_genome2weight(dataset="prok1402"):
    _genome2weight = {}

    if dataset == "prok1402":
        _genome2weight = {l.split()[0]: float(l.split()[1]) for l in
                          open(os.path.join(prok1402_path, 'Prok1402.weight.tab')).readlines()}

        _genome2weight["Arthrobacter_Rue61a_uid174511"] = 3.0814e-01
        _genome2weight["Legionella_pneumophila_LPE509_uid193710"] = 3.9054e-02
        _genome2weight["Magnetospirillum_gryphiswaldense_MSR_1_uid232249"] = 1.9391e+00
        _genome2weight["Oscillibacter_valericigenes_Sjm18_20_uid73895"] = 4.0300e+00
        _genome2weight["Klebsiella_pneumoniae_rhinoscleromatis_SB3432_uid203334"] = 1.5744e-01

        _genome2weight["Nitrosoarchaeum_koreensis_MY1_MY1"] = 0  # missing
        _genome2weight["Nitrosoarchaeum_limnia_SFB1"] = 0  # missing
        _genome2weight["Streptomyces_cattleya_NRRL_8057___DSM_46488_uid77117"] = 3.1135e-01
        _genome2weight["nanoarchaeote_Nst1"] = 0  # missing

    return _genome2weight


def map_gi2profiles():
    ccp_path = "/panfs/pan1/compgen/cog/"
    gi2profiles = defaultdict(set)

    for f in os.listdir(ccp_path):

        if not f.endswith("ccp.csv"):
            continue

        for l in open(os.path.join(ccp_path, f)):
            parts = l.split(",")
            _gi = parts[0]
            _profile = parts[6]
            update_dictionary_set(gi2profiles, _gi, _profile)

    return gi2profiles


def map_archaea_genome2weight():
    return {l.split()[0]: float(l.split()[2]) for l in
            open(os.path.join(gv.data_path, 'Archea/arCOG/ar14.genclade.tab')).readlines()}


def map_gid2src(map_file):
    out_map = {}
    with open(map_file) as f:
        for l in f:
            terms = l.split('\t')
            src = terms[1]
            gids = terms[2]
            for gid in gids.strip().split():
                out_map[gid] = src
    return out_map


def map_gid2arcog():
    return {l.split(',')[0]: l.split(',')[6] for l in
            open(os.path.join(gv.data_path, 'Archea/arCOG/ar14.arCOG.csv')).readlines() if 'arCOG' in l}


def map_gid2arcog_cdd():
    _gid2arcog = map_gid2arcog()
    _gid2cdd = map_gid2cdd()

    return dict(_gid2cdd.items() + _gid2arcog.items())


def get_weighted_profiles_from_neighborhoods(neighborhoods_path, exclude_target=True):
    if exclude_target:
        _target_profiles = target_profiles()
    _src2org = map_src2org()
    _gnm2weight = map_genome2weight()
    neighborhoods = [cl.Neighborhood(os.path.join(neighborhoods_path, f)) for f in os.listdir(neighborhoods_path)]

    profile_stats = {}
    for nbr in neighborhoods:
        src_name = nbr.genes[0].src
        org_name = _src2org[src_name]
        org_weight = _gnm2weight[org_name] if org_name in _gnm2weight else 1
        for g in nbr.genes:
            if g.cogid == "":
                continue
            for tmpCog in g.cogid.split():
                if exclude_target and tmpCog in _target_profiles:
                    continue
                if tmpCog in profile_stats:
                    if g.gid not in profile_stats[tmpCog].gids:
                        profile_stats[tmpCog].weight += org_weight
                        profile_stats[tmpCog].count += 1
                        profile_stats[tmpCog].gids += [g.gid]
                else:
                    profile_stats[tmpCog] = cl.ProfileCount(1, org_weight, g.gid)

    profile_weights = [(k, v.weight) for k, v in profile_stats.items()]
    profile_weights = sorted(profile_weights, key=itemgetter(1), reverse=True)

    return profile_weights


def load_neighborhoods(path, target_files=None):
    import classes as cl

    if not target_files:
        files = [os.path.join(path, f) for f in os.listdir(path)]
    else:
        files = [os.path.join(path, f) for f in target_files]

    return [cl.Neighborhood(f) for f in files]


def load_compressed_pickle(fname, compression='bz2'):

    if (isinstance(fname, str) or isinstance(fname, unicode)) and compression == 'bz2':
        f = bz2.BZ2File(fname, 'rb')
        retval = cPickle.load(f)
    else:
        raise NotImplementedError

    return retval


def dump_compressed_pickle(fname, data, compression='bz2'):

    if (isinstance(fname, str) or isinstance(fname, unicode)) and compression == 'bz2':
        f = bz2.BZ2File(fname, 'wb')
        cPickle.dump(data, f, protocol=cPickle.HIGHEST_PROTOCOL)
    else:
        raise NotImplementedError


def map_file2organism():
    retval = dict()
    fname = os.path.join(gv.project_data_path, 'Archea/genes_and_flanks/win_10/filename_source_organism.tab')
    for l in open(fname):
        terms = l.strip().split()
        retval[terms[0]] = terms[2]

    fname = os.path.join(gv.project_data_path, 'Bacteria/genes_and_flanks/win_10/filename_source_organism.tab')
    for l in open(fname):
        terms = l.strip().split()
        retval[terms[0]] = terms[2]

    return retval


def map_wgs_profile_count(dataset):
    # Before filtration, After filtration
    bf, af = {}, {}
    if dataset == 1:
        fname = os.path.join(gv.project_data_path, 'CRISPR/datasets/cas1_1/wgs_profile_count.tab')
        bf = {l.strip().split()[1]: l.strip().split()[0] for l in open(fname).readlines()}
        fname = os.path.join(gv.project_data_path, 'CRISPR/datasets/cas1_1/wgs_profile_count_af.tab')
        bf = {l.strip().split()[1]: l.strip().split()[0] for l in open(fname).readlines()}
    elif dataset == 2:
        fname = os.path.join(gv.project_data_path, 'CRISPR/datasets/cas1_2/wgs_profile_count.tab')
        bf = {l.strip().split()[1]: l.strip().split()[0] for l in open(fname).readlines()}
        fname = os.path.join(gv.project_data_path, 'CRISPR/datasets/cas1_2/wgs_profile_count_af.tab')
        af = {l.strip().split()[1]: l.strip().split()[0] for l in open(fname).readlines()}
    elif dataset == 3:
        fname = os.path.join(gv.project_data_path, 'CRISPR/datasets/crispr/wgs_profile_count.tab')
        bf = {l.strip().split()[1]: l.strip().split()[0] for l in open(fname).readlines()}
        fname = os.path.join(gv.project_data_path, 'CRISPR/datasets/crispr/wgs_profile_count_af.tab')
        af = {l.strip().split()[1]: l.strip().split()[0] for l in open(fname).readlines()}

    return bf, af


def map_global_cdd_profile_count():
    map_file = os.path.join(os.path.expanduser('~'), 'data/CDD/profile2weight.tab')
    profile2weight = {l.split()[0]: float(l.split()[1]) for l in open(map_file).readlines() if not l.startswith("#")}

    return profile2weight


def write_kplets_to_csv(kplets, out_name, compression=False):
    if compression == 'True':
        fout = bz2.BZ2File(out_name, "w")
    else:
        fout = open(out_name, "w")

    csv_writer = csv.writer(fout)
    csv_writer.writerow(('Id', 'K', 'Count', 'Codes', 'Files'))
    for kplet in kplets:
        row = []
        row += [kplet.id, kplet.k, kplet.count, ]
        row += [' '.join(list(kplet.codes)), ' '.join(kplet.files)]
        csv_writer.writerow(row)


def write_genes_to_pty(genes, file_name, header_lines=None):

    src = genes[0].src
    org = genes[0].organism

    column_names = "#GI\tCoordinates\tStrand\tOrganism\tChromosome\tCDD\tGene\n"

    with open(file_name, 'wb') as outf:

        if header_lines:
            [outf.write("#%s\n" % l) for l in header_lines]

        outf.write(column_names)

        for gene in genes:

            coords = "%d..%d" % (gene.pFrom, gene.pTo)

            if not gene.profiles:
                _profiles = ""

            elif type(gene.profiles) == str:
                _profiles = gene.profiles
            else:
                _profiles = " ".join(gene.profiles)

            if gene.gene_names:
                _gene_names = gene.gene_names.strip()
            else:
                _gene_names = ""

            pty_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene.gid, coords, gene.strand, org, src, _profiles, _gene_names)

            outf.write(pty_line)


def parse_crispr_loci(source_file, separator="===", filter_set=[]):
    def build_locus(header_line, block_lines):

        parts = header_line.rstrip().split('\t')

        organism = parts[3]
        source = parts[4]
        type = parts[5]
        status = parts[7]
        complete = True if parts[6].lower() == "complete" else False
        id = parts[8].split()[1]

        genes = []

        for _line in block_lines[1:]:
            parts = _line.strip().split('\t')

            gi = parts[0]
            pFrom, pTo = parts[1].split('..')
            strand = parts[2]
            profile = parts[6]
            gene_name = parts[8]
            is_seed = True if parts[10] == "Seed" else False
            cluster_id = parts[17]

            _gene = CrisprGene(source, pFrom, pTo, gid=gi, strand=strand, cogid=profile,
                               gene_name=gene_name, cluster_id=cluster_id, is_seed=is_seed)

            genes.append(_gene)

        _locus = CrisprLocus(organism=organism,
                             source=source,
                             genes=genes,
                             status=status,
                             type=type,
                             complete=complete,
                             id=id)

        return _locus

    loci = []

    with open(source_file) as inf:

        cnt = 0

        l = inf.readline()
        block_lines = []
        while l:

            if l.startswith(separator):

                if not block_lines:
                    header_line = l.strip()
                    block_lines.append(l)

                else:

                    parts = header_line.rstrip().split('\t')

                    organism = parts[3]

                    if filter_set and organism not in filter_set:
                        block_lines = [l]
                        header_line = l
                        l = inf.readline()
                        continue

                    loci.append(build_locus(header_line, block_lines))

                    header_line = l
                    block_lines = [l]

                    cnt += 1
                    # if cnt == 10:
                    #     break

            else:

                block_lines.append(l)

            l = inf.readline()

        loci.append(build_locus(header_line, block_lines))

    return loci


def parse_plain_pty_file(_f):

    """
    Parses plain pty files generated trom prok_extract. Unless the name includes "plain_pty", the genes should be searched for annotation information.

    :param file with pty lines:
    :return list of BioClass.Gene objects:
    """
    _genes = []

    for l in open(_f):

        if l.startswith("#"):
            continue

        parts = l.rstrip().split("\t")

        gi = parts[6]
        pFrom, pTo = parts[1].split('..')
        strand = parts[2]
        organism = parts[3]
        source = parts[4]

        _genes.append(Gene(source, pFrom, pTo, gi=gi, organism=organism, strand=strand))

    return _genes


def parse_pty_file(_f):
    _genes = []

    for l in open(_f):

        if l.startswith("#"):
            continue

        parts = l.rstrip().split("\t")

        gi = parts[0]
        pFrom, pTo = parts[1].split('..')
        strand = parts[2]
        organism = parts[3]
        source = parts[4]

        profiles = None
        gene_names = None

        if len(parts) > 5:
            profiles = parts[5]
        if len(parts) > 6:
            gene_names = parts[6]

        _genes.append(Gene(source, pFrom, pTo, gi=gi, organism=organism, strand=strand, profiles=profiles,
                           gene_names=gene_names))

    return _genes


def write_pty_file(genes, file_name):
    OUTFMT = "{gi}\t{coords}\t{strand}\t{genome}\t{source}\t{profiles}\n"

    with open(file_name, "w") as outf:
        for gene in genes:
            out_line = OUTFMT.format(gi=gene.gid,
                                     coords="%d..%d" % (gene.pFrom, gene.pTo),
                                     strand=gene.strand,
                                     genome=gene.organism,
                                     source=gene.src,
                                     profiles=" ".join(gene.profiles))
            outf.write(out_line)


def profile_matrix():

    THR_PROB = 99.0
    THR_EVAL = 1E-10
    THR_COV = 0.8

    source_dir = "/panfs/pan1/patternquest/data/CDD/selfmat/"
    save_file = "/panfs/pan1/patternquest/Projects/NewSystems/data/profiles/profile2profile_mat.txt.bz2"

    pr2pr = {}
    print "Scanning files"

    cnt = 1
    for f in os.listdir(source_dir):

        if cnt % 1000 == 0:
            print cnt

        cnt += 1
        
        _p = os.path.splitext(f)[0]
        hits = hhpred.hhsearch_parse(os.path.join(source_dir, f), THR_PROB, THR_EVAL, THR_COV, THR_COV)

        if f == "COG1074.hhr" or f == "COG0210.hhr":
            print f
            for hit in hits:
                print hit.line

        if not hits:
            continue

        pr2pr[_p] = {}
        for hit in hits:
            pr2pr[_p][hit.profile] = (hit.prob, hit.evalue, hit.score)

    if os.path.exists(save_file):
        backup_file = save_file+".backup"
        print "The file to be overritten is backed up here:", backup_file
        shutil.copy(save_file, backup_file)

    print "Saving to file", save_file
    dump_compressed_pickle(save_file, pr2pr)

    return pr2pr


def profiles_detect_synonyms():

    """The following file is generated by the function profile_matrix() above. """
    _src_file = "/panfs/pan1/patternquest/Projects/NewSystems/data/profiles/profile2profile_mat.txt.bz2"
    pr2pr = load_compressed_pickle(_src_file)

    G = nx.MultiGraph()

    for profile in pr2pr.keys():
        for synonym in pr2pr[profile]:
            if profile != synonym:
                G.add_edge(profile, synonym)

    save_file = "/panfs/pan1/patternquest/Projects/NewSystems/data/profiles/synonyms.txt"
    if os.path.exists(save_file):
        backup_file = save_file+".backup"
        print "The file to be overritten is backed up here:", backup_file
        shutil.copy(save_file, backup_file)

    print "Saving to file:", save_file
    with open(save_file, "w") as outf:

        for component in sorted(nx.connected_components(G), key=lambda x: len(x), reverse=True):
            outf.write("%d\t%s\n" % (len(component), " ".join(component)))

    return pr2pr


if __name__ == '__main__':
    # file2organism = map_file2organism()

    # pass
    # source_path = '/Users/hudaiber/Projects/NewSystems/data/Archea/genes_and_flanks/win_10/pty/'
    # target_file = 'src_to_gids.txt'
    #
    # build_src2gids_map(source_path, target_file)

    # pr2pr = profile_matrix()

    profiles_detect_synonyms()