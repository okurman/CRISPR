__author__ = 'Sanjarbek Hudaiberdiev'

import sys
import os
if sys.platform=='darwin':
    sys.path.append('/Users/hudaiber/Projects/SystemFiles/')
elif sys.platform=='linux2':
    sys.path.append('/home/hudaiber/Projects/SystemFiles/')
import global_variables as gv

from lib.utils import tools as t
from lib.utils.classes import Neighborhood


def get_flank_distributions(kplets_2d_list, neighborhood_path, target_profiles):

    org2weights = t.map_genome2weight()
    flanking_genes_count = []

    cog2gids = []

    gid2weight = dict()

    for kplets_list in kplets_2d_list:
        cur_flanking_genes_count = dict()

        cur_cog2gids = dict()

        for kplet in kplets_list:
            neighborhoods = [Neighborhood(os.path.join(neighborhood_path, f)) for f in kplet.files]

            for neighborhood in neighborhoods:
                for gene in neighborhood.genes:

                    gid2weight[int(gene.gid)] = org2weights[gene.organism]

                    for cogid in gene.cogid.split():
                        # if cogid in target_profiles:
                        #     continue
                        t.update_dictionary(cur_flanking_genes_count,cogid,org2weights[gene.organism])
                        t.update_dictionary_set(cur_cog2gids, cogid, set([int(gene.gid)]))

        flanking_genes_count.append(cur_flanking_genes_count)
        cog2gids.append(cur_cog2gids)

    return flanking_genes_count, cog2gids, gid2weight