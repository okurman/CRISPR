


class Locus(object):
    def __init__(self, genome, contig, locus_from, locus_to, genes, cas8_cluster=None,
                  strand=None):
        self.genome = genome
        self.contig = contig
        self.locus_from = locus_from
        self.locus_to = locus_to
        self.genes = genes
        self.orders = set([gene.order for gene in genes])
        self.length = abs(locus_from - locus_to)
        self.label = "%s %s %d..%d" % (genome, contig, locus_from, locus_to)

        self.cas8_cluster = cas8_cluster
        self.strand=strand


    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}".format(self.genome, self.contig, self.l_from, self.l_to, self.length)


class Gene(object):
    def __init__(self, gi, p_from, p_to, strand, name, profile, order=None):
        self.gi = gi
        self.p_from = p_from
        self.p_to = p_to
        self.strand = strand
        self.name = name
        self.profile = profile
        self.length = abs(self.p_from - self.p_to)
        self.order = order

    def shift_back(self, dist):
        self.p_from = self.p_from - dist
        self.p_to = self.p_to - dist

    def __str__(self):
        terms = [self.gi, self.p_from, self.p_to, self.length, self.strand, self.profile, self.name]
        return "\t".join(str(t) for t in terms)

    def __cmp__(self, other):
        if self.p_from>other.p_from:
            return 1
        elif self.p_from<other.p_from:
            return -1
        else:
            return 0

    # def overlaps(self, other):
    #     dist1 = other.p_from - self.p_to
    #     dist2 = self.p_from - other.p_to
    #     return True if dist1*dist2>=0 else False


def load_locus(lines):

    parts = lines[0].split("\t")

    genome = parts[3]
    contig = parts[4]

    if len(parts)>8:
        cas8_cl_id = parts[9]
        locus_strand = parts[10].strip()

    genes = []

    for l in lines[1:]:
        parts = l.rstrip().split("\t")

        gi = parts[0]
        coords = parts[1]
        p_from, p_to = coords.split("..")
        p_from, p_to = int(p_from), int(p_to)
        strand = parts[2]
        profile = parts[6]
        name = parts[8]

        order = parts[18] if len(parts)>18 else None

        _gene = Gene(gi=gi,
                     p_from=p_from,
                     p_to=p_to,
                     strand=strand,
                     profile=profile,
                     name=name,
                     order=order)
        genes.append(_gene)

    locus_from = genes[0].p_from
    locus_to = genes[-1].p_to
    
    locus = Locus(genome=genome,
                  contig=contig,
                  locus_from=locus_from,
                  locus_to=locus_to,
                  genes=genes,
                  cas8_cluster=cas8_cl_id,
                  strand=locus_strand)

    return locus


def load_islands(islands_file):

    all_lines = open(islands_file).readlines()
    loci = []

    start_points = [i for i in range(len(all_lines)) if all_lines[i].startswith("===")]

    for i in range(len(start_points)):

        if i < len(start_points)-1:
            block = all_lines[start_points[i]:start_points[i+1]]
        else:
            block = all_lines[start_points[i]:]

        loci.append(load_locus(block))

    return loci