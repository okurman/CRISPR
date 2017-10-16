import subprocess as sp

class Hit(object):
    def __init__(self, hit_line):

        parts = hit_line.split("\t")
        self.hit_line = hit_line
        self.perc_id = float(parts[2])
        self.length = int(parts[3])
        self.q_start = int(parts[6])
        self.q_end = int(parts[7])
        self.q_len = abs(self.q_end - self.q_start)
        self.s_start = int(parts[8])
        self.s_end = int(parts[9])
        self.s_len = abs(self.s_end - self.s_start)

        self.slope = self.s_len / float(self.q_len)

    def overlaps(self, other):

        # of the hit is on the same direction
        if other.s_start < other.s_end:

            dist1 = other.q_start - self.q_end
            dist2 = self.q_start - other.q_end
            q_overlaps = True if dist1 * dist2 >= 0 else False

            dist1 = other.s_start - self.s_end
            dist2 = self.s_start - other.s_end
            s_overlaps = True if dist1 * dist2 >= 0 else False

            return True if q_overlaps and s_overlaps else False

        # inverse hits.
        else:

            dist1 = other.q_start - self.q_end
            dist2 = self.q_start - other.q_end
            q_overlaps = True if dist1 * dist2 >= 0 else False

            other_s_start = other.s_end
            other_s_end = other.s_start

            self_s_start = self.s_end
            self_s_end = self.s_start

            dist1 = other_s_start - self_s_end
            dist2 = self_s_start - other_s_end

            s_overlaps = True if dist1 * dist2 >= 0 else False

            return True if q_overlaps and s_overlaps else False

    def add_padding(self, dist):

        self.q_start += dist
        self.q_end += dist
        self.s_start += dist
        self.s_end += dist


def read_blast_hits(blast_file, THR_ID=70, SHORT_HIT_LENGTH=100):

    all_hits = [Hit(l) for l in open(blast_file) if not l.startswith("#")]
    distinct_hits = []

    for hit in all_hits:

        if hit.perc_id < THR_ID:
            continue

        if hit.length < SHORT_HIT_LENGTH:
            continue

        overlap_found = False

        for i, distinct_hit in enumerate(distinct_hits):
            if distinct_hit.overlaps(hit):
                distinct_hits[i] = distinct_hit if distinct_hit.length >= hit.length else hit
                overlap_found = True

        if not overlap_found:
            distinct_hits.append(hit)

    return distinct_hits


def megablast_query_subject(query_file, subject_file, out_file):

    cmd_fmt = "blastn -task megablast -query {query_file} -subject {subject_file} -outfmt 7 -out {out_file}"
    cmd = cmd_fmt.format(query_file=query_file,
                         subject_file=subject_file,
                         out_file=out_file)
    sp.call(cmd, shell=True)

