import gzip

from cart2enst_.helper import Transcript, Exon


def read_transcript_db(fn):

    ret = {}

    for line in gzip.open(fn):
        line = line.strip()
        if line == '' or line[0] == '#':
            continue

        t = _parse_transcript_line(line)

        if t.gene_symbol not in ret:
            ret[t.gene_symbol] = []
        ret[t.gene_symbol].append(t)

    return ret


def _parse_transcript_line(line):

    cols = line.split('\t')
    exons = []
    id_ = cols[0]
    gene_symbol = cols[1]
    gene_id = cols[2]
    chrom = cols[4]
    strand = int(cols[5])
    transcript_start = int(cols[6])
    transcript_end = int(cols[7])
    coding_start_genomic = int(cols[9]) - 1

    for i in range(1, len(cols) - 11, 2):
        exons.append(
            Exon(int(cols[10 + i]), int(cols[11 + i]))
        )

    coding_end_genomic = int(cols[10]) - 1

    return Transcript(
        id_=id_,
        gene_symbol=gene_symbol,
        gene_id=gene_id,
        chrom=chrom,
        strand=strand,
        transcript_start=transcript_start,
        transcript_end=transcript_end,
        coding_start_genomic=coding_start_genomic,
        coding_end_genomic=coding_end_genomic,
        exons=exons
    )


def read_gene_synonyms(fn):

    ret = {}

    for line in open(fn):
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        cols = line.split()
        gene_names = [cols[0]] + cols[1].split(',')
        gene_names = [x for x in gene_names if x != '.']
        gene_names = list(set(gene_names))

        for k in gene_names:
            ret[k] = gene_names

    return ret


def read_ensembl_transcripts(ensembl_tabix_file, cart):

    ret = []
    for line in ensembl_tabix_file.fetch(cart.chrom, cart.transcript_start, cart.transcript_end):
        line = line.strip()
        if line == '' or line[0] == '#':
            continue

        t = _parse_transcript_line(line)
        ret.append(t)
    return ret


def read_input_file(fn):

    ret = []
    for line in open(fn):
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        ret.append(line)
    return ret
