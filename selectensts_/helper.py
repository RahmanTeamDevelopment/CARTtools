from __future__ import division
import parser
import utr


class Transcript(object):

    def __init__(
            self,
            id_,
            gene_symbol,
            gene_id,
            chrom,
            strand,
            transcript_start,
            transcript_end,
            coding_start_genomic,
            coding_end_genomic,
            exons
    ):

        self.id_ = id_
        self.gene_symbol = gene_symbol
        self.gene_id = gene_id
        self.chrom = chrom
        self.strand = '+' if strand == 1 else '-'
        self.transcript_start = transcript_start
        self.transcript_end = transcript_end
        self.coding_start_genomic = coding_start_genomic
        self.coding_end_genomic = coding_end_genomic
        self.exons = exons

        self.cds_exons = self.cds_regions()
        self.utr5_exons = self.utr5_regions()
        self.utr3_exons = self.utr3_regions()
        self.utr_exons = self.utr5_exons + self.utr3_exons

    def __str__(self):

        return '{}({})_{}:{}-{}[{}]'.format(self.id_, self.gene_symbol, self.chrom, self.transcript_start,
                                            self.transcript_end, self.strand)

    def cds_regions(self):

        ret = []
        for exon in self.exons:
            cds_region = exon.get_cds(self.coding_start_genomic, self.coding_end_genomic)
            if cds_region is not None:
                ret.append(cds_region)
        return ret

    def utr5_regions(self):

        ret = []
        for exon in self.exons:
            utr5_region = exon.get_utr5(self.coding_start_genomic, self.strand)
            if utr5_region is not None:
                ret.append(utr5_region)
        return ret

    def utr3_regions(self):

        ret = []
        for exon in self.exons:
            utr3_region = exon.get_utr3(self.coding_end_genomic, self.strand)
            if utr3_region is not None:
                ret.append(utr3_region)
        return ret


class Exon(object):

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def get_cds(self, coding_start, coding_end):
        """Return CDS interval or None if there is no CDS in the exon"""

        coding_min = min(coding_start, coding_end)
        coding_max = max(coding_start, coding_end)

        if self.end - 1 < coding_min or self.start > coding_max:
            return

        cds_start = max(self.start, coding_min)
        cds_end = min(self.end - 1, coding_max) + 1
        return cds_start, cds_end

    def get_utr5(self, coding_start, strand):
        """Return UTR5 interval or None if there is no UTR5 in the exon"""

        if strand == '+':
            if self.start >= coding_start:
                return
            return self.start, min(coding_start - 1, self.end - 1) + 1
        else:
            if self.end - 1 <= coding_start:
                return
            return max(self.start, coding_start + 1), self.end

    def get_utr3(self, coding_end, strand):
        """Return UTR3 interval or None if there is no UTR3 in the exon"""

        if strand == '+':
            if self.end - 1 <= coding_end:
                return
            return max(self.start, coding_end + 1), self.end
        else:
            if self.start >= coding_end:
                return
            return self.start, min(self.end - 1, coding_end - 1) + 1

    def __str__(self):
        return '{}-{}'.format(self.start, self.end)

    def __contains__(self, position):
        return self.start <= position < self.end


#######################################################################################################################


def perfect_and_cds_matches_synonyms(gene_synonyms, ensembl_data, nm, gene, out):

    perfect_matches = []
    identical_cds = []

    if gene in gene_synonyms:
        synonyms_in = [g for g in gene_synonyms[gene] if g in ensembl_data]
    else:
        synonyms_in = [gene] if gene in ensembl_data else []

    if len(synonyms_in) == 0:
        output(
            out,
            nm,
            gene,
            None,
            'gene_not_found',
            '.',
            '.',
            '.'
        )
        return None, None

    for g in synonyms_in:

        for enst in ensembl_data[g]:

            if are_perfect_match(nm, enst):
                perfect_matches.append(enst)
                continue

            if have_identical_cds(nm, enst):
                identical_cds.append(enst)

    return perfect_matches, identical_cds


def perfect_and_cds_matches_any_gene(ensembl_tabix_file, nm, gene, out):

    perfect_matches = []
    identical_cds = []

    overlapping_ensts = parser.read_ensembl_transcripts(ensembl_tabix_file, nm)

    if len(overlapping_ensts) == 0:
        output(
            out,
            nm,
            gene,
            None,
            'not_found',
            '.',
            '.',
            '.'
        )
        return None, None

    for enst in overlapping_ensts:

        if are_perfect_match(nm, enst):
            perfect_matches.append(enst)
            continue

        if have_identical_cds(nm, enst):
            identical_cds.append(enst)

    return perfect_matches, identical_cds


def output_header(out):

    header = [
        'NM',
        'GENE_SYMBOL',
        'HGNC_ID',
        'ENST',
        'ENST_GENE_SYMBOL',
        'ENSG',
        'UTR_DIFF',
        'UTR_EXON_NUM_DIFF',
        'FLAG',
        'UTR_SELECTION',
        'PHASE1_DIFFERENCE_TYPE',
        'PHASE1_DECISIVE_CRITERIA'
    ]

    out.write('#' + '\t'.join(header) + '\n')


def output(out, nm, gene, ensts, flag, utr_selection_reason, phase1_difference_type, phase1_decisive_criteria):

    record = [nm.id_, gene, nm.gene_id]

    if ensts is not None:
        record += [
            ';'.join([x.id_ for x in ensts]),
            ';'.join([x.gene_symbol for x in ensts]),
            ';'.join([x.gene_id for x in ensts]),
            ';'.join([utr.utr_difference(nm, x) for x in ensts]),
            ';'.join([utr.utr_exon_number_difference(nm, x) for x in ensts])
        ]
    else:
        record += ['.'] * 5

    record += [flag, utr_selection_reason, phase1_difference_type, phase1_decisive_criteria]

    out.write('\t'.join(record) + '\n')


def are_perfect_match(t1, t2):

    cds_eq = t1.cds_exons == t2.cds_exons
    utr5_eq = t1.utr5_exons == t2.utr5_exons
    utr3_eq = t1.utr3_exons == t2.utr3_exons

    return cds_eq and utr5_eq and utr3_eq


def have_identical_cds(t1, t2):

    return t1.cds_exons == t2.cds_exons


def separator_line(txt):

    print '{} {} {}\n'.format('-' * 3, txt, '-' * (97 - len(txt)))
