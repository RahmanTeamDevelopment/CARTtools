import gzip


#######################################################################################################################

# Class representing an exon sequence
class ExonSequence(object):

    # Constructor
    def __init__(self, index, seq):
        self.index = index
        self.seq = seq


# Class representing a single Ensembl transcript
class EnsemblTranscript(object):

    # Constructor
    def __init__(self, line):
        self.exons = []
        cols = line.split('\t')
        cols = [cols[0], '.'] + cols[1:]
        self.ID = cols[0]
        self.gene = cols[1]
        self.geneID = cols[3]
        self.chrom = cols[5]
        self.strand = int(cols[6])
        self.transcriptStart = int(cols[7])
        self.transcriptEnd = int(cols[8])
        self.codingStart = int(cols[9])
        self.codingStartGenomic = int(cols[10])
        self.codingEndGenomic = int(cols[11])
        # Initializing and adding exons
        for i in range(1, len(cols) - 12, 2):
            self.exons.append(EnsemblExon(int((i + 1) / 2), int(cols[11 + i]), int(cols[12 + i])))

    # Return CDS, UTR5 and UTR3 dictionaries (keys: exon index, value: sequence)
    def getSeqs(self, ref):
        cds_seqs = []
        utr5_seqs = []
        utr3_seqs = []
        # int_seqs = dict()

        for i in range(len(self.exons)):
            exon = self.exons[i]

            if self.strand == 1:
                if self.codingStartGenomic > exon.end:
                    utr5_seqs.append(ExonSequence(i+1, ref.read(self.chrom, exon.start+1, exon.end)))
                    continue
                if self.codingEndGenomic < exon.start+1:
                    utr3_seqs.append(ExonSequence(i+1, ref.read(self.chrom, exon.start+1, exon.end)))
                    continue
                if exon.start+1 <= self.codingStartGenomic <= exon.end:
                    utr5_seqs.append(ExonSequence(i+1,ref.read(self.chrom, exon.start+1, self.codingStartGenomic-1)))
                    cdsstart = self.codingStartGenomic-1
                else:
                    cdsstart = exon.start
                if exon.start+1 <= self.codingEndGenomic <= exon.end:
                    utr3_seqs.append(ExonSequence(i+1, ref.read(self.chrom, self.codingEndGenomic+1, exon.end)))
                    cdsend = self.codingEndGenomic
                else:
                    cdsend = exon.end
                cds_seqs.append(ExonSequence(i+1, ref.read(self.chrom, cdsstart+1, cdsend)))

            else:
                if self.codingStartGenomic < exon.start+1:
                    utr5_seqs.append(ExonSequence(i+1, reverseComplement(ref.read(self.chrom, exon.start+1, exon.end))))
                    continue
                if self.codingEndGenomic > exon.end:
                    utr3_seqs.append(ExonSequence(i+1, reverseComplement(ref.read(self.chrom, exon.start+1, exon.end))))
                    continue
                if exon.start+1 <= self.codingEndGenomic <= exon.end:
                    utr3_seqs.append(ExonSequence(i+1, reverseComplement(ref.read(self.chrom, exon.start+1, self.codingEndGenomic-1))))
                    cdsstart = self.codingEndGenomic-1
                else:
                    cdsstart = exon.start
                if exon.start+1 <= self.codingStartGenomic <= exon.end:
                    utr5_seqs.append(ExonSequence(i+1, reverseComplement(ref.read(self.chrom, self.codingStartGenomic+1, exon.end))))
                    cdsend = self.codingStartGenomic
                else:
                    cdsend = exon.end
                cds_seqs.append(ExonSequence(i+1, reverseComplement(ref.read(self.chrom, cdsstart+1, cdsend))))

        return cds_seqs, utr5_seqs, utr3_seqs


# Class representing a single exon of an Ensembl transcript
class EnsemblExon(object):

    # Constructor
    def __init__(self, index, start, end):
        self.index = index
        self.start = start
        self.end = end
        self.length = end - start


#######################################################################################################################


# Return reverse complement sequence
def reverseComplement(seq):

    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "a": "t", "t": "a", "c": "g", "g": "c", "n": "n", '-': '-'}
    ret = ''
    for base in seq[::-1]:
        ret += complement[base]
    return ret


# Read Ensembl transcript database
def readEnsemblData(filename):
    ret = dict()
    n_transcripts = 0
    for line in gzip.open(filename,'r'):
        line = line.strip()
        if line=='' or line.startswith('#'): continue
        t = EnsemblTranscript(line)
        n_transcripts += 1
        ret[t.ID] = t

    return ret, n_transcripts


# Read Ensembl info txt file
def readEnsemblInfo(path):
    for line in open(path[:-3]+'.txt'):
        line = line.strip()
        if line == '': continue
        if line[0] == '#':
            return line[line.find('GRCh3'):line.find('GRCh3')+6], line[line.find('Ensembl release'):]


# Compare query transcript with a set of transcripts
def compare(g, transcript, comparewith, ref_GRCh37, ref_GRCh38, options, out):
    identical = False
    cds_identical = False
    for x in comparewith:
        flags = compareEnsemblVsEnsembl(transcript, x, ref_GRCh37, ref_GRCh38, options)

        if all([('CDS' not in flag) for flag in flags]): cds_identical = True

        if len(flags) == 0:
            identical = True
            if not options.discr: flags.append('.')
            else: continue

        out.write('\t'.join(['HGNC:' + g, transcript.ID, x.ID, ';'.join(flags)])+'\n')

    return identical, cds_identical


# Compare two Ensembl transcripts
def compareEnsemblVsEnsembl(enst_transcript1, enst_transcript2, ref_GRCh37, ref_GRCh38, options):

    flags = []

    # Extract CDS, UTR5 and UTR3 sequences for both releases using the appropriate genome build
    cds_1, utr5_1, utr3_1 = enst_transcript1.getSeqs(ref_GRCh37)
    cds_2, utr5_2, utr3_2 = enst_transcript2.getSeqs(ref_GRCh38)

    # Adding CDS flags
    if not len(cds_1) == len(cds_2):
        if options.simple: flags.append('CDS')
        else: flags.append('CDS:DN')
    else:
        discrv = []
        for i in range(len(cds_1)):
            if not cds_1[i].seq == cds_2[i].seq:
                discrv.append(str(cds_1[i].index))
        if len(discrv) > 0:
            if options.simple: flags.append('CDS')
            else: flags.append('CDS:'+','.join(discrv))

    # Adding UTR5 flags
    if not len(utr5_1) == len(utr5_2):
        if options.simple: flags.append('UTR5')
        else: flags.append('UTR5:DN')
    else:
        discrv = []
        for i in range(len(utr5_1)):
            if not utr5_1[i].seq == utr5_2[i].seq:
                if i == 0:
                    N = min(len(utr5_1[0].seq),len(utr5_2[0].seq))
                    if utr5_1[0].seq[-N:] == utr5_2[0].seq[-N:]: discrv.append('END')
                    else: discrv.append(str(utr5_1[0].index))
                else:
                    discrv.append(str(utr5_1[i].index))
        if len(discrv) > 0:
            if options.simple: flags.append('UTR5')
            else: flags.append('UTR5:'+','.join(discrv))

    # Adding UTR3 flags
    if not len(utr3_1) == len(utr3_2):
        if options.simple: flags.append('UTR3')
        else: flags.append('UTR3:DN')
    else:
        discrv = []
        for i in range(len(utr3_1)):
            if not utr3_1[i].seq == utr3_2[i].seq:
                if i == len(utr3_1)-1:
                    N = min(len(utr3_1[i].seq), len(utr3_2[i].seq))
                    if utr3_1[i].seq[:N] == utr3_2[i].seq[:N]: discrv.append('END')
                    else: discrv.append(str(utr3_1[i].index))
                else:
                    discrv.append(str(utr3_1[i].index))
        if len(discrv) > 0:
            if options.simple: flags.append('UTR3')
            else: flags.append('UTR3:'+','.join(discrv))

    return flags


def read_input_genes(fn):

    ret = []
    for line in open(fn):
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        ret.append(line[5:])
    return ret


def read_enst_file(fn):

    ret = {}
    for line in open(fn):
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        cols = line.split()
        if cols[3] == '.':
            continue
        ret[cols[2]] = cols[3]
    return ret

