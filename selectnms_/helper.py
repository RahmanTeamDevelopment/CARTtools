import os
from transcripts import Transcript, TranscriptDB
import utr


def read_excluded_transcripts(fn):
    """Read excluded transcripts from txt file"""

    ret = dict()
    for line in open(fn):
        line = line.strip()
        if line.startswith('#'):
            continue
        cols = line.split()
        ret[cols[0]] = cols[1]
    return ret


def readInputData(options):
    """Read input datasets"""

    # Read in gene name dictionaries
    gene_dict = read_gene_dict(options.genes)

    # Read in APPRIS data
    appris, appris_principal = read_appris_file(options.appr)

    # Read in RefSeq transcript database
    refseq_db = TranscriptDB(options.refsdb)
    refseq_db.read()

    # Read excluded transcripts from file
    excluded = read_excluded_transcripts(options.refsdb[:-3] + '_excluded.txt')

    # Read in refseq_scan output
    refseqscan = read_refseqscan_output(options.refss)

    return gene_dict, appris, appris_principal, refseq_db, excluded, refseqscan


def convert_file_names_to_absolute_path(options):
    """Convert file names specified to their absolute path"""

    if '/' in options.hgnc: options.hgnc = os.path.abspath(options.hgnc)
    if '/' in options.appr: options.appr = os.path.abspath(options.appr)
    if '/' in options.refsdb: options.refsdb = os.path.abspath(options.refsdb)
    if '/' in options.refss: options.refss = os.path.abspath(options.refss)
    return options


def initialize_output_files(options):
    """Initialize output files"""

    # Output files
    selected = open(options.out + '_auto_selected.txt', 'w')
    missing = open(options.out + '_auto_missing.txt', 'w')
    log = open(options.out + '_auto_log.txt', 'w')

    # Write headers
    selected.write('#' + '\t'.join(
        ['HGNCID',
         'RelatedTranscript',
         'Version',
         'PRINCIPAL',
         'UTRDifferenceType',
         'UTRDecisiveCriteria',
         'GenomeDifference',
         'Alignment'
         ]
    ) + '\n')

    missing.write('#' + '\t'.join(['HGNCID', 'Reason']) + '\n')
    return selected, missing, log


def to_missing_list(missing, log, hgncid, missing_msg, log_msg):
    """Output to missing list"""

    missing.write(hgncid + '\t' + missing_msg + '\n')
    log.write(log_msg + '\n')


def translate_gene_id(gene_dict, hgncid, missing, log):
    """Translate HGNC ID to NCBI ID"""

    ncbi_geneid = gene_dict.get(hgncid)
    if ncbi_geneid is None:
        to_missing_list(missing, log, hgncid, 'no_NCBI_GeneID', '\nThe gene has no NCBI Gene ID\n\nGene added to Missing List')
        return
    log.write('-> NCBI Gene ID: ' + ncbi_geneid + '\n')
    log.write('\n')
    return ncbi_geneid


def get_appris_principal_isoforms_for_gene(appris, ncbi_geneid, missing, log, hgncid):
    """Get APPRIS PRINCIPAL isoforms for gene"""

    if ncbi_geneid not in appris:
        to_missing_list(missing, log, hgncid, 'not_in_APPRIS', '\nThe gene has no APPRIS principal isoforms \n\nGene added to Missing List')
        return
    log.write('The gene has the following APPRIS principal isoforms: ' + str(appris[ncbi_geneid]) + '\n\n')
    return appris[ncbi_geneid]


def get_nms(appris_principal, missing, log, hgncid):
    """Get only NM principal isoform transcripts (XMs are excluded)"""

    nms = [id for id in appris_principal if id.startswith('NM_')]
    if len(nms) == 0:
        to_missing_list(missing, log, hgncid, 'only_XM', 'APPRIS contains only XM principal isoforms for this gene; gene added to Missing List')
        return
    return nms


def check_nms_in_refseq_db(refseq_db, nms, missing, log, hgncid):
    """Get NMs that are in the RefSeq db"""

    nms_in_refseq = [id for id in nms if refseq_db.contains(id)]
    if len(nms_in_refseq) == 0:
        to_missing_list(missing, log, hgncid, 'not_in_RefSeq', 'None of the APPRIS principal NMs for this gene are found in RefSeq db; gene added to Missing List')
        return
    log.write(str(len(nms_in_refseq)) + ' APPRIS principal NMs present in RefSeq db: ' + str(nms_in_refseq) + '\n\n')
    return nms_in_refseq


def create_transcript_objects(refseq_db, nms_in_refseq):
    """Retrieve Transcript objects from RefSeq db"""

    transcripts = []
    for id in nms_in_refseq:
        t = refseq_db.by_id(id)
        transcripts.append(t)
    return transcripts


def get_transcripts_with_consistent_hgncid(transcripts, hgncid, missing, log):
    """Get only transcripts that have consistent HGNC ID"""

    transcripts = [t for t in transcripts if t.hgnc_id == hgncid]
    if len(transcripts) == 0:
        to_missing_list(missing, log, hgncid, 'RefSeq_gene_mismatch', 'All APPRIS principal NMs that are present in RefSeq db have discrepant HGNC:ID; gene added to Missing List')
        return
    log.write(str(len(transcripts)) + ' NMs to select from: ' + str([t.id for t in transcripts]) + '\n')
    return transcripts


def output_single_transcript(transcripts, refseqscan, appris_principal, selected, log, hgncid):
    """Output single selected transcript"""

    log.write('\nSELECTED TRANSCRIPT: ' + transcripts[0].id + ' (The only transcript to select)\n')
    diff = refseqscan[transcripts[0].id]

    if transcripts[0].exon_cigars == '.':
        cigar_str = '.'
    else:
        cigar_info = []
        tmp = transcripts[0].exon_cigars.split(',')
        for i in range(len(tmp)):
            if 'X' in tmp[i] or 'I' in tmp[i] or 'D' in tmp[i]:
                cigar_info.append('Ex{}:{}'.format(i + 1, tmp[i]))
        cigar_str = ','.join(cigar_info) if len(cigar_info) > 0 else '.'

    log.write('Added.\n')
    selected.write('\t'.join([hgncid, transcripts[0].id, transcripts[0].version, appris_principal[transcripts[0].id], 'not_required', 'not_required', diff, cigar_str]) + '\n')


def select_from_multiple_candidates(transcripts, refseqscan, appris_principal, selected, missing, log, hgncid):
    """Select from multiple transcript candidates"""

    if not all_have_same_cds(transcripts):
        to_missing_list(missing, log, hgncid, 'coding_mismatch', 'Not all of these transcripts have the same CDS;  gene added to Missing List')
        return

    sel, differenceType, decisiveCriteria = utr.utr_selection(transcripts, log)

    if sel is None:
        to_missing_list(missing, log, hgncid, 'non_unique_mapping', 'Candidate templates have identical mapping; gene added to Missing List')
        return

    diff = refseqscan[sel.id]

    if sel.exon_cigars == '.':
        cigar_str = '.'
    else:
        cigar_info = []
        tmp = sel.exon_cigars.split(',')
        for i in range(len(tmp)):
            if 'X' in tmp[i] or 'I' in tmp[i] or 'D' in tmp[i]:
                cigar_info.append('Ex{}:{}'.format(i + 1, tmp[i]))
        cigar_str = ','.join(cigar_info) if len(cigar_info) > 0 else '.'


    log.write('SELECTED TRANSCRIPT: ' + sel.id + ' (Selected by UTR criteria)\n')

    log.write('Added.\n')

    selected.write('\t'.join([hgncid, sel.id, sel.version, appris_principal[sel.id], differenceType, decisiveCriteria, diff, cigar_str]) + '\n')


def all_have_same_cds(transcripts):
    """Check if all transcripts have the same CDS (same CDS exon boundaries)"""

    for i in range(1, len(transcripts)):
        if not same_cds(transcripts[0], transcripts[i]):
            return False
    return True


def same_cds(transcript1, transcript2):
    """Check if two transcripts have the same CDS (same CDS exon boundaries)"""

    cds_exons1 = transcript1.cds_regions()
    cds_exons2 = transcript2.cds_regions()

    if len(cds_exons1) != len(cds_exons2):
        return False
    for i in range(len(cds_exons1)):
        exon1 = cds_exons1[i]
        exon2 = cds_exons2[i]
        if exon1[0] != exon2[0] or exon1[1] != exon2[1]:
            return False
    return True


def read_appris_file(fn):
    """Read file downloaded from the APPRIS website"""

    ret = dict()
    ret_principal = dict()
    for line in open(fn):
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        cols = line.split()
        if 'PRINCIPAL' not in cols[4]:
            continue
        if cols[1] not in ret:
            ret[cols[1]] = []
        id = cols[2]
        if '.' in id:
            id = id[:id.find('.')]
        ret[cols[1]].append(id)
        ret_principal[id] = cols[4]
    return ret, ret_principal


def read_refseqscan_output(fn):
    """Read RefSeqScan output file"""

    ret = dict()
    for line in open(fn):
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        cols = line.split()
        ret[cols[0]] = cols[2]
    return ret


def read_gene_dict(fn):
    """Read gene name dictionary file"""

    ret_NCBIGeneID = dict()
    for line in open(fn):
        line = line.strip()
        if line == '' or line.startswith('#') or line.startswith('HGNC ID'):
            continue
        cols = line.split()
        hgnc = ncbi = ''
        for x in cols:
            x = x.strip()
            if x.startswith('HGNC:'):
                hgnc = x
            else:
                ncbi = x
        if ncbi != '':
            ret_NCBIGeneID[hgnc] = ncbi
    return ret_NCBIGeneID
