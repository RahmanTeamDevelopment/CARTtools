import sys
import datetime
from transcripts import TranscriptDB
from reference import Reference
from main.version import __version__


def run(options):

    if options.input is None:
        sys.exit('\nInput RefSeq database file not specified.')
    if options.reference is None:
        sys.exit('\nReference genome file not specified.')

    # Welcome message
    print '\n' + '=' * 100
    print 'RefSeqScan v' + __version__ + ' started: ' + str(datetime.datetime.now()) + '\n'

    # Print reference genome name
    print 'RefSeq transcript db: ' + options.input
    print 'Reference genome: ' + options.reference
    print ''

    # Initialize output file
    outfile = open(options.output, 'w')
    outfile.write(
        '#Created: ' + datetime.datetime.now().strftime("%d-%m-%Y") + '; Reference genome: ' + options.reference + '\n')
    outfile.write('#' + '\t'.join(['ID', 'VERSION', 'DIFFS']) + '\n')

    # Transcript database
    tdb = TranscriptDB(options.input)

    # Initialize Reference object
    reference = Reference(options.reference)

    sys.stdout.write('Processing transcripts ... ')
    sys.stdout.flush()

    # Iterate through transcripts
    for transcript in tdb.generator():
        # Retrieve reference sequence corresponding to transcript
        reference_sequence = get_reference_sequence(transcript, reference)

        # Compare transcript sequence with reference sequence
        diffs = compare_sequences(reference_sequence, transcript.sequence, transcript.cdna_coding_start,
                                  transcript.cdna_coding_end)

        # Write results to output file
        outfile.write('\t'.join([transcript.id, transcript.version, diffs]) + '\n')

    print '- Done'

    # Close output files
    outfile.close()

    # Print output file name
    print '\nOutput file (' + options.output + ') created'

    # Goodbye message
    print '\nFinished: ' + str(datetime.datetime.now())
    print '=' * 100 + '\n'


def reverse_complement(s):
    """Return reverse complement sequence"""

    ret = ''
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "a": "t", "t": "a", "c": "g", "g": "c", "n": "n"}
    for base in s[::-1]:
        ret += complement[base]
    return ret


def get_reference_sequence(transcript, reference):
    """Return reference cDNA sequence of transcript"""

    ret = ''
    for exon in transcript.exons:
        exon_seq = reference.get_sequence(transcript.chrom, exon.start, exon.end)
        if transcript.strand == '-':
            exon_seq = reverse_complement(exon_seq)
        ret += exon_seq
    return ret


def compare_sequences(ref, seq, cdna_coding_start, cdna_coding_end):
    """Compare transcript sequence with reference"""

    if len(ref) != len(seq):
        while seq[-1] == 'A':
            seq = seq[:-1]
        while ref[-1] == 'A':
            ref = ref[:-1]

    if ref == seq: return '.'

    if len(ref) != len(seq):
        return 'lengthChange'

    ret = []
    for i in range(len(ref)):
        if ref[i] != seq[i]:
            pos = i + 1
            if pos < cdna_coding_start:
                coord = str(pos - cdna_coding_start)
            elif pos > cdna_coding_end:
                coord = '+' + str(pos - cdna_coding_end)
            else:
                coord = str(pos - cdna_coding_start + 1)
            ret.append('c.' + coord + ref[i] + '>' + seq[i])

    return ','.join(ret)



