from __future__ import division
from transcripts import TranscriptDBWriter
from urllib2 import urlopen
import urllib
import os
import sys
import refseq
import mapping
import urllib
import datetime
from main.version import __version__


def run(options):

    if options.mapping not in ['ucsc', 'ncbi']:
        print '\nAllowed values for option -m are \"ucsc\" or \"ncbi\".\n'
        sys.exit()

    if options.build not in ['GRCh37', 'GRCh38']:
        print '\nAllowed values for option -b are \"GRCh37\" or \"GRCh38\".\n'
        sys.exit()


    if options.mapping == 'ncbi':
        print '- RefSeq interim release, 2017-01-13 (incl. mapping) -\n'
    else:
        print '- Latest RefSeq release + UCSC mapping -\n'



    columns = ['ID', 'VERSION', 'HGNC_ID', 'INFO', 'STRAND', 'CHROM', 'START', 'END', 'EXONS', 'CODING_START',
               'CODING_END', 'SEQUENCE', 'CDNA_CODING_START', 'CDNA_CODING_END', 'EXON_CIGARS']

    tdb_writer = TranscriptDBWriter(options.output, source='refseq_db ' + __version__, build='GRCh37', columns=columns)

    # Initialize output files and write headers
    out_incl = open(options.output + '_included.txt', 'w')
    out_incl.write('\t'.join(['#ID', 'HGNCID']) + '\n')
    out_excl = open(options.output + '_excluded.txt', 'w')
    out_excl.write('\t'.join(['#ID', 'Reason']) + '\n')

    if options.mapping == 'ucsc':
        urls = refseq.access_refseq_files()
    else:
        if options.build == 'GRCh37':
            urls = [
                'ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/GRCh37.p13_interim_annotation/interim_GRCh37.p13_rna.gbk.gz']
        else:
            urls = [
                'ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.108/GRCh38.p10_interim_annotation/interim_GRCh38.p10_rna.gbk.gz']

    # Download and read mapping data
    sys.stdout.write('\rProcessing {} mapping data ... '.format(options.build))
    sys.stdout.flush()
    if options.build == 'GRCh37':
        if options.mapping == 'ucsc':
            mapping.download_ucsc_mapping('GRCh37', 'ucsc_mapping.gz')
            mappings = mapping.read_ucsc_mapping('ucsc_mapping.gz')
            os.remove('ucsc_mapping.gz')
        else:
            mapping.download_ncbi_mapping('GRCh37')
            mappings = mapping.read_ncbi_mapping('interim_GRCh37.p13_knownrefseq_alignments_2017-01-13.bam',
                                                 'interim_GRCh37.p13_top_level_2017-01-13.gff3.gz')
            os.remove('interim_GRCh37.p13_knownrefseq_alignments_2017-01-13.bam')
            os.remove('interim_GRCh37.p13_top_level_2017-01-13.gff3.gz')
        print '- done'
    else:
        if options.mapping == 'ucsc':
            mapping.download_ucsc_mapping('GRCh38', 'ucsc_mapping.gz')
            mappings = mapping.read_ucsc_mapping('ucsc_mapping.gz')
            os.remove('ucsc_mapping.gz')
        else:
            mapping.download_ncbi_mapping('GRCh38')
            mappings = mapping.read_ncbi_mapping('interim_GRCh38.p10_knownrefseq_alignments_2017-01-13.bam',
                                                 'interim_GRCh38.p10_top_level_2017-01-13.gff3.gz')
            os.remove('interim_GRCh38.p10_knownrefseq_alignments_2017-01-13.bam')
            os.remove('interim_GRCh38.p10_top_level_2017-01-13.gff3.gz')
        print '- done'

    # Iterate through available RefSeq data files
    counter_incl = 0
    counter_excl = 0
    for i in range(len(urls)):

        print_progress_info(i + 1, len(urls))

        # Download RefSeq data file and avoid timeout error
        while True:
            try:
                urllib.urlretrieve(urls[i], 'refseqdata.gz')
                break
            except:
                pass

        # Process RefSeq data file
        counter_incl, counter_excl = refseq.process_refseq_file('refseqdata.gz', mappings, tdb_writer, out_incl, out_excl, counter_incl, counter_excl)

        # Remove RefSeq data file
        os.remove('refseqdata.gz')

    print '- done'

    # Finalize transcript database
    tdb_writer.finalize()

    # Close output files
    out_incl.close()
    out_excl.close()

    print '\nSummary:'
    print ' - NMs included in the database: {}'.format(counter_incl)
    print ' - NMs excluded from the database: {}'.format(counter_excl)

    print '\nOutput files created:'
    print ' - {}.gz (+ .tbi)'.format(options.output)
    print ' - {}_included.txt'.format(options.output)
    print ' - {}_excluded.txt'.format(options.output)


def print_progress_info(counter, N):
    """Print out progress information"""

    sys.stdout.write('\rProcessing RefSeq data (file {}/{}) ... '.format(counter, N))
    sys.stdout.flush()

