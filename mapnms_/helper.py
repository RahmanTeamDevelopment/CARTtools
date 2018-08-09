from __future__ import division

from operator import itemgetter
import pysam
import datetime
import os

from main.version import __version__
from transcripts import TranscriptDB


allowed_chroms = map(str, range(1, 24)) + ['X', 'Y', 'MT']

class TranscriptDBWriter(object):
    """Class for creating new transcript database"""

    def __init__(self, fn, source='', build='', columns=[]):
        """Constructor of the TranscriptDBWriter class"""

        self._fn = fn
        self._source = source
        self._build = build
        self._columns = [x.lower() for x in columns]
        self._records = {c: [] for c in allowed_chroms}
        self.idx_chrom = self._columns.index('chrom')
        self.idx_start = self._columns.index('start')
        self.idx_end = self._columns.index('end')

    def add(self, transcript):
        """Add transcript to DB"""

        record = []
        for c in self._columns:
            if c in ['exons', 'cdna_exons']:
                record.append(','.join([str(e.start) + '-' + str(e.end) for e in getattr(transcript, c.lower())]))
            elif c in ['start', 'end', 'coding_start', 'coding_end', 'cdna_coding_start', 'cdna_coding_end']:
                record.append(int(getattr(transcript, c.lower())))
            else:
                record.append(str(getattr(transcript, c.lower())))
        self._records[transcript.chrom].append(record)

    def _sort_records(self):
        """Sort records by chrom, start, end"""

        idx_start = self._columns.index('start')
        idx_end = self._columns.index('end')
        for c in allowed_chroms:
            if c in self._records:
                self._records[c] = sorted(self._records[c], key=itemgetter(idx_start, idx_end))

    def _index_with_tabix(self):
        """Compress and index output file by Tabix"""

        pysam.tabix_compress(self._fn + '_tmp', self._fn + '.gz', force=True)
        pysam.tabix_index(self._fn + '.gz', seq_col=4, start_col=6, end_col=7, meta_char='#', force=True)


    def finalize(self, options):
        """Write to file, compress and index, clean up"""

        # Sort records by CHROM, START, END
        self._sort_records()

        # Initialize file and write header
        out = open(self._fn + '_tmp', 'w')
        out.write('#date: ' + str(datetime.datetime.today()).split()[0] + '\n')
        out.write('#build: ' + self._build + '\n')
        out.write('#ncbi_source_db: ' + options.ncbi + '\n')
        out.write('#ucsc_source_db: ' + options.ucsc + '\n')
        out.write('#hgnc_biomart_file: ' + options.hgnc + '\n')

        # Write records to file
        for c in allowed_chroms:
            if c in self._records:
                for record in self._records[c]:
                    record = map(str, record)

                    old_record = [
                        record[0],
                        record[2],
                        record[1],
                        record[3],
                        record[5],
                        '1' if record[4] == '+' else '-1',
                        record[6],
                        record[7],
                        record[-2],
                        str(int(record[9]) + 1),
                        str(int(record[10]) + 1)
                    ]

                    for e in record[8].split(','):
                        [start, end] = e.split('-')
                        old_record.append(start)
                        old_record.append(end)

                    out.write('\t'.join(old_record) + '\n')

        out.close()

        # Compress and index by Tabix
        self._index_with_tabix()

        # Remove temporary file
        os.remove(self._fn + '_tmp')


def check_for_missing_hgnc_ids(fn, db_ncbi, db_ucsc, genes_symbols, symbols):
    not_found = []
    for line in open(fn):
        line = line.strip()
        if line=='' or line.startswith('#'):
            continue
        nm = line.split()[1]

        if db_ncbi.contains(nm):
            transcript = db_ncbi.by_id(nm)
        elif db_ucsc.contains(nm):
            transcript = db_ucsc.by_id(nm)
        else:
            continue

        if transcript.hgnc_id not in genes_symbols and transcript.hgnc_id not in symbols:
            not_found.append(transcript.hgnc_id)

    if len(not_found) > 0:
        print '----------------------------------------------------------------------------------------'
        print 'Oops. The following {} HGNC ID{} not found in HGNC BioMart file:'.format(len(not_found), ' is' if len(not_found) == 1 else 's are')
        for x in not_found:
            print '- {}'.format(x)
        print 'Use command line option -s to supply a txt file providing gene {}.'.format('symbol for this HGNC ID' if len(not_found) == 1 else 'symbols for these HGNC IDs')
        print 'The txt file must have two columns: HGNC ID and gene symbol.'
        print '----------------------------------------------------------------------------------------\n'
        print 'No outputs generated!'
        print '=' * 100 + '\n'
        quit()


def read_gene_symbol_file(fn):
    ret = {}
    for line in open(fn):
        line = line.strip()
        if line=='' or line.startswith('#'):
            continue
        cols = line.split()
        ret[cols[0]] = cols[1]
    return ret


def number_of_input_records(fn):

    ret = 0
    for line in open(fn):
        line = line.strip()
        if line=='' or line.startswith('#'):
            continue
        ret += 1
    return ret


def read_excluded_list(fn):

    ret = {}
    for line in open(fn):
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        cols = line.split()
        ret[cols[0]] = cols[1]
    return ret


def read_database(fn, txt):

    ret = TranscriptDB(fn)
    ret.read()
    ret_excluded = read_excluded_list(fn[:-3] + '_excluded.txt')
    print 'Transcript database ({} mapping): {} ({} transcripts)'.format(txt, fn, len(ret._data))
    return ret, ret_excluded


def initialize_transcript_db_writer(options):

    columns = [
        'ID',
        'HGNC_ID',
        'GENE_SYMBOL',
        'INFO',
        'STRAND',
        'CHROM',
        'START',
        'END',
        'EXONS',
        'CODING_START',
        'CODING_END',
        'CDNA_CODING_START',
        'CDNA_CODING_END'
    ]
    tdb_writer = TranscriptDBWriter(
        options.output,
        source='?',
        build='?',
        columns=columns
    )
    return tdb_writer


def print_summary_info(options, counter, counter_ncbi, counter_ucsc, counter_missing):

    print 'done'
    print '\nSummary:'
    print '{} mapped NMs added to output database ({} with NCBI mapping, {} with UCSC mapping)'.format(
        counter,
        counter_ncbi,
        counter_ucsc
    )
    print '{} NMs missing from output database'.format(counter_missing)

    print '\nOutput files:'
    print '----------------------------------------------------'
    print ' - Mapped NM database: {}.gz (+.tbi)'.format(options.output)
    print ' - Mapping source file: {}_source.txt'.format(options.output)
    print ' - Missing NMs file: {}_missing.txt'.format(options.output)
    print '----------------------------------------------------'


def create_transcript_object(db_ncbi, db_ucsc, nm, out_source, counter_ncbi, counter_ucsc, db_ncbi_excluded,
                             db_ucsc_excluded, out_missing, counter_missing, genes_symbols, symbols):

    if db_ncbi.contains(nm):
        transcript = db_ncbi.by_id(nm)
        out_source.write('{}\tNCBI\n'.format(nm))
        counter_ncbi += 1

    elif db_ucsc.contains(nm):
        transcript = db_ucsc.by_id(nm)
        out_source.write('{}\tUCSC\n'.format(nm))
        counter_ucsc += 1

    else:
        issues = [nm]
        if nm in db_ncbi_excluded:
            issues.append(db_ncbi_excluded[nm])
        else:
            issues.append('not_found')
        if nm in db_ucsc_excluded:
            issues.append(db_ucsc_excluded[nm])
        else:
            issues.append('not_found')
        out_missing.write('\t'.join(issues) + '\n')
        counter_missing += 1
        return None, counter_ncbi, counter_ucsc, counter_missing

    transcript.id = nm
    set_gene_symbol(transcript, genes_symbols, symbols)
    transcript.hgnc_id = transcript.hgnc_id[5:]

    return transcript, counter_ncbi, counter_ucsc, counter_missing


def set_gene_symbol(transcript, genes_symbols, symbols):

    if transcript.hgnc_id in genes_symbols:
        transcript.gene_symbol = genes_symbols[transcript.hgnc_id]
    elif transcript.hgnc_id in symbols:
        transcript.gene_symbol = symbols[transcript.hgnc_id]
    else:
        transcript.gene_symbol = 'hahahahaha'
        return

        # This should never happen
        print '\nError 1\n'
        quit()


def read_gene_symbols(options):

    # Reading gene symbol file
    genes_symbols = read_gene_symbol_file(options.hgnc)
    print 'HGNC BioMart file: {}'.format(options.hgnc)

    # Reading additional gene symbol file
    if options.symbols:
        symbols = read_gene_symbol_file(options.symbols)
        print 'Txt file supplying missing gene symbols: {}'.format(options.symbols)
    else:
        symbols = {}

    return genes_symbols, symbols
