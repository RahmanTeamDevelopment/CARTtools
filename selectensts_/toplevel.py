from __future__ import division
import sys
import pysam
import parser
import helper
from utr import utr_selection
from main.version import __version__


def read_data(options):

    # Read mapped NM data
    sys.stdout.write('Reading mapped NM data ... ')
    sys.stdout.flush()
    NMs_data = parser.read_transcript_db(options.mapped_nms)
    print 'Done.'

    # Connect to Ensembl tabix file or read Ensembl data into memory
    sys.stdout.write('Reading Ensembl data ... ')
    sys.stdout.flush()
    if options.any_gene:
        ensembl_tabix_file = pysam.Tabixfile(options.ensembl_data)
        ensembl_data = None
        gene_synonyms = None
    else:
        ensembl_tabix_file = None
        ensembl_data = parser.read_transcript_db(options.ensembl_data)
        gene_synonyms = parser.read_gene_synonyms(options.gene_synonyms)
    print 'Done.'

    return NMs_data, ensembl_tabix_file, ensembl_data, gene_synonyms


def run(options):

    # Initialize output files
    out = open(options.output + '.txt', 'w')
    helper.output_header(out)
    log = open(options.output + '_phase1_log.txt', 'w')

    # Read NMs and Ensembl data files
    NMs_data, ensembl_tabix_file, ensembl_data, gene_synonyms = read_data(options)

    if options.input:
        nm_list = parser.read_input_file(options.input)

    sys.stdout.write('\nProcessing NMs ... ')
    sys.stdout.flush()

    # Iterate through the NMs
    counter = 0
    counter_selected = 0
    counter_problem = 0
    for gene, NMs in NMs_data.iteritems():
        for nm in NMs:

            if options.input:
                if nm.id_ not in nm_list:
                    continue

            # Find perfect matches or ENSTs with identical CDS
            if options.any_gene:
                perfect_matches, identical_cds = helper.perfect_and_cds_matches_any_gene(
                    ensembl_tabix_file,
                    nm,
                    gene,
                    out
                )
            else:
                perfect_matches, identical_cds = helper.perfect_and_cds_matches_synonyms(
                    gene_synonyms,
                    ensembl_data,
                    nm,
                    gene,
                    out
                )

            if perfect_matches is None:
                counter += 1
                counter_problem += 1
                continue

            # If there are multiple perfect matches
            if len(perfect_matches) > 1:
                helper.output(
                    out,
                    nm,
                    gene,
                    perfect_matches,
                    'multiple_exact_matches',
                    '.',
                    '.',
                    '.'
                )
                counter += 1
                counter_selected += 1
                continue

            # If there is only one perfect match
            if len(perfect_matches) == 1:
                helper.output(
                    out,
                    nm,
                    gene,
                    [perfect_matches[0]],
                    'single_exact_match',
                    '.',
                    '.',
                    '.'
                )
                counter += 1
                counter_selected += 1
                continue

            # If there is no ENST with identical CDS
            if len(identical_cds) == 0:
                helper.output(
                    out,
                    nm,
                    gene,
                    None,
                    'no_cds_match',
                    '.',
                    '.',
                    '.'
                )
                counter += 1
                counter_problem += 1
                continue

            # If there is only one ENST with identical CDS
            if len(identical_cds) == 1:
                helper.output(
                    out,
                    nm,
                    gene,
                    [identical_cds[0]],
                    'single_cds_match',
                    '.',
                    '.',
                    '.'
                )
                counter += 1
                counter_selected += 1
                continue

            # If there are multiple ENSTs with identical CDS, apply UTR selection
            selected, reason, phase1_difference_type, phase1_decisive_criteria = utr_selection(nm, identical_cds, log)

            helper.output(
                out,
                nm,
                gene,
                selected,
                'multiple_cds_matches',
                reason,
                phase1_difference_type,
                phase1_decisive_criteria
            )
            counter += 1
            counter_selected += 1

    print 'Done.'

    print '\nA total of {} NMs have been processed (-> {})'.format(counter, options.output + '.txt')
    print 'ENST selected for {} NMs, no selection made for {} NMs'.format(counter_selected, counter_problem)

    # Close output files
    out.close()
    log.close()

