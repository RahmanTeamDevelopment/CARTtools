from __future__ import division
import sys
import pysam
import parser
import helper
from utr import utr_selection
from main.version import __version__


def read_data(options):

    # Read CARTs data
    sys.stdout.write('Reading CARTs data ... ')
    sys.stdout.flush()
    cart_data = parser.read_transcript_db(options.cart_data)
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

    return cart_data, ensembl_tabix_file, ensembl_data, gene_synonyms


def run(options):

    print '\n'
    helper.separator_line('cart2enst {} running'.format(__version__))

    # Initialize output files
    out = open(options.output + '.txt', 'w')
    helper.output_header(out)
    log = open(options.output + '_phase1_log.txt', 'w')

    # Read CART and Ensembl data files
    cart_data, ensembl_tabix_file, ensembl_data, gene_synonyms = read_data(options)

    if options.input:
        cart_list = parser.read_input_file(options.input)

    sys.stdout.write('\nProcessing CARTs... ')
    sys.stdout.flush()

    # Iterate through the CARTs
    for gene_id, carts in cart_data.iteritems():
        for cart in carts:

            if options.input:
                if cart.id_ not in cart_list:
                    continue

            # Find perfect matches or ENSTs with identical CDS
            if options.any_gene:
                perfect_matches, identical_cds = helper.perfect_and_cds_matches_any_gene(
                    ensembl_tabix_file,
                    cart,
                    gene_id,
                    out
                )
            else:
                perfect_matches, identical_cds = helper.perfect_and_cds_matches_synonyms(
                    gene_synonyms,
                    ensembl_data,
                    cart,
                    gene_id,
                    out
                )

            if perfect_matches is None:
                continue

            # If there are multiple perfect matches
            if len(perfect_matches) > 1:
                helper.output(
                    out,
                    cart,
                    gene_id,
                    perfect_matches,
                    'multiple_exact_matches',
                    '.',
                    '.',
                    '.'
                )
                continue

            # If there is only one perfect match
            if len(perfect_matches) == 1:
                helper.output(
                    out,
                    cart,
                    gene_id,
                    [perfect_matches[0]],
                    'single_exact_match',
                    '.',
                    '.',
                    '.'
                )
                continue

            # If there is no ENST with identical CDS
            if len(identical_cds) == 0:
                helper.output(
                    out,
                    cart,
                    gene_id,
                    None,
                    'no_cds_match',
                    '.',
                    '.',
                    '.'
                )
                continue

            # If there is only one ENST with identical CDS
            if len(identical_cds) == 1:
                helper.output(
                    out,
                    cart,
                    gene_id,
                    [identical_cds[0]],
                    'single_cds_match',
                    '.',
                    '.',
                    '.'
                )
                continue

            # If there are multiple ENSTs with identical CDS, apply UTR selection
            selected, reason, phase1_difference_type, phase1_decisive_criteria = utr_selection(cart, identical_cds, log)

            helper.output(
                out,
                cart,
                gene_id,
                selected,
                'multiple_cds_matches',
                reason,
                phase1_difference_type,
                phase1_decisive_criteria
            )

    print 'Done.'

    # Close output files
    out.close()
    log.close()

    print ''
    helper.separator_line('cart2enst finished')
