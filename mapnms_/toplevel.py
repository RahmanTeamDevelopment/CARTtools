import helper
import sys


def count_input_nms(fn):

    count = 0
    for line in open(fn):
        if line[0] == '#':
            continue
        line = line.strip()
        cols = line.split('\t')
        if cols[-1] != '.':
            count += 1
    return count


def run(options):

    num_of_nms = count_input_nms(options.input)

    print 'Input file: {} ({} transcripts)\n'.format(options.input, num_of_nms)

    # Read gene symbols
    genes_symbols, symbols = helper.read_gene_symbols(options)

    # Initialize transcript database writer
    tdb_writer = helper.initialize_transcript_db_writer(options)

    # Read NCBI and UCSC databases
    db_ncbi, db_ncbi_excluded = helper.read_database(options.ncbi, 'NCBI')
    db_ucsc, db_ucsc_excluded = helper.read_database(options.ucsc, 'UCSC')

    # Check for missing HGNC IDs
    #helper.check_for_missing_hgnc_ids(options.input, db_ncbi, db_ucsc, genes_symbols, symbols)


    out_source = open(options.output + '_source.txt', 'w')
    out_source.write('#NM\tsource_db\n')

    out_missing = open(options.output + '_missing.txt', 'w')
    out_missing.write('#NM\treason_NCBI_db\treason_UCSC_db\n')


    # Iterating through input records

    sys.stdout.write('\nProcessing data ... ')
    sys.stdout.flush()
    counter = counter_ncbi = counter_ucsc = counter_missing = 0

    for line in open(options.input):
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue

        cols = line.split()
        nm = cols[-1]

        if nm == '.':
            continue

        # Create transcript object and write to _source and _missing output files
        transcript, counter_ncbi, counter_ucsc, counter_missing = helper.create_transcript_object(
            db_ncbi,
            db_ucsc,
            nm,
            out_source,
            counter_ncbi,
            counter_ucsc,
            db_ncbi_excluded,
            db_ucsc_excluded,
            out_missing,
            counter_missing,
            genes_symbols,
            symbols
        )

        # If nm is missing from the databases
        if transcript is None:
            continue

        # Adding transcript to database writer
        tdb_writer.add(transcript)

        counter += 1

    tdb_writer.finalize(options)

    out_source.close()
    out_missing.close()

    # Print out summary info
    helper.print_summary_info(options, counter, counter_ncbi, counter_ucsc, counter_missing)
