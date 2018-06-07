import helper
import reference
import sys


def run(options):

    print 'Input file: {} -> {} CARTs\n'.format(options.input, helper.number_of_input_carts(options.input))

    # Read gene symbols
    genes_symbols, symbols = helper.read_gene_symbols(options)

    # Initialize reference sequence reader
    ref = reference.Reference(options.ref)
    print 'Reference genome file: {}\n'.format(options.ref)

    # Initialize transcript database writer
    tdb_writer = helper.initialize_transcript_db_writer(options)

    # Read NCBI and UCSC databases
    db_ncbi, db_ncbi_excluded = helper.read_database(options.ncbi, 'NCBI')
    db_ucsc, db_ucsc_excluded = helper.read_database(options.ucsc, 'UCSC')

    # Check for missing HGNC IDs
    helper.check_for_missing_hgnc_ids(options.input, db_ncbi, db_ucsc, genes_symbols, symbols)

    # Initialize output files
    out_source, out_missing, out_genepred, out_fasta, \
    out_genepred_annovar, out_fasta_annovar, gbk_dir = helper.initialize_output_files(options)

    # Iterating through input records

    sys.stdout.write('Processing data ... ')
    sys.stdout.flush()
    counter = counter_ncbi = counter_ucsc = counter_missing = 0
    gff3_lines = {}

    for line in open(options.input):
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue

        cols = line.split()
        cart_id = cols[0]
        nm = cols[1]

        # Create transcript object and write to _source and _missing output files
        transcript, counter_ncbi, counter_ucsc, counter_missing = helper.create_transcript_object(
            db_ncbi,
            db_ucsc,
            nm,
            out_source,
            counter_ncbi,
            counter_ucsc,
            cart_id,
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

        # Creating content of gff3 file
        gff3_lines = helper.create_gff3_lines(transcript, gff3_lines)

        # Writing to gp file
        helper.output_genepred(transcript, out_genepred)

        # Writing to gbk output
        if options.gbk:
            helper.output_gbk(transcript, ref, gbk_dir)

        # Writing to fasta file
        helper.output_fasta(transcript, out_fasta, ref)

        # Writing annovar files
        if options.annovar:
            helper.output_genepred(transcript, out_genepred_annovar)
            helper.output_fasta_annovar(transcript, out_fasta_annovar, ref)

        counter += 1

    # Creating bgzipped, Tabix-index GFF3 output
    helper.output_gff3(gff3_lines, options.output + '.gff')

    # Finalize outputs
    helper.finalize_outputs(
        options,
        tdb_writer,
        out_source,
        out_missing,
        out_fasta,
        out_genepred,
        out_genepred_annovar,
        out_fasta_annovar,
        gbk_dir
    )

    # Print out summary info
    helper.print_summary_info(options, counter, counter_ncbi, counter_ucsc, counter_missing)
