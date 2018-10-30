import helper
import reference
import sys
import transcripts
from main.version import __version__


def run(options):

    if not ((options.series.startswith('CART37') or options.series.startswith('CART38')) and len(options.series) == 7):
        print '\nSeries code incorrect!\n'
        quit()

    # ...
    selected_ensts = helper.read_selected_ensts(options.selected_ensts)

    # ...
    canonical_ensts = helper.read_canonical_ensts(options.canonical)

    # Initialize reference sequence reader
    ref = reference.Reference(options.ref)

    # Initialize transcript database writer
    tdb_writer = helper.initialize_transcript_db_writer(options)

    # Read Ensembl database
    ensembl_db = transcripts.read_ensembl_db(options.ensembl)
    ensembl_by_symbol = transcripts.read_ensembl_db_by_symbol(options.ensembl)

    # Read previous CAVA db output and reference genome if required
    if options.prev_cava_db:
        prev_ref = reference.Reference(options.prev_ref)
        prev_cava_db = helper.read_prev_cava_db(options.prev_cava_db, prev_ref)
    else:
        prev_cava_db = None

    # Initialize output files
    out_genepred, out_fasta, out_genepred_annovar, out_fasta_annovar, gbk_dir, out_id, out_excl = helper.initialize_output_files(options)

    # Initialize progress info
    sys.stdout.write('Processing {} genes ... '.format(helper.number_of_genes(options.selected_nms)))
    sys.stdout.flush()

    # Initialize CART numbering
    cartidx = 10000 if options.prev_cava_db is None else helper.get_last_cartidx(options.prev_cava_db)

    # Iterate through input records
    count_excluded = 0
    count_selected = 0
    count_canonical_or_longest = 0
    gff2_lines = {}
    gff3_lines = {}
    for line in open(options.selected_nms):
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        cols = line.split()

        symbol = cols[0]
        hgnc_id = cols[1]
        assoc_nm = cols[-1]

        if hgnc_id in selected_ensts and selected_ensts[hgnc_id] != '.':
            enst = selected_ensts[hgnc_id]
            selected = True
        elif symbol in canonical_ensts and canonical_ensts[symbol] in ensembl_db:
            enst = canonical_ensts[symbol]
            selected = False
        elif symbol in ensembl_by_symbol:
            enst = transcripts.find_longest_transcript(ensembl_by_symbol[symbol]).id
            selected = False
        else:
            out_excl.write('{}\t{}\t{}\n'.format(hgnc_id, symbol, 'no_selection_or_canonical_or_longest'))
            count_excluded += 1
            continue

        # Add to the list of excluded genes if ENST not found in Ensembl database
        if enst not in ensembl_db:
            out_excl.write('{}\t{}\t{}\n'.format(hgnc_id, symbol, 'not_in_ensembl_db'))
            count_excluded += 1
            continue

        try:

            # Retrieve data about ENST
            transcript = ensembl_db[enst]

            if selected:

                # Calculating CART ID
                if options.prev_cava_db is None:
                    cartidx += 1
                    cart_id = '{}{}'.format(options.series, cartidx)
                else:
                    content = (
                        transcript.strand,
                        len(transcript.exons),
                        helper.read_mrna_sequence(transcript, ref)
                    )

                    if hgnc_id in prev_cava_db and content == prev_cava_db[hgnc_id]['content']:
                        cart_id = '{}{}'.format(options.series, prev_cava_db[hgnc_id]['cartidx'])
                    else:
                        cartidx += 1
                        cart_id = '{}{}'.format(options.series, cartidx)

                template_id = cart_id

            else:
                template_id = enst

            # Add template ID and HGNC ID to transcript
            transcript.id = template_id
            transcript.hgnc_id = hgnc_id

            transcript.assoc_nm = assoc_nm
            transcript.assoc_enst = enst

            # Write IDs to file
            helper.output_ids(out_id, hgnc_id, template_id)

            # Add transcript to database writer
            tdb_writer.add(transcript)

            # Create content of gff3 file
            gff2_lines = helper.create_gff2_lines(transcript, gff2_lines)
            gff3_lines = helper.create_gff3_lines(transcript, gff3_lines)

            # Write to gp file
            helper.output_genepred(transcript, out_genepred)

            # Write to gbk output
            if options.gbk:
                helper.output_gbk(transcript, ref, gbk_dir)

            # Write to fasta file
            helper.output_fasta(transcript, out_fasta, ref)

            # Write annovar files
            if options.annovar:
                helper.output_genepred(transcript, out_genepred_annovar)
                helper.output_fasta_annovar(transcript, out_fasta_annovar, ref)

            if selected:
                count_selected += 1
            else:
                count_canonical_or_longest += 1

        except:
            out_excl.write('{}\t{}\t{}\n'.format(hgnc_id, symbol, 'output_error'))
            count_excluded += 1

    # Create bgzipped, Tabix-index GFF2 and GFF3 outputs
    helper.output_gff2(gff2_lines, options.output + '.gff2')
    helper.output_gff3(gff3_lines, options.output + '.gff3')

    # Finalize outputs
    helper.finalize_outputs(
        options,
        tdb_writer,
        out_fasta,
        out_genepred,
        out_genepred_annovar,
        out_fasta_annovar,
        gbk_dir,
        out_id,
        out_excl
    )

    # Print out summary info
    helper.print_summary_info(options, count_selected, count_canonical_or_longest, count_excluded)


