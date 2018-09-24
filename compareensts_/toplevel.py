import os
import sys
import helper
import reference


def run(options):


    # Checks if Ensembl TXT files exist if required
    if not os.path.isfile(options.dataX[:-3] + '.txt'):
        print 'Error: Dataset X txt file (' + options.dataX[:-3] + '.txt) cannot be found.\n'
        quit()
    if not os.path.isfile(options.dataY[:-3] + '.txt'):
        print 'Error: Dataset Y txt file (' + options.dataY[:-3] + '.txt) cannot be found.\n'
        quit()

    input_genes = helper.read_input_genes(options.input)

    # Read transcript database X
    sys.stdout.write('GRCh37 Ensembl db [' + options.dataX + ']: READING...')
    sys.stdout.flush()
    dataX, N = helper.readEnsemblData(options.dataX)

    sys.stdout.write('\rGRCh37 Ensembl db [' + options.dataX + ']: ' + str(N) + ' transcripts')
    print ''

    # Read transcript database Y
    sys.stdout.write('GRCh38 Ensembl db [' + options.dataY + ']: READING...')
    sys.stdout.flush()
    dataY, N = helper.readEnsemblData(options.dataY)

    sys.stdout.write('\rGRCh38 Ensembl db [' + options.dataY + ']: ' + str(N) + ' transcripts')
    print '\n'


    # Initialize GRCh37 reference genome
    if not os.path.isfile(options.ref37):
        print '\nError: GRCh37 reference genome file (' + options.ref37 + ') cannot be found.\n'
        quit()
    ref_GRCh37 = reference.Reference(options.ref37)


    # Initialize GRCh38 reference genome
    if not os.path.isfile(options.ref38):
        print '\nError: GRCh38 reference genome file (' + options.ref38 + ') cannot be found.\n'
        quit()
    ref_GRCh38 = reference.Reference(options.ref38)

    ensts_37 = helper.read_enst_file(options.enstsx)
    ensts_38 = helper.read_enst_file(options.enstsy)

    # Initialize output file
    out = open(options.output, 'w')
    out.write('\t'.join(['#GENE', 'ENST_37', 'ENST_38', 'DIFFERENCE']) + '\n')


    # Iterate through the input list of transcripts

    count_gene_in_both = 0

    n_identical = 0
    n_cds_identical = 0
    i = 1
    for g in input_genes:

        i += 1

        if g not in ensts_37 or g not in ensts_38:
            continue

        count_gene_in_both += 1

        enst1 = ensts_37[g]
        enst2 = ensts_38[g]

        flags = []
        comparewith = []

        sys.stdout.write('\rAnalysing gene ' + str(i) + '/' + str(len(input_genes)))
        sys.stdout.flush()

        if enst1 not in dataX.keys():
            flags.append('NFX')
        else:
            transcript = dataX[enst1]

        if enst2 not in dataY.keys():
            flags.append('NFY')
        else:
            comparewith = [dataY[enst2]]

        if len(flags) > 0:
            out.write('\t'.join(['HGNC:' + g, enst1, enst2, ';'.join(flags)]) + '\n')
            continue

        identical, cds_identical = helper.compare(
            g, transcript, comparewith, ref_GRCh37, ref_GRCh38, options, out)
        if identical:
            n_identical += 1
        if cds_identical:
            n_cds_identical += 1

    print ' - Done.'

    # Close output file
    out.close()

    # Goodbye message
    print '\nSummary:'
    print '- {} of the {} genes are on both ENSTs lists'.format(count_gene_in_both, len(input_genes))
    print '- ' + str(n_identical) + ' genes have identical ENSTs'
    print '- ' + str(n_cds_identical) + ' genes have CDS-identical ENSTs'
    print '\nOutput written to file: ' + options.output


