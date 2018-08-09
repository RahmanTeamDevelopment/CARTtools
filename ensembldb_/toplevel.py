from __future__ import division
import os
import helper


def run(options):

    # Checking if all required options are specified
    if options.ensembl is None:
        print '\nError: no Ensembl release specified. Use option -h to get help!\n'
        quit()

    try:
        options.ensembl = int(options.ensembl)
    except:
        print '\nError: Ensembl release specified is not an integer. Use option -h to get help!\n'
        quit()

    if options.output is None:
        print '\nError: no output file name specified. Use option -h to get help!\n'
        quit()

    # Must use Ensembl release >= 70 or v65
    if not (options.ensembl >= 70 or options.ensembl == 65):
        print '\nError: This version works with Ensembl v65 or >= v70.\n'
        quit()

    # Genome build
    genome_build = 'GRCh37' if options.ensembl <= 75 else 'GRCh38'

    # Print info
    print 'Ensembl version:  ' + str(options.ensembl)
    print 'Reference genome: ' + genome_build

    # Creating compressed output file
    Nretrieved = helper.process_data(options, genome_build)
    print '\nA total of ' + str(Nretrieved) + ' transcripts have been included\n'

    # Indexing output file with Tabix
    helper.indexFile(options)

    # Removing uncompressed output file
    os.remove(options.output)

    # Printing out summary information
    print ''
    print '---------------------'
    print 'Output files created:'
    print '---------------------'
    print options.output + '.gz (transcript database)'
    print options.output + '.gz.tbi (index file)'
    print options.output + '.txt (list of transcripts)'

