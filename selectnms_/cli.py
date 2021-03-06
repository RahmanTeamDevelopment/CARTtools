#!env/bin/python

from optparse import OptionParser
import toplevel
import datetime
from main.version import __version__



def welcome(name):

    print '\n' + '=' * 100
    now = str(datetime.datetime.now())
    now = now[:now.find('.')]
    print '{} started: {}\n'.format(name, now)


def goodbye(name):

    now = str(datetime.datetime.now())
    now = now[:now.find('.')]
    print '\n{} finished: {}'.format(name, now)
    print '=' * 100 + '\n'



def start_cli():

    parser = OptionParser(
        description='SelectNMs v{}'.format(__version__),
        usage='CARTtools/selectnms <options>',
        version=__version__
    )

    parser.add_option(
        '--input_genes',
        default=None,
        dest='hgnc',
        action='store',
        help="Input file containing HGNC IDs"
    )

    parser.add_option(
        '--appris',
        default=None,
        dest='appr',
        action='store',
        help="APPRIS data file"
    )

    parser.add_option(
        '--refsdb',
        default=None,
        dest='refsdb',
        action='store',
        help="RefSeq transcript database file (output of RefSeqDB)"
    )

    parser.add_option(
        '--refschk',
        default=None,
        dest='refss',
        action='store',
        help="RefSeqCheck output file"
    )

    parser.add_option(
        '--genes_dict',
        default=None,
        dest='genes',
        action='store',
        help="Gene ID dictionary file"
    )

    parser.add_option(
        '--build',
        default=None,
        dest='build',
        action='store',
        help="Genome build (GRCh37 or GRCh38)"
    )

    parser.add_option(
        '--out_auto',
        dest='out',
        action='store',
        help="Output file name prefix for automatic selection"
    )

    parser.add_option(
        '--out',
        dest='out_final',
        action='store',
        help="Output file name prefix for final selection"
    )

    parser.add_option(
        '--refsdbinc',
        dest='refsdbinc',
        action='store',
        help="List of transcripts included in the RefSeq DB (output of RefSeqDB)"
    )


    (options, args) = parser.parse_args()

    welcome('SelectNMs')

    toplevel.run(options)

    goodbye('SelectNMs')
