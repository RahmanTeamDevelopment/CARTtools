#!env/bin/python

from optparse import OptionParser
import toplevel
from main.version import __version__
import datetime


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
        description='MapNMs v{}'.format(__version__),
        usage='CARTtools/mapnms <options>',
        version=__version__
    )

    parser.add_option(
        '--input',
        default=None,
        dest='input',
        action='store',
        help="Input file (output of SelectNMs)"
    )

    parser.add_option(
        '--ncbi',
        default=None,
        dest='ncbi',
        action='store',
        help="RefSeqDB output file with NCBI interim mapping data"
    )

    parser.add_option(
        '--ucsc',
        default=None,
        dest='ucsc',
        action='store',
        help="RefSeqDB output file with UCSC mapping data"
    )

    parser.add_option(
        '--hgncid_to_symbol',
        default=None,
        dest='hgnc',
        action='store',
        help="HGNC ID to Gene Symbol dictionary file"
    )

    parser.add_option(
        '--output',
        default='output',
        dest='output',
        action='store',
        help="Output file name prefix"
    )

    parser.add_option(
        '--more_symbols',
        default=None,
        dest='symbols',
        action='store',
        help="TXT file for specifying gene symbols for missing HGNC IDs"
    )

    (options, args) = parser.parse_args()

    welcome('MapNMs')

    toplevel.run(options)

    goodbye('MapNMs')

