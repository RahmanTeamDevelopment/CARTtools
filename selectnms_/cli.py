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
        help="APPRIS file"
    )

    parser.add_option(
        '--refsdb',
        default=None,
        dest='refsdb',
        action='store',
        help="RefSeq transcript database file"
    )

    parser.add_option(
        '--refss',
        default=None,
        dest='refss',
        action='store',
        help="refseq_scan output file"
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
        help="Genome build"
    )

    parser.add_option(
        '--series',
        dest='series',
        default=None,
        action='store',
        help="CART series"
    )

    parser.add_option(
        '--out',
        dest='out',
        action='store',
        help="Output file name"
    )

    (options, args) = parser.parse_args()

    welcome('SelectNMs')

    toplevel.run(options)

    goodbye('SelectNMs')
