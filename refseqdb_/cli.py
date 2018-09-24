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
        description='RefSeqDB v{}'.format(__version__),
        usage='CARTtools/refseqdb <options>',
        version=__version__
    )

    parser.add_option(
        '--build',
        default='GRCh37',
        dest='build',
        action='store',
        help="Genome build (GRCh37 or GRCh38)"
    )

    parser.add_option(
        '--mapping',
        default='ncbi',
        dest='mapping',
        action='store',
        help='Mapping source (ncbi or ucsc)'
    )

    parser.add_option(
        '--output',
        default='output',
        dest='output',
        action='store',
        help='Output file name prefix'
    )

    (options, args) = parser.parse_args()

    welcome('RefSeqDB')

    toplevel.run(options)

    goodbye('RefSeqDB')