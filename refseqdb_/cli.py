#!env/bin/python

from optparse import OptionParser
import toplevel
from main.version import __version__


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
        help="Genome build [default value: %default]"
    )

    parser.add_option(
        '--mapping',
        default='ncbi',
        dest='mapping',
        action='store',
        help='Mapping to use (ncbi or ucsc) [default value: %default]'
    )

    parser.add_option(
        '--output',
        default='output',
        dest='output',
        action='store',
        help='Output file name prefix [default value: %default]'
    )

    (options, args) = parser.parse_args()
    toplevel.run(options)
