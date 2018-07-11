from optparse import OptionParser
import toplevel
from main.version import __version__


def start_cli():

    parser = OptionParser(
        description='RefSeqScan v{}'.format(__version__),
        usage='CARTtools/refseqscan <options>',
        version=__version__
    )

    parser.add_option(
        '--input',
        default=None,
        dest='input',
        action='store',
        help='Input RefSeq database file'
    )

    parser.add_option(
        '--reference',
        default=None,
        dest='reference',
        action='store',
        help='Reference genome file'
    )

    parser.add_option(
        '--output',
        default='output.txt',
        dest='output',
        action='store',
        help='Output file name [default value: %default]'
    )

    (options, args) = parser.parse_args()
    toplevel.run(options)
