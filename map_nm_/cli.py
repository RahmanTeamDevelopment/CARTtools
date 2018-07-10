#!env/bin/python

from optparse import OptionParser
import toplevel
from main.version import __version__


def start_cli():

    parser = OptionParser(
        description='map_NM v{}'.format(__version__),
        usage='CARTtools/map_nm <options>',
        version=__version__
    )

    parser.add_option(
        '--input',
        default=None,
        dest='input',
        action='store',
        help="Input file name"
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
        '--hgnc',
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
        help="Output file name prefix [default value: %default]"
    )

    parser.add_option(
        '--symbols',
        default=None,
        dest='symbols',
        action='store',
        help="Txt file for specifying gene symbols for missing HGNC IDs"
    )



    (options, args) = parser.parse_args()
    toplevel.run(options)
