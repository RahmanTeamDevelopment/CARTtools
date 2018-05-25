#!env/bin/python

from optparse import OptionParser
import toplevel
from .version import __version__


def start_cli():

    parser = OptionParser(
        description='CART2ENST v{}'.format(__version__),
        usage='CART2ENST/cart2enst <options>',
        version=__version__
    )

    parser.add_option(
        '--cart_data',
        default=None,
        dest='cart_data',
        action='store',
        help="CART data file"
    )

    parser.add_option(
        '--ensembl_data',
        default=None,
        dest='ensembl_data',
        action='store',
        help="Ensembl data file"
    )

    parser.add_option(
        '--gene_synonyms',
        default=None,
        dest='gene_synonyms',
        action='store',
        help="Gene synonyms file"
    )

    parser.add_option(
        '--any_gene',
        default=False,
        dest='any_gene',
        action='store_true',
        help="Match a CART to ENSTs with any gene symbol"
    )

    parser.add_option(
        '--input',
        default=None,
        dest='input',
        action='store',
        help="Input file (list of CART IDs)"
    )

    parser.add_option(
        '--output',
        default=None,
        dest='output',
        action='store',
        help="Output file name prefix"
    )

    (options, args) = parser.parse_args()
    toplevel.run(options)
