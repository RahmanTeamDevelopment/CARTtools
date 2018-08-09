#!env/bin/python

from optparse import OptionParser
import toplevel
from main.version import __version__


def start_cli():

    parser = OptionParser(
        description='CART_pipeline v{}'.format(__version__),
        usage='CARTtools/cart_pipeline <options>',
        version=__version__
    )

    parser.add_option(
        '--config',
        default=None,
        dest='config',
        action='store',
        help="Configuration file"
    )

    (options, args) = parser.parse_args()
    toplevel.run(options)
