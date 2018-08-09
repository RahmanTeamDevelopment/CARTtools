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
        description='SelectENSTs v{}'.format(__version__),
        usage='CARTtools/selectensts <options>',
        version=__version__
    )

    parser.add_option(
        '--mapped_nms',
        default=None,
        dest='mapped_nms',
        action='store',
        help="Mapped NMs data file"
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
        help="Match a mapped NM to ENSTs with any gene symbol"
    )

    parser.add_option(
        '--input',
        default=None,
        dest='input',
        action='store',
        help="Input file (list of NMs)"
    )

    parser.add_option(
        '--output',
        default=None,
        dest='output',
        action='store',
        help="Output file name prefix"
    )

    (options, args) = parser.parse_args()

    welcome('SelectENSTs')

    toplevel.run(options)

    goodbye('SelectENSTs')