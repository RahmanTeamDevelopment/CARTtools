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
        description='FormatCARTs v{}'.format(__version__),
        usage='CARTtools/formatcarts <options>',
        version=__version__
    )

    parser.add_option(
        '--selected_nms',
        dest='selected_nms',
        action='store',
        help="Final selected NMs file (output of SelectNMs)"
    )

    parser.add_option(
        '--selected_ensts',
        dest='selected_ensts',
        action='store',
        help="Selected ENSTs file (output of SelectENSTs)"
    )

    parser.add_option(
        '--canonical',
        dest='canonical',
        action='store',
        help="Canonical ENSTs file"
    )

    parser.add_option(
        '--ensembl',
        default=None,
        dest='ensembl',
        action='store',
        help="Ensembl transcript database (output of EnsemblDB)"
    )

    parser.add_option(
        '--series',
        default=None,
        dest='series',
        action='store',
        help="CARTs series code (e.g. CART37A)"
    )

    parser.add_option(
        '--output',
        default='output',
        dest='output',
        action='store',
        help="Output file name prefix"
    )

    parser.add_option(
        '--gbk',
        default=False,
        dest='gbk',
        action='store_true',
        help="Create GBK output"
    )

    parser.add_option(
        '--ref',
        default=None,
        dest='ref',
        action='store',
        help="Reference genome file"
    )

    parser.add_option(
        '--annovar',
        default=False,
        dest='annovar',
        action='store_true',
        help="Create GenePred and FASTA files for Annovar"
    )

    parser.add_option(
        '--prev_cava_db',
        default=None,
        dest='prev_cava_db',
        action='store',
        help="CAVA db output of previous run from which CART numbering will be continued"
    )

    parser.add_option(
        '--prev_ref',
        default=None,
        dest='prev_ref',
        action='store',
        help="Reference genome of previous run from which CART numbering will be continued"
    )

    (options, args) = parser.parse_args()

    welcome('FormatCARTs')

    toplevel.run(options)

    goodbye('FormatCARTs')