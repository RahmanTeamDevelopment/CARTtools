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
        description='CompareENSTs v{}'.format(__version__),
        usage='CARTtools/compareensts <options>',
        version=__version__
    )

    parser.add_option(
        '--input',
        default=None,
        dest='input',
        action='store',
        help="Input file"
    )

    parser.add_option(
        '--ensts37',
        default=None,
        dest='enstsx',
        action='store',
        help="Table of selected ENSTs for GRCh37 (output of SelectedENSTs)"
    )

    parser.add_option(
        '--ensts38',
        default=None,
        dest='enstsy',
        action='store',
        help="Table of selected ENSTs for GRCh38 (output of SelectedENSTs)"
    )

    parser.add_option(
        '--data37',
        default=None,
        dest='dataX',
        action='store',
        help="Ensembl database for GRCh37 (output of EnsemblDB)"
    )

    parser.add_option(
        '--data38',
        default=None,
        dest='dataY',
        action='store',
        help="Ensembl database for GRCh38 (output of EnsemblDB)"
    )

    parser.add_option(
        '--ref37',
        default=None,
        dest='ref37',
        action='store',
        help="Reference genome fasta file (GRCh37)"
    )

    parser.add_option(
        '--ref38',
        default=None,
        dest='ref38',
        action='store',
        help="Reference genome fasta file (GRCh38)"
    )

    parser.add_option(
        '--output',
        default=None,
        dest='output',
        action='store',
        help="Output file name"
    )

    parser.add_option(
        '--discr',
        default=False,
        dest='discr',
        action='store_true',
        help="Output only ENST pairs where there are differences"
    )

    parser.add_option(
        '--simple',
        default=False,
        dest='simple',
        action='store_true',
        help="Create simplified output"
    )

    (options, args) = parser.parse_args()

    welcome('CompareENSTs')

    toplevel.run(options)

    goodbye('CompareENSTs')

