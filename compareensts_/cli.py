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
        help="..."
    )

    parser.add_option(
        '--enstsx',
        default=None,
        dest='enstsx',
        action='store',
        help="..."
    )

    parser.add_option(
        '--enstsy',
        default=None,
        dest='enstsy',
        action='store',
        help="..."
    )

    parser.add_option(
        '--datax',
        default=None,
        dest='dataX',
        action='store',
        help="Ensembl database file, release 1"
    )

    parser.add_option(
        '--datay',
        default=None,
        dest='dataY',
        action='store',
        help="Ensembl database file, release 2"
    )

    parser.add_option(
        '--ref37',
        default=None,
        dest='ref37',
        action='store',
        help="..."
    )

    parser.add_option(
        '--ref38',
        default=None,
        dest='ref38',
        action='store',
        help="..."
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
        help="Show only discrepant transcripts"
    )

    parser.add_option(
        '--simple',
        default=False,
        dest='simple',
        action='store_true',
        help="Simple output"
    )

    (options, args) = parser.parse_args()

    welcome('CompareENSTs')

    toplevel.run(options)

    goodbye('CompareENSTs')

