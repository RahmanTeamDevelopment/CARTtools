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
        description='EnsemblDB v{}'.format(__version__),
        usage='CARTtools/ensembldb <options>',
        version=__version__
    )

    parser.add_option(
        '--input',
        default=None,
        dest='input',
        action='store',
        help="Input filename (list of ENST IDs)"
    )

    parser.add_option(
        '--output',
        default=None,
        dest='output',
        action='store',
        help="Output filename prefix"
    )

    parser.add_option(
        '--release',
        default=None,
        dest='ensembl',
        action='store',
        help="Ensembl release version"
    )

    (options, args) = parser.parse_args()

    welcome('EnsemblDB')

    toplevel.run(options)

    goodbye('EnsemblDB')