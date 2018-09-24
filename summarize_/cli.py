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
        description='Summarize v{}'.format(__version__),
        usage='CARTtools/summarize <options>',
        version=__version__
    )

    parser.add_option(
        '--prefix',
        default=None,
        dest='prefix',
        action='store',
        help="Input file prefix"
    )

    parser.add_option(
        '--start',
        default=None,
        dest='start',
        action='store',
        help="Start date/time"
    )

    parser.add_option(
        '--end',
        default=None,
        dest='end',
        action='store',
        help="End date/time"
    )

    parser.add_option(
        '--inputfn',
        default=None,
        dest='inputfn',
        action='store',
        help="Gene input file name used in the pipeline"
    )

    parser.add_option(
        '--configfn',
        default=None,
        dest='configfn',
        action='store',
        help="Config file name used in the pipeline"
    )

    parser.add_option(
        '--ens37',
        default=None,
        dest='ens37',
        action='store',
        help="Ensembl transcript database (GRCh37)"
    )

    parser.add_option(
        '--ens38',
        default=None,
        dest='ens38',
        action='store',
        help="Ensembl transcript database (GRCh38)"
    )

    parser.add_option(
        '--output',
        default=None,
        dest='output',
        action='store',
        help="Output file name prefix"
    )




    (options, args) = parser.parse_args()

    welcome('Summarize')

    toplevel.run(options)

    goodbye('Summarize')

