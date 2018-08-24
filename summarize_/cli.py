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
        '--nms',
        default=None,
        dest='nms',
        action='store',
        help="..."
    )

    parser.add_option(
        '--ensts37',
        default=None,
        dest='ensts37',
        action='store',
        help="..."
    )

    parser.add_option(
        '--ensts38',
        default=None,
        dest='ensts38',
        action='store',
        help="..."
    )

    parser.add_option(
        '--output',
        default=None,
        dest='output',
        action='store',
        help="..."
    )




    (options, args) = parser.parse_args()

    welcome('Summarize')

    toplevel.run(options)

    goodbye('Summarize')

