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
        description='RefSeqScan v{}'.format(__version__),
        usage='CARTtools/refseqscan <options>',
        version=__version__
    )

    parser.add_option(
        '--input',
        default=None,
        dest='input',
        action='store',
        help='Input RefSeq database file'
    )

    parser.add_option(
        '--reference',
        default=None,
        dest='reference',
        action='store',
        help='Reference genome file'
    )

    parser.add_option(
        '--output',
        default='output.txt',
        dest='output',
        action='store',
        help='Output file name [default value: %default]'
    )

    (options, args) = parser.parse_args()

    welcome('RefSeqScan')

    toplevel.run(options)

    goodbye('RefSeqScan')
