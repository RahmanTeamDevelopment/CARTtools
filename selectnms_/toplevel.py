from __future__ import division
import sys
import helper
import automatic
import community
from main.version import __version__


def count_records(fn):

    count = 0
    for line in open(fn):
        if not line[0] == '#':
            count += 1
    return count


def run(options):

    count_input = count_records(options.hgnc)
    print 'Input file: {} ({} genes)'.format(options.hgnc, count_input)

    automatic.run(options)

    community.run(options)

    count_included = count_records(options.out + '_selected.txt')
    count_missing = count_records(options.out + '_missing.txt')

    print '\nSummary:'
    print ' - Genes included in the output: {} (-> {})'.format(count_included, options.out + '_selected.txt')
    print ' - Missing genes: {} (-> {})'.format(count_missing, options.out + '_missing.txt')