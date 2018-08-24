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


def count_final_records(fn):

    count = count_missing = 0
    for line in open(fn):
        if not line[0] == '#':
            count += 1
            line = line.strip()
            cols = line.split('\t')
            if cols[-2] == '.':
                count_missing += 1
    return count, count_missing


def run(options):

    count_input = count_records(options.hgnc)
    print 'Input file: {} ({} genes)'.format(options.hgnc, count_input)

    sys.stdout.write('\nAutomatic selection process running ... ')
    sys.stdout.flush()

    automatic.run(options)

    print 'done'

    count_included = count_records(options.out + '_auto_selected.txt')
    count_missing = count_records(options.out + '_auto_missing.txt')

    print '\nAutomatic selection summary:'
    print ' - Genes included in the output: {} (-> {})'.format(count_included, options.out + '_auto_selected.txt')
    print ' - Missing genes: {} (-> {})'.format(count_missing, options.out + '_auto_missing.txt')

    sys.stdout.write('\nFinal selection process running ... ')
    sys.stdout.flush()

    community.run(options, options.out + '_auto_selected.txt', options.out + '_auto_missing.txt')

    print 'done'

    count, count_missing = count_final_records(options.out_final + '.txt')

    print '\nFinal selection summary:'
    print ' - Genes with NM selection: {}'.format(count - count_missing)
    print ' - Genes with no NM selection: {}'.format(count_missing)

    print '\nFinal selection output file: {}'.format(options.out_final + '.txt')