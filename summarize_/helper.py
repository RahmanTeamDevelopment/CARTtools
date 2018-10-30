from datetime import datetime
import os


def create_summary_table(options):

    full = os.path.dirname(os.path.realpath(__file__))
    rootdir = full[:full.rfind('/env')]

    nms_data, header_original = read_nms('{}_nms_GRCh37_final_selected.txt'.format(options.prefix))
    ensts_37 = read_ensts('{}_selected_ensts_GRCh37.txt'.format(options.prefix))
    ensts_38 = read_ensts('{}_selected_ensts_GRCh38.txt'.format(options.prefix))
    compared_ensts = read_compared_ensts('{}_compared_ensts.txt'.format(options.prefix))
    ids_37 = read_ids('{}_GRCh37_ids.txt'.format(options.output))
    ids_38 = read_ids('{}_GRCh38_ids.txt'.format(options.output))

    out = open('{}_summary.txt'.format(options.output), 'w')

    header = header_original.split('\t') + [
        'CART_ASSOCIATED_ENST_37',
        'UTR_DIFF_37',
        'UTR_EXON_NUM_DIFF_37',
        'CART_ASSOCIATED_ENST_38',
        'UTR_DIFF_38',
        'UTR_EXON_NUM_DIFF_38',
        'CART37_and_CART38_DIFF',
        'CART_ID_37',
        'CART_ID_38'
    ]
    out.write('\t'.join(header) + '\n')

    for k, v in nms_data.iteritems():

        hgnc_id = nms_data[k][0]

        if hgnc_id in ensts_37:
            enst37 = ensts_37[hgnc_id]
        else:
            enst37 = None

        if hgnc_id in ensts_38:
            enst38 = ensts_38[hgnc_id]
        else:
            enst38 = None

        if enst37 is not None:
            columns37 = [enst37['ENST'], enst37['UTR_DIFF'], enst37['UTR_EXON_NUM_DIFF']]
        else:
            columns37 = ['No_ENST'] * 3

        if enst38 is not None:
            columns38 = [enst38['ENST'], enst38['UTR_DIFF'], enst38['UTR_EXON_NUM_DIFF']]
        else:
            columns38 = ['No_ENST'] * 3

        if hgnc_id in compared_ensts:
            enst_diff = compared_ensts[hgnc_id]
        else:
            enst_diff = 'No_ENST_pair'

        if hgnc_id in ids_37:
            cart_id_37 = ids_37[hgnc_id]
        else:
            cart_id_37 = '.'

        if hgnc_id in ids_38:
            cart_id_38 = ids_38[hgnc_id]
        else:
            cart_id_38 = '.'

        out.write('\t'.join([k] + nms_data[k] + columns37 + columns38 + [enst_diff] + [cart_id_37, cart_id_38]) + '\n')

    out.close()

    print '- Summary table of CARTs created: {}'.format('{}_summary.txt'.format(options.output))


def create_report(options):

    start = datetime.strptime(options.start, '%Y-%m-%d_%H:%M:%S')
    end = datetime.strptime(options.end, '%Y-%m-%d_%H:%M:%S')
    runtime = str(end - start)
    cwd = os.getcwd()

    out = open('{}_report.txt'.format(options.output), 'w')

    title = 'CART Pipeline Report (created: {})'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    out.write('\n{}\n'.format(title))
    out.write('='*len(title)+'\n')

    out_section(out, '1 - General info')
    out.write('Start time: {}\n'.format(start))
    out.write('End time: {}\n'.format(end))
    out.write('Total runtime: {}\n'.format(runtime))
    out.write('Config file: {}\n'.format(options.configfn))
    out.write('Input file name (gene list): {} ({} genes)\n'.format(options.inputfn, count_lines(options.inputfn)))


    out_section(out, '2 - Mapping of all NMs in interim RefSeq release')
    out.write('NCBI mapping to GRCh37: {} NMs ({} NMs excluded)\n'.format(
            count_lines('{}_refseqdb_GRCh37_ncbi_included.txt'.format(options.prefix)),
            count_lines('{}_refseqdb_GRCh37_ncbi_excluded.txt'.format(options.prefix))
        )
    )
    out.write('UCSC mapping to GRCh37: {} NMs ({} NMs excluded)\n'.format(
        count_lines('{}_refseqdb_GRCh37_ucsc_included.txt'.format(options.prefix)),
        count_lines('{}_refseqdb_GRCh37_ucsc_excluded.txt'.format(options.prefix))
        )
    )
    out.write('NCBI mapping to GRCh38: {} NMs ({} NMs excluded)\n'.format(
        count_lines('{}_refseqdb_GRCh38_ncbi_included.txt'.format(options.prefix)),
        count_lines('{}_refseqdb_GRCh38_ncbi_excluded.txt'.format(options.prefix))
        )
    )
    out.write('UCSC mapping to GRCh38: {} NMs ({} NMs excluded)\n'.format(
        count_lines('{}_refseqdb_GRCh38_ucsc_included.txt'.format(options.prefix)),
        count_lines('{}_refseqdb_GRCh38_ucsc_excluded.txt'.format(options.prefix))
        )
    )

    out_section(out, '3 - Selected NMs')

    out.write('Automatical selection: {} genes ({} missing genes)\n'.format(
            count_lines('{}_nms_GRCh37_auto_selected.txt'.format(options.prefix)),
            count_lines('{}_nms_GRCh37_auto_missing.txt'.format(options.prefix))
        )
    )
    selected, missing = count_genes('{}_nms_GRCh37_final_selected.txt'.format(options.prefix))
    out.write('Final selection: {} genes ({} missing genes)\n'.format(selected, missing))

    out_section(out, '4 - Mapping of selected NMs')
    ncbi, ucsc = count_mapped_nms('{}_mapped_nms_GRCh37_source.txt'.format(options.prefix))
    missing = count_lines('{}_mapped_nms_GRCh37_missing.txt'.format(options.prefix))
    out.write('Mapped to GRCh37: {} genes (NCBI: {}, UCSC: {}) - {} genes missing\n'.format(
            ncbi+ucsc,
            ncbi,
            ucsc,
            missing
        )
    )

    ncbi, ucsc = count_mapped_nms('{}_mapped_nms_GRCh38_source.txt'.format(options.prefix))
    missing = count_lines('{}_mapped_nms_GRCh38_missing.txt'.format(options.prefix))
    out.write('Mapped to GRCh38: {} genes (NCBI: {}, UCSC: {}) - {} genes missing\n'.format(
            ncbi + ucsc,
            ncbi,
            ucsc,
            missing
        )
    )

    out_section(out, '5 - Ensembl transcript databases')

    out.write('Release {} (GRCh37): {} ENSTs\n'.format(
            options.ens37,
            count_lines('{}_ensembl_{}.txt'.format(options.prefix, options.ens37))
        )
    )
    out.write('Release {} (GRCh38): {} ENSTs\n'.format(
            options.ens38,
            count_lines('{}_ensembl_{}.txt'.format(options.prefix, options.ens38))
        )
    )

    out_section(out, '6 - Selected ENSTs')
    enst, noenst = count_selected_ensts('{}_selected_ensts_GRCh37.txt'.format(options.prefix))
    out.write('GRCh37: ENST selected for {} genes (no ENST for {} genes)\n'.format(
            enst,
            noenst
        )
    )
    enst, noenst = count_selected_ensts('{}_selected_ensts_GRCh38.txt'.format(options.prefix))
    out.write('GRCh38: ENST selected for {} genes (no ENST for {} genes)\n'.format(
            enst,
            noenst
        )
    )

    out_section(out, '7 - ENST comparison across builds')
    count_all, count_identical, count_cds_identical = count_compared_ensts('{}_compared_ensts.txt'.format(options.prefix))
    out.write('Genes with ENST for both builds: {}\n'.format(count_all))
    out.write('Identical ENSTs: {}\n'.format(count_identical))
    out.write('ENSTs with identical CDS: {}\n'.format(count_cds_identical))

    out_section(out, '8 - Genes outputted')
    count_all_37, count_selected_37, count_canonical_37 = count_outputted_genes('{}_GRCh37_ids.txt'.format(options.output))
    count_excluded_37 = count_excluded_genes('{}_GRCh37_excl.txt'.format(options.output))
    count_all_38, count_selected_38, count_canonical_38 = count_outputted_genes('{}_GRCh38_ids.txt'.format(options.output))
    count_excluded_38 = count_excluded_genes('{}_GRCh38_excl.txt'.format(options.output))
    out.write('GRCh37: {} genes outputted (selected CART: {}, canonical/longest ENST: {}) - {} genes excluded\n'.format(
        count_all_37,
        count_selected_37,
        count_canonical_37,
        count_excluded_37
    ))
    out.write('GRCh38: {} genes outputted (selected CART: {}, canonical/longest ENST: {}) - {} genes excluded\n'.format(
        count_all_38,
        count_selected_38,
        count_canonical_38,
        count_excluded_38
    ))

    out_section(out, '9 - Outputs')
    out.write('Output files: {}/\n'.format(cwd))
    out.write('Intermediate files: {}/cart_pipeline_files/\n'.format(cwd))


    out.write('\n')
    out.close()

    print '- Final report created: {}'.format('{}_report.txt'.format(options.output))


# Helper functions:

def read_nms(fn):

    ret = {}
    header = ''
    for line in open(fn):
        line = line.strip()
        if line[0] == '#':
            header = line
            continue
        if line == '':
            continue
        cols = line.split('\t')
        ret[cols[0]] = cols[1:]
    return ret, header


def read_ensts(fn):

    ret = {}
    for line in open(fn):
        line = line.strip()
        if line[0] == '#' or line == '':
            continue
        cols = line.split('\t')
        ret['HGNC:' + cols[2]] = {
            'ENST': cols[3] if cols[3] != '.' else cols[8],
            'UTR_DIFF': cols[6],
            'UTR_EXON_NUM_DIFF': cols[7]
        }

    return ret


def read_ids(fn):

    ret = {}
    for line in open(fn):
        line = line.strip()
        if line[0] == '#' or line == '':
            continue
        cols = line.split('\t')
        ret[cols[0]] = cols[1]
    return ret


def out_section(out, s):

    out.write('\n{}\n'.format(s))
    out.write('-'*len(s)+'\n')


def count_lines(fn):

    ret = 0
    for line in open(fn):
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        ret += 1
    return ret


def count_genes(fn):

    selected = 0
    missing = 0
    for line in open(fn):
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        cols = line.split('\t')
        if cols[14] == '.':
            missing += 1
        else:
            selected += 1
    return selected, missing


def count_mapped_nms(fn):

    ncbi = 0
    ucsc = 0
    for line in open(fn):
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        cols = line.split('\t')
        if cols[1] == 'NCBI':
            ncbi += 1
        else:
            ucsc += 1
    return ncbi, ucsc


def count_selected_ensts(fn):

    enst = 0
    noenst = 0
    for line in open(fn):
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        cols = line.split('\t')
        if cols[3].startswith('ENST'):
            enst += 1
        else:
            noenst += 1
    return enst, noenst


def count_compared_ensts(fn):

    count_all = 0
    count_identical = 0
    count_cds_identical = 0
    for line in open(fn):
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        cols = line.split('\t')
        count_all += 1

        if cols[3] == '.':
            count_identical += 1

        if 'NFX' in cols[3] or 'NFY' in cols[3]:
            continue

        if 'CDS' not in cols[3]:
            count_cds_identical += 1

    return count_all, count_identical, count_cds_identical


def count_outputted_genes(fn):

    count_all = 0
    count_selected = 0
    count_canonical = 0
    for line in open(fn):
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        cols = line.split('\t')
        count_all += 1
        if 'CART' in cols[1]:
            count_selected += 1
        else:
            count_canonical += 1
    return count_all, count_selected, count_canonical


def count_excluded_genes(fn):

    count = 0
    for line in open(fn):
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        count += 1
    return count


def read_compared_ensts(fn):

    ret = {}
    for line in open(fn):
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        cols = line.split('\t')
        ret[cols[0]] = cols[3]
    return ret


def read_gdm_categories(fn):

    ret = {}
    for line in open(fn):
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        cols = line.split('\t')
        ret[cols[3]] = cols[2]
    return ret
