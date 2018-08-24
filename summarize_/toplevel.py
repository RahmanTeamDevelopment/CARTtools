
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


def run(options):

    nms_data, header_original = read_nms(options.nms)
    ensts_37 = read_ensts(options.ensts37)
    ensts_38 = read_ensts(options.ensts38)

    out = open(options.output, 'w')

    header = header_original.split('\t') + [
        'Selected_ENST_37',
        'UTR_DIFF_37',
        'UTR_EXON_NUM_DIFF_37',
        'Selected_ENST_38',
        'UTR_DIFF_38',
        'UTR_EXON_NUM_DIFF_38'
    ]
    out.write('\t'.join(header)+'\n')

    for k, v in nms_data.iteritems():

        if k in ensts_37:
            enst37 = ensts_37[k]
        else:
            enst37 = None

        if k in ensts_38:
            enst38 = ensts_38[k]
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

        out.write('\t'.join([k] + nms_data[k] + columns37 + columns38) + '\n')

    out.close()






