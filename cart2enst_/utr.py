import phase1_utr


def utr_selection(cart, candidates, log):

    filtered, reason, phase1_difference_type, phase1_decisive_criteria = _utr_selection_utr5(cart, candidates, log)
    if len(filtered) == 1:
        if filtered[0] is None:
            print 'return value is None'
        else:
            return filtered[0], reason, phase1_difference_type, phase1_decisive_criteria

    return cart, 'requires_utr3_selection', '.', '.'


def _utr_selection_utr5(cart, candidates, log):

    filtered = [enst for enst in candidates if len(enst.utr5_exons) == len(cart.utr5_exons)]
    if len(filtered) == 1:
        return [filtered[0]], 'single_same_num_of_utr5_exons', '.', '.'
    if len(filtered) > 0:
        candidates = filtered

    filtered = [enst for enst in candidates if _enst_utr5_encompass_cart_utr5(cart, enst)]
    if len(filtered) == 1:
        return [filtered[0]], 'single_encompassing_utr5', '.', '.'
    if len(filtered) > 0:
        candidates = filtered

    sel, difference_type, decisive_criteria = phase1_utr.utr_selection(candidates, log)

    return [sel], 'phase1_utr_selection', difference_type, decisive_criteria


def _enst_utr5_encompass_cart_utr5(cart, enst):

    if cart.strand == '+':
        return enst.transcript_start <= cart.transcript_start
    else:
        return enst.transcript_start >= cart.transcript_start


def utr_difference(cart, enst):

    ret = []
    if cart.utr5_exons != enst.utr5_exons:
        ret.append('UTR5')
    if cart.utr3_exons != enst.utr3_exons:
        ret.append('UTR3')

    if len(ret) == 0:
        return '.'
    else:
        return ','.join(ret)


def utr_exon_number_difference(cart, enst):

    ret = []

    delta_utr5 = len(enst.utr5_exons) - len(cart.utr5_exons)
    delta_utr3 = len(enst.utr3_exons) - len(cart.utr3_exons)

    if delta_utr5 != 0:
        s = 'UTR5:'
        if delta_utr5 > 0:
            s += '+'
        s += str(delta_utr5)
        ret.append(s)

    if delta_utr3 != 0:
        s = 'UTR3:'
        if delta_utr3 > 0:
            s += '+'
        s += str(delta_utr3)
        ret.append(s)

    if len(ret) == 0:
        return '.'
    else:
        return ','.join(ret)
