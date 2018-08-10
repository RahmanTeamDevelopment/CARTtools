import subprocess
import os


def parse_config_file(fn):

    ret = {}
    section = ''
    with open(fn) as infile:
        for line in infile:
            line = line.strip()
            if line == '' or line.startswith('#'):
                continue

            if line[0] == '[' and line[-1] == ']':
                s = line[1:-1]
                if '[' in s or ']' in s or '=' in s:
                    continue
                section = s

            elif line.count('=') == 1:
                key, value = line.split('=')
                key = key.strip()
                value = value.strip()
                if section != '':
                    key = '{}.{}'.format(section, key)
                ret[key] = value
    return ret


def map_keys(config):

    ret = {}

    map_of_keys = {
        'reference_37': 'ref37',
        'reference_38': 'ref38',
        'select_nm.input_genes': 'inputgenes',
        'select_nm.appris': 'apprisfile',
        'select_nm.genes_dict': 'genesdict',
        'ensembl_db.release_37': 'ens37',
        'ensembl_db.release_38': 'ens38',
        'map_nm.hgncid_to_symbol': 'hgncidtosymbol',
        'select_enst.gene_synonyms': 'genesynonyms',
        'format_cart.series': 'series',
        'format_cart.output': 'output'
    }

    for k in map_of_keys.keys():
        if k in config:
            ret[map_of_keys[k]] = config[k]

    return ret


def add_default_values(c, rootdir):

    default_values = {
        'apprisfile': '{}/default/appris_refseq107.txt'.format(rootdir),
        'genesdict': '{}/default/hgnc_biomart-05052017.txt'.format(rootdir),
        'hgncidtosymbol': '{}/default/20180104_HGNCID_to_GeneSymbol_from_HGNC_BioMart.txt'.format(rootdir),
        'ens37': '75',
        'ens38': '92',
        'genesynonyms': '{}/default/gene_synonyms.txt'.format(rootdir),
        'series': 'CART37A',
        'output': 'CART37A'
    }

    for k, v in default_values.iteritems():
        if k not in c or c[k] == '.':
            c[k] = v

    return c


def run(options):

    rootdir = os.path.dirname(os.path.realpath(__file__)).split('/env')[0]

    config = parse_config_file(options.config)
    mapped_args = map_keys(config)
    mapped_args = add_default_values(mapped_args, rootdir)

    args_to_pass = []
    for k, v in mapped_args.iteritems():
        args_to_pass.append('{}:{}'.format(k, v))

    subprocess.call(['{}/pipeline/pipeline.sh'.format(rootdir)] + args_to_pass)

