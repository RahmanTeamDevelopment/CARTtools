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
                section = s.lower()

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
        'selectnms.input_genes': 'inputgenes',
        'selectnms.appris': 'apprisfile',
        'selectnms.genes_dict': 'genesdict',
        'ensembldb.release_37': 'ens37',
        'ensembldb.release_38': 'ens38',
        'mapnms.hgncid_to_symbol': 'hgncidtosymbol',
        'selectensts.gene_synonyms': 'genesynonyms',
        'formatcarts.series_37': 'series37',
        'formatcarts.series_38': 'series38',
        'formatcarts.canonical_37': 'canonical37',
        'formatcarts.canonical_38': 'canonical38',
        'output_prefix': 'prefix'
    }

    for k in map_of_keys.keys():
        if k in config:
            ret[map_of_keys[k]] = config[k]

    return ret


def add_default_values(c, rootdir):

    default_values = {
        'apprisfile': '{}/default/appris_refseq107.txt'.format(rootdir),
        'genesdict': '{}/default/20180104_hgnc_biomart.txt'.format(rootdir),
        'hgncidtosymbol': '{}/default/20180104_HGNCID_to_GeneSymbol_from_HGNC_BioMart.txt'.format(rootdir),
        'ens37': '75',
        'ens38': '92',
        'genesynonyms': '{}/default/gene_synonyms.txt'.format(rootdir),
        'series37': 'CART37A',
        'series38': 'CART38A',
        'canonical37': '{}/default/canonical_enst_v75.txt'.format(rootdir),
        'canonical38': '{}/default/canonical_enst_v92.txt'.format(rootdir),
        'prefix': 'output'
    }

    for k, v in default_values.iteritems():
        if k not in c or c[k] == '.':
            c[k] = v

    return c


def run(options):

    full = os.path.dirname(os.path.realpath(__file__))
    rootdir = full[:full.rfind('/env')]

    config = parse_config_file(options.config)
    mapped_args = map_keys(config)
    mapped_args = add_default_values(mapped_args, rootdir)

    args_to_pass = []
    for k, v in mapped_args.iteritems():
        args_to_pass.append('{}:{}'.format(k, v))

    args_to_pass.append('conffn:{}'.format(options.config))

    subprocess.call(['{}/pipeline/pipeline.sh'.format(rootdir)] + args_to_pass)
