import re
import time
import os
import gzip


class Gene(object):

    def __init__(self, hgnc_id, nm, nm_version, utr_diff, utr_criteria, genome_diff, alignment):

        self.hgnc_id = hgnc_id
        self.automated_nm = nm
        self.automated_nm_version = nm_version
        self.utr_diff = utr_diff
        self.utr_criteria = utr_criteria
        self.genome_diff = genome_diff
        self.alignment = alignment
        self.symbol = "."
        self.gdm_colour = "."
        self.ncbi_gene = "."
        self.nms_in_refseq = "."
        self.refseqgene = "."
        self.clinvar = "."
        self.community_nm = "."
        self.final_nm = "."


    def __str__(self):

        return '\t'.join(
            [
                self.symbol,
                self.hgnc_id,
                self.ncbi_gene,
                self.gdm_colour,
                self.automated_nm,
                self.automated_nm_version,
                self.utr_diff,
                self.utr_criteria,
                self.genome_diff,
                self.alignment,
                self.nms_in_refseq,
                self.refseqgene,
                self.clinvar,
                self.community_nm,
                self.final_nm
            ]
        )


def read_automated_selection(fn):

    ret = {}
    for line in open(fn, 'r'):
        if line.startswith("#") or line == '':
            continue
        line = line.rstrip()
        fields = line.split("\t")
        ret[fields[0]] = Gene(fields[0], fields[1], fields[2], fields[4], fields[5], fields[6], fields[7])
    return ret


def read_automated_missing_genes(genes, fn):

    for line in open(fn, 'r'):
        if line.startswith("#") or line == '':
            continue
        line = line.rstrip()
        fields = line.split("\t")
        if fields[0] in genes:
            print "WARNING: HGNC %s was in both the CART selection and missing files" % (fields[0])
        genes[fields[0]] = Gene(fields[0], ".", ".", ".", ".", ".", ".")


def read_symbol_to_hgnc_file():

    ret = {}
    for line in open("{}/default/20180104_HGNCID_to_GeneSymbol_from_HGNC_BioMart.txt".format(rootdir), 'r'):
        line = line.rstrip()
        if line.startswith("#") or line == '':
            continue
        fields = line.split("\t")
        ret[fields[1]] = fields[0]
    return ret


def read_ncbi_to_hgnc_file():

    ret = {}
    for line in open("{}/default/20180104_hgnc_biomart.txt".format(rootdir), 'r'):
        line = line.rstrip()
        if line.startswith('HGNC ID') or line == '':
            continue
        fields = line.split("\t")
        if len(fields) > 1:
            ret[fields[1]] = fields[0]
    return ret


def read_refseqgene_to_gene(fn, ncbi_to_hgnc):

    ret = {}
    for line in open(fn, 'r'):
        line = line.rstrip()
        if line.startswith("#") or line == '':
            continue
        fields = line.split("\t")
        rsg = fields[3].split(".", 1)[0]
        if fields[1] in ncbi_to_hgnc:
            ret[rsg] = ncbi_to_hgnc[fields[1]]
    return ret


def add_gdm_color(fn, genes):

    for line in open(fn, 'r'):
        line = line.rstrip()
        if line.startswith("#") or line == '':
            continue
        fields = line.split("\t")
        if fields[3] in genes:
            genes[fields[3]].gdm_colour = fields[2]


def read_refseq_nms(fn, genes):

    ret = {}

    refseq_nms = {}
    for line in open(fn, 'r'):
        line = line.rstrip()
        if line.startswith("#") or line == '':
            continue
        fields = line.split("\t")
        if fields[1] not in refseq_nms:
            refseq_nms[fields[1]] = []
        refseq_nms[fields[1]].append(fields[0])
        ret[fields[0]] = fields[1]

    for hgnc in refseq_nms.keys():
        if hgnc in genes:
            genes[hgnc].nms_in_refseq = ",".join(list(set(refseq_nms[hgnc])))

    return ret


def read_refseqgene_selection(ncbi_to_hgnc, nm_to_hgnc, refseqgene_to_hgnc, symbol_to_hgnc, genes, missing_file):

    rsg_nms = {}
    for line in open("{}/default/20170502_refseqgene_LRG_RefSeqGene.txt".format(rootdir), 'r'):
        line = line.rstrip()
        if line.startswith("#"):
            continue
        fields = line.split("\t")

        ncbi = fields[1]
        symbol = fields[2]
        nm = fields[5].split(".", 1)[0]
        rsg = fields[3].split(".", 1)[0]
        category = fields[9]

        if category != "reference standard":
            continue

        if not nm.startswith("NM_"):
            continue

        hgnc = "."
        if ncbi in ncbi_to_hgnc:
            hgnc = ncbi_to_hgnc[ncbi]
        elif nm in nm_to_hgnc:
            hgnc = nm_to_hgnc[nm]
        elif rsg in refseqgene_to_hgnc:
            hgnc = refseqgene_to_hgnc[rsg]
        elif symbol in symbol_to_hgnc:
            hgnc = symbol_to_hgnc[symbol]

        if hgnc == ".":
            missing_file.write("WARNING: RSG missing HGNC ID for %s\n" % line)
        else:
            if hgnc not in rsg_nms:
                rsg_nms[hgnc] = []
            rsg_nms[hgnc].append(nm)

    for hgnc in rsg_nms.keys():
        if hgnc in genes:
            genes[hgnc].refseqgene = ",".join(list(set(rsg_nms[hgnc])))


def read_clinvar_selection(ncbi_to_hgnc, nm_to_hgnc, symbol_to_hgnc, genes, missing_file):

    clinvar_nms = {}
    # Fix to only use NM associated with variant
    for line in gzip.open("{}/default/20170502_clinvar_variant_summary_extract.txt.gz".format(rootdir), 'r'):
        line = line.rstrip()
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        name = fields[0]
        ncbi = fields[1]
        symbol = fields[2]
        assembly = fields[3]

        if assembly != "GRCh37":
            continue
        if not name.startswith("NM_"):
            continue

        tmps = name.split(":")
        nm = re.sub("[\(\[].*?[\)\]]", "", tmps[0])
        var_symbol = "."
        if len(tmps[0].split("(")) > 1:
            var_symbol = re.sub("\)", "", tmps[0].split("(")[1])
        nm = nm.split(".", 1)[0]

        hgnc = "."
        if nm in nm_to_hgnc:
            hgnc = nm_to_hgnc[nm]
        elif var_symbol == symbol:
            hgnc = fields[4]
            if ncbi in ncbi_to_hgnc:
                hgnc = ncbi_to_hgnc[ncbi]
            elif symbol in symbol_to_hgnc:
                hgnc = symbol_to_hgnc[symbol]

        # Hard coded a fix for this HGNC ID. The NCBI_ID that used to go to HGNC:18432 now points to HGNC:7057.
        if hgnc == "HGNC:18432":
            hgnc = "HGNC:7057"

        if hgnc == "-":
            missing_file.write("WARNING: ClinVar missing HGNC ID for %s\n" % line)
        else:
            if hgnc not in clinvar_nms:
                clinvar_nms[hgnc] = []
            clinvar_nms[hgnc].append(nm)

    for hgnc in clinvar_nms.keys():
        if hgnc in genes:
            genes[hgnc].clinvar = ",".join(list(set(clinvar_nms[hgnc])))


def select_community_nm(gene):

    if gene.refseqgene != ".":
        gene.community_nm = gene.refseqgene
    else:
        if gene.clinvar != ".":
            gene.community_nm = gene.clinvar
        else:
            gene.community_nm = '.'


def select_final_nm(gene):

    community_nms = gene.community_nm.split(",")
    if gene.automated_nm != "." and gene.community_nm == '.':
        gene.final_nm = gene.automated_nm
    elif gene.automated_nm != "." and gene.automated_nm in community_nms:
        gene.final_nm = gene.automated_nm
    elif gene.automated_nm != "." and gene.automated_nm not in community_nms and len(community_nms) == 1:
        gene.final_nm = gene.community_nm
    elif gene.automated_nm == "." and len(community_nms) == 1:
        gene.final_nm = gene.community_nm
    else:
        gene.final_nm = '.'


########################################################################################################################


full = os.path.dirname(os.path.realpath(__file__))
rootdir = full[:full.rfind('/env')]

def run(options, automated_selection_fn, automated_missing_fn):

    # Initialize the Genes
    genes = read_automated_selection(automated_selection_fn)
    read_automated_missing_genes(genes, automated_missing_fn)

    # Generate dictionaries for converting to HGNC IDs
    symbol_to_hgnc = read_symbol_to_hgnc_file()
    ncbi_to_hgnc = read_ncbi_to_hgnc_file()
    for ncbi_gene in ncbi_to_hgnc.keys():
        if ncbi_to_hgnc[ncbi_gene] in genes:
            genes[ncbi_to_hgnc[ncbi_gene]].ncbi_gene = ncbi_gene
    refseqgene_to_hgnc = read_refseqgene_to_gene("{}/default/20171026_gene_RefSeqGene.txt".format(rootdir), ncbi_to_hgnc)

    # Assign gene symbols to the genes
    for symbol in symbol_to_hgnc.keys():
        if symbol_to_hgnc[symbol] in genes:
            genes[symbol_to_hgnc[symbol]].symbol = symbol

    # Add GDM colors
    add_gdm_color("{}/default/GDM_beta_201807.txt".format(rootdir), genes)

    # Get set of NMs in RefSeq interim alignment file
    nm_to_hgnc = read_refseq_nms(options.refsdbinc, genes)

    # Read community databases
    missing_file = open('{}_missingCommunityGenes.txt'.format(options.out_final), 'w')
    read_refseqgene_selection(ncbi_to_hgnc, nm_to_hgnc, refseqgene_to_hgnc, symbol_to_hgnc, genes, missing_file)
    read_clinvar_selection(ncbi_to_hgnc, nm_to_hgnc, symbol_to_hgnc, genes, missing_file)
    missing_file.close()

    out_file = open('{}.txt'.format(options.out_final), 'w')
    header = [
        'GeneSymbol',
        'HGNC_ID',
        'NCBI.Gene',
        'GDM_Colour',
        'Automated_NM',
        'Automated_NM_Version',
        'UTRDifferenceType',
        'UTRDecisiveCriteria',
        'GenomeDifference',
        'AlignmentDifference',
        'NMs.in.RefSeq',
        'RefSeqGene_NMs',
        'ClinVar_NMs',
        'Community_NM',
        'CART_NM',
    ]
    out_file.write('#{}\n'.format('\t'.join(header)))

    for hgnc_id in sorted(genes.keys()):
        gene = genes[hgnc_id]
        select_community_nm(gene)
        select_final_nm(gene)
        out_file.write('{}\n'.format(gene))

    out_file.close()
