import re
import time
import os
import gzip


class CART(object):
    def __init__(self, hgnc_id, nm, nm_version, principal, utrDif, utrCriteria, genomeDif, alignment):
        self.hgnc_id = hgnc_id
        self.cart_id = '?' if nm_version != '.' else '.'
        self.nm = nm
        self.nm_version = nm_version
        self.principal = principal
        self.utrDif = utrDif
        self.utrCriteria = utrCriteria
        self.genomeDif = genomeDif
        self.alignment = alignment
        self.symbol = "."
        self.gdm_colour = "."
        self.ncbi_gene = "."
        self.nms_in_refseq = "."
        self.nms_in_appris = "."
        self.icr100 = "."
        self.refseqgene = "."
        self.clinvar = "."
        self.lrg = "."
        self.lovd = "."
        self.community = "."
        self.community_appris = "."
        self.tgmi_vs_rsg = "."
        self.tgmi_vs_clinvar = "."
        self.tgmi_vs_lrg = "."
        self.tgmi_vs_lovd = "."
        self.tgmi_vs_community = "."
        self.automatedcolour = "."
        self.communitycolour = "."
        self.final_nm = "."
        self.final_colour_group = "."
        self.all_nms_in_refseq = "."
        self.all_xms_in_refseq = "."

    def return_string(self):
        return (
        "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
        self.hgnc_id,self.nm, self.nm_version, self.principal, self.utrDif, self.utrCriteria,
        self.genomeDif, self.alignment, self.symbol, self.gdm_colour, self.ncbi_gene, self.nms_in_refseq,
        self.nms_in_appris, self.all_nms_in_refseq, self.all_xms_in_refseq, self.icr100, self.refseqgene, self.clinvar,
        self.lrg, self.lovd, self.community, self.community_appris, self.tgmi_vs_rsg, self.tgmi_vs_clinvar,
        self.tgmi_vs_lrg, self.tgmi_vs_lovd, self.tgmi_vs_community, self.automatedcolour, self.communitycolour,
        self.final_nm, self.final_colour_group))


def open_cart_selection(carts, infile):

    for line in open(infile, 'r'):
        if line.startswith("#"): continue
        if line == "": continue
        line = line.rstrip()
        fields = line.split("\t")
        carts[fields[0]] = CART(fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[6], fields[7])


def open_cart_missing(carts, infile):

    for line in open(infile, 'r'):
        if line.startswith("#"): continue
        if line == "": continue
        line = line.rstrip()
        fields = line.split("\t")
        if fields[0] in carts:
            print "WARNING: HGNC %s was in both the CART selection and missing files" % (fields[0])
        carts[fields[0]] = CART(fields[0], fields[1], ".", ".", ".", ".", ".", ".")


def open_symbol_to_hgnc_file(symbol_to_hgnc):

    for line in open("{}/default/20180104_HGNCID_to_GeneSymbol_from_HGNC_BioMart.txt".format(rootdir), 'r'):
        line = line.rstrip()
        if line == "": continue
        if line.startswith("#"): continue
        fields = line.split("\t")
        symbol_to_hgnc[fields[1]] = fields[0]


def open_ncbi_to_hgnc_file(ncbi_to_hgnc):

    for line in open("{}/default/20180104_hgnc_biomart.txt".format(rootdir), 'r'):
        line = line.rstrip()
        if line == "": continue
        if line.startswith("HGNC ID"): continue
        fields = line.split("\t")
        if len(fields) > 1:
            ncbi_to_hgnc[fields[1]] = fields[0]


def open_refseqGene_to_gene(infile, ncbi_to_hgnc, dict):

    for line in open(infile, 'r'):
        line = line.rstrip()
        if line == "": continue
        if line.startswith("#"): continue
        fields = line.split("\t")
        rsg = fields[3].split(".", 1)[0]
        if fields[1] in ncbi_to_hgnc: dict[rsg] = ncbi_to_hgnc[fields[1]]


def fill_in_gdm_color(infile, carts):
    for line in open(infile, 'r'):
        line = line.rstrip()
        if line == "": continue
        if line.startswith("#"): continue
        fields = line.split("\t")
        if len(fields) != 2:
            continue
        if fields[0] in carts: carts[fields[0]].gdm_colour = fields[1]


def fill_in_refseq_nms(infile, carts, nm_to_hgnc):

    refseq_nms = {}
    for line in open(infile, 'r'):
        line = line.rstrip()
        if line == "": continue
        if line.startswith("#"): continue
        fields = line.split("\t")
        if fields[1] not in refseq_nms: refseq_nms[fields[1]] = []
        refseq_nms[fields[1]].append(fields[0])
        nm_to_hgnc[fields[0]] = fields[1]

    for hgnc in refseq_nms.keys():
        if hgnc in carts: carts[hgnc].nms_in_refseq = ", ".join(list(set(refseq_nms[hgnc])))


def open_appris_file(infile, appris_hgnc_to_nms, appris_nms_to_princ, ncbi_to_hgnc):

    for line in open(infile, 'r'):
        line = line.rstrip()
        if line == "": continue
        if line.startswith("#"): continue
        fields = line.split("\t")
        nm = fields[2].split(".", 1)[0]
        if fields[1] in ncbi_to_hgnc:
            hgnc = ncbi_to_hgnc[fields[1]]
            if hgnc not in appris_hgnc_to_nms: appris_hgnc_to_nms[hgnc] = []
            if nm not in appris_hgnc_to_nms[hgnc]: appris_hgnc_to_nms[hgnc].append(nm)
            appris_nms_to_princ[nm] = fields[4]


def open_icr100_genes_selection(symbol_to_hgnc, carts):

    for line in open("{}/default/MCG_transcripts_NMstem.txt".format(rootdir), 'r'):
        line = line.rstrip()
        if line.startswith("Gene"): continue
        fields = line.split("\t")
        symbol = fields[0]
        nm = fields[1]
        if symbol in symbol_to_hgnc and symbol_to_hgnc[symbol] in carts: carts[symbol_to_hgnc[symbol]].icr100 = nm


def open_lrgs_selection(nm_to_hgnc, symbol_to_hgnc, refseqgene_to_hgnc, carts, missing_file):

    lrg_nms = {}
    for line in open("{}/default/20170502_lrg_list_LRGs_transcripts_xrefs.txt".format(rootdir), 'r'):
        line = line.rstrip()
        if line.startswith("#"): continue
        fields = line.split("\t")
        symbol = fields[1]
        nm = fields[4].split(".", 1)[0]
        rsg = fields[2].split(".", 1)[0]

        hgnc = "."
        if nm in nm_to_hgnc:
            hgnc = nm_to_hgnc[nm]
        elif rsg in refseqgene_to_hgnc:
            hgnc = refseqgene_to_hgnc[rsg]
        elif symbol in symbol_to_hgnc:
            hgnc = symbol_to_hgnc[symbol]

        if hgnc == ".":
            missing_file.write("WARNING: LRG missing HGNC ID for %s\n" % (line))
        else:
            if hgnc not in lrg_nms: lrg_nms[hgnc] = []
            lrg_nms[hgnc].append(nm)

    for hgnc in lrg_nms.keys():
        if hgnc in carts: carts[hgnc].lrg = ", ".join(list(set(lrg_nms[hgnc])))


def open_refseqgene_selection(ncbi_to_hgnc, nm_to_hgnc, refseqgene_to_hgnc, symbol_to_hgnc, carts, missing_file):

    rsg_nms = {}
    for line in open("{}/default/20170502_refseqgene_LRG_RefSeqGene.txt".format(rootdir), 'r'):
        line = line.rstrip()
        if line.startswith("#"): continue
        fields = line.split("\t")
        ncbi = fields[1]
        symbol = fields[2]
        nm = fields[5].split(".", 1)[0]
        rsg = fields[3].split(".", 1)[0]
        category = fields[9]
        if category != "reference standard": continue
        if not (nm.startswith("NM_")): continue

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
            missing_file.write("WARNING: RSG missing HGNC ID for %s\n" % (line))
        else:
            if hgnc not in rsg_nms: rsg_nms[hgnc] = []
            rsg_nms[hgnc].append(nm)

    for hgnc in rsg_nms.keys():
        if hgnc in carts: carts[hgnc].refseqgene = ", ".join(list(set(rsg_nms[hgnc])))


def open_clinvar_selection(ncbi_to_hgnc, nm_to_hgnc, symbol_to_hgnc, carts, missing_file):

    clinvar_nms = {}
    # Fix to only use NM associated with variant
    for line in gzip.open("{}/default/20170502_clinvar_variant_summary_extract.txt.gz".format(rootdir), 'r'):
        line = line.rstrip()
        if line.startswith("#"): continue
        fields = line.split("\t")
        name = fields[0]
        ncbi = fields[1]
        symbol = fields[2]
        assembly = fields[3]

        if assembly != "GRCh37": continue
        if not name.startswith("NM_"): continue

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

        # Hard coded a fix for this HGNC ID.  The NCBI_ID that used to go to HGNC:18432 now points to HGNC:7057.
        if hgnc == "HGNC:18432": hgnc = "HGNC:7057"

        if hgnc == "-":
            missing_file.write("WARNING: ClinVar missing HGNC ID for %s\n" % (line))
        else:
            if hgnc not in clinvar_nms: clinvar_nms[hgnc] = []
            clinvar_nms[hgnc].append(nm)

    for hgnc in clinvar_nms.keys():
        if hgnc in carts: carts[hgnc].clinvar = ", ".join(list(set(clinvar_nms[hgnc])))


def open_lovd_selection(nm_to_hgnc, symbol_to_hgnc, carts, missing_file):

    lovd_nms = {}
    for line in open("{}/default/20170504_LOVD_Transcripts.txt".format(rootdir), 'r'):
        line = line.rstrip()
        if line.startswith("#"): continue
        line = line.replace('"', '')
        if line.startswith('{{chromosome'): continue
        fields = line.split("\t")
        symbol = fields[1]
        nm = fields[3].split(".", 1)[0]

        if not nm.startswith("NM_"): continue

        hgnc = "."
        if nm in nm_to_hgnc:
            hgnc = nm_to_hgnc[nm]
        elif symbol in symbol_to_hgnc:
            hgnc = symbol_to_hgnc[symbol]

        if hgnc == ".":
            missing_file.write("WARNING: LOVD missing HGNC ID for %s\n" % line)
        else:
            if hgnc not in lovd_nms: lovd_nms[hgnc] = []
            lovd_nms[hgnc].append(nm)

    for hgnc in lovd_nms.keys():
        if hgnc in carts: carts[hgnc].lovd = ", ".join(list(set(lovd_nms[hgnc])))


def select_community_nm(cart):

    if "," not in cart.refseqgene and cart.refseqgene != ".":
        return cart.refseqgene
    elif cart.refseqgene == "." and "," not in cart.clinvar and cart.clinvar != ".":
        return cart.clinvar
    elif cart.refseqgene == "." and cart.clinvar == "." and "," not in cart.lrg and cart.lrg != ".":
        return cart.lrg
    else:
        return "."


def open_all_nm_and_xms_in_refseq(ncbi_to_hgnc, nms_given_hgnc, xms_given_hgnc):

    for line in gzip.open("{}/default/20171211_refSeq_release85.accession2geneid_extract.txt.gz".format(rootdir), 'r'):
        line = line.rstrip()
        (tax, NCBI_ID, transcript_id, protien_id) = line.split("\t")
        if tax != "9606": continue
        transcript_id = transcript_id.split(".")[0]
        if not (transcript_id.startswith("NM_") or transcript_id.startswith("XM_")): continue
        if NCBI_ID not in ncbi_to_hgnc: continue
        hgnc_id = ncbi_to_hgnc[NCBI_ID]
        if hgnc_id not in nms_given_hgnc:
            nms_given_hgnc[hgnc_id] = []
            xms_given_hgnc[hgnc_id] = []

        if transcript_id.startswith("NM_"):
            nms_given_hgnc[hgnc_id].append(transcript_id)
        elif transcript_id.startswith("XM_"):
            xms_given_hgnc[hgnc_id].append(transcript_id)




#############################################€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€###########################################


full = os.path.dirname(os.path.realpath(__file__))
rootdir = full[:full.rfind('/env')]


def run(options, cart_selection_fn, cart_missing_fn):

    # Initialize the CARTs
    carts = {}
    open_cart_selection(carts, cart_selection_fn)
    open_cart_missing(carts, cart_missing_fn)

    # Generate dictionaries to convert to HGNCI IDs
    symbol_to_hgnc = {}
    ncbi_to_hgnc = {}
    refseqgene_to_hgnc = {}
    open_symbol_to_hgnc_file(symbol_to_hgnc)
    open_ncbi_to_hgnc_file(ncbi_to_hgnc)
    for ncbi_gene in ncbi_to_hgnc.keys():
        if ncbi_to_hgnc[ncbi_gene] in carts: carts[ncbi_to_hgnc[ncbi_gene]].ncbi_gene = ncbi_gene
    open_refseqGene_to_gene("{}/default/20171026_gene_RefSeqGene.txt".format(rootdir),
                            ncbi_to_hgnc, refseqgene_to_hgnc)

    # Assign the symbol to the carts
    for symbol in symbol_to_hgnc.keys():
        if symbol_to_hgnc[symbol] in carts: carts[symbol_to_hgnc[symbol]].symbol = symbol

    # fill in gene colors
    fill_in_gdm_color("{}/default/GDM_Score_20180712.txt".format(rootdir), carts)

    # Get set of NMs in RefSeq interim alignment file
    nm_to_hgnc = {}
    fill_in_refseq_nms(options.refsdbinc, carts, nm_to_hgnc)

    # Get appris info
    appris_hgnc_to_nms = {}
    appris_nms_to_princ = {}
    open_appris_file("{}/default/20170113_marton_appris_refseq107.txt".format(rootdir),
                     appris_hgnc_to_nms, appris_nms_to_princ, ncbi_to_hgnc)

    # get nms in refseq
    nms_given_hgnc = {}
    xms_given_hgnc = {}
    open_all_nm_and_xms_in_refseq(ncbi_to_hgnc, nms_given_hgnc, xms_given_hgnc)

    # read in community databases
    missing_file = open('{}_missingCommunityGenes.txt'.format(options.out_final), 'w')
    open_icr100_genes_selection(symbol_to_hgnc, carts)
    open_lrgs_selection(nm_to_hgnc, symbol_to_hgnc, refseqgene_to_hgnc, carts, missing_file)
    open_refseqgene_selection(ncbi_to_hgnc, nm_to_hgnc, refseqgene_to_hgnc, symbol_to_hgnc, carts, missing_file)
    open_clinvar_selection(ncbi_to_hgnc, nm_to_hgnc, symbol_to_hgnc, carts, missing_file)
    open_lovd_selection(nm_to_hgnc, symbol_to_hgnc, carts, missing_file)
    missing_file.close()

    # Add appris info
    for hgnc in appris_hgnc_to_nms.keys():
        if hgnc in carts:
            cur_ap = ""
            for nm in appris_hgnc_to_nms[hgnc]:
                if cur_ap == "":
                    cur_ap = nm + " (" + appris_nms_to_princ[nm] + ")"
                else:
                    cur_ap += ", " + nm + " (" + appris_nms_to_princ[nm] + ")"
            carts[hgnc].nms_in_appris = cur_ap

    # Add community and appris
    for hgnc in carts:
        cur_com_nms = []
        if carts[hgnc].refseqgene != ".": cur_com_nms = cur_com_nms + (carts[hgnc].refseqgene.split(", "))
        if carts[hgnc].clinvar != ".": cur_com_nms = cur_com_nms + carts[hgnc].clinvar.split(", ")
        if carts[hgnc].lrg != ".": cur_com_nms = cur_com_nms + carts[hgnc].lrg.split(", ")
        if carts[hgnc].lovd != ".": cur_com_nms = cur_com_nms + carts[hgnc].lovd.split(", ")
        cur_com_nms = list(set(cur_com_nms))

        if len(cur_com_nms) != 0:
            cur_app_nms = ""
            carts[hgnc].community = ", ".join(cur_com_nms)
            for nm in cur_com_nms:
                if cur_app_nms == "":
                    if nm in appris_nms_to_princ:
                        cur_app_nms = nm + " (" + appris_nms_to_princ[nm] + ")"
                    else:
                        cur_app_nms = nm + " (.)"
                else:
                    if nm in appris_nms_to_princ:
                        cur_app_nms += ", " + nm + " (" + appris_nms_to_princ[nm] + ")"
                    else:
                        cur_app_nms += ", " + nm + " (.)"
            carts[hgnc].community_appris = cur_app_nms

    out_file = open('{}.txt'.format(options.out_final), 'w')
    out_file.write(
        "#HGNC_ID\tAutomatedRelatedTranscript\tAutomatedRelatedTranscriptVersion\tAutomatedRelatedTranscriptAPPRIS\tUTRDifferenceType\tUTRDecisiveCriteria\tGenomeDifference\tAlignmentDifference\tSymbol\tGDM_Colour\tNCBI.Gene\tNMs.in.RefSeq\tNMs.in.APPRIS\tNMs_in_all_of_RefSeq\tXMs_in_all_of_RefSeq\tICR100.NM\tRefSeqGene\tClinVar\tLRG\tLOVD\tCommunity\tCommunity_APPRIS\tTGMI_vs_RSG\tTGMI_vs_ClinVar\tTGMI_vs_LRG\tTGMI_vs_LOVD\tTGMI_vs_Community\tAutomatedColour\tCommunityColour\tCART_NM\tSetFlag\n")

    tmp_colours1 = ["green1", "green2", "grey", "."]
    tmp_colours2 = ["green1", "grey1", "grey2", "grey3", "grey4"]
    colour_table = {}
    auto_total_table = {}
    comm_total_table = {}
    final_flag_total = {}
    for c in tmp_colours1:
        colour_table[c] = {}
        auto_total_table[c] = 0
        for c2 in tmp_colours2:
            colour_table[c][c2] = 0
            comm_total_table[c2] = 0

    for hid in sorted(carts.keys()):
        cur_community = carts[hid].community.split(", ")
        cur_rs_nms = carts[hid].nms_in_refseq.split(", ")

        if hid in nms_given_hgnc:
            carts[hid].all_nms_in_refseq = ", ".join(nms_given_hgnc[hid])
        if hid in xms_given_hgnc:
            carts[hid].all_xms_in_refseq = ", ".join(xms_given_hgnc[hid])

        if carts[hid].community != ".":
            if carts[hid].nm in cur_community:
                carts[hid].tgmi_vs_community = "1"
            else:
                carts[hid].tgmi_vs_community = "0"

        if carts[hid].refseqgene != ".":
            if carts[hid].nm in carts[hid].refseqgene.split(", "):
                carts[hid].tgmi_vs_rsg = "1"
            else:
                carts[hid].tgmi_vs_rsg = "0"

        if carts[hid].clinvar != ".":
            if carts[hid].nm in carts[hid].clinvar.split(", "):
                carts[hid].tgmi_vs_clinvar = "1"
            else:
                carts[hid].tgmi_vs_clinvar = "0"

        if carts[hid].lrg != ".":
            if carts[hid].nm in carts[hid].lrg.split(", "):
                carts[hid].tgmi_vs_lrg = "1"
            else:
                carts[hid].tgmi_vs_lrg = "0"

        if carts[hid].lovd != ".":
            if carts[hid].nm in carts[hid].lovd.split(", "):
                carts[hid].tgmi_vs_lovd = "1"
            else:
                carts[hid].tgmi_vs_lovd = "0"

        if carts[hid].cart_id != ".":
            if len(cur_rs_nms) == 1 and (
                    carts[hid].principal == "PRINCIPAL:1" or carts[hid].principal == "PRINCIPAL:2") and \
                            carts[hid].nm == carts[hid].nms_in_refseq:
                carts[hid].automatedcolour = "green1"
            elif (carts[hid].principal == "PRINCIPAL:1" or carts[hid].principal == "PRINCIPAL:2") and carts[
                hid].nm in cur_rs_nms:
                carts[hid].automatedcolour = "green2"
            else:
                carts[hid].automatedcolour = "grey"

        if carts[hid].community == ".":
            carts[hid].communitycolour = "grey4"
        elif carts[hid].community == carts[hid].nm:
            carts[hid].communitycolour = "green1"
        elif carts[hid].nm in cur_community:
            carts[hid].communitycolour = "grey1"
        elif len(cur_community) == 1:
            carts[hid].communitycolour = "grey2"
        else:
            carts[hid].communitycolour = "grey3"

        if carts[hid].communitycolour == "green1" and (
                            carts[hid].automatedcolour == "green1" or carts[hid].automatedcolour == "green2" or carts[
                    hid].automatedcolour == "grey"):
            carts[hid].final_colour_group = "Algorithmic"
            carts[hid].final_nm = carts[hid].nm
        elif carts[hid].communitycolour == "grey1" and carts[hid].automatedcolour == "green1":
            carts[hid].final_colour_group = "Algorithmic"
            carts[hid].final_nm = carts[hid].nm
        elif carts[hid].communitycolour == "grey1" and (
                        carts[hid].automatedcolour == "green2" or carts[hid].automatedcolour == "grey"):
            if carts[hid].nm == carts[hid].refseqgene or carts[hid].nm in carts[hid].refseqgene or (
                            carts[hid].refseqgene == "." and (
                                    carts[hid].nm == carts[hid].clinvar or carts[hid].nm in carts[hid].clinvar)) or (
                                carts[hid].refseqgene == "." and carts[hid].clinvar == "." and (
                                    carts[hid].nm == carts[hid].clinvar or carts[hid].nm in carts[hid].clinvar)) or (
                                carts[hid].refseqgene == "." and carts[hid].clinvar == "." and carts[hid].lrg == "."):
                carts[hid].final_colour_group = "AlgorithmicFlag"
                carts[hid].final_nm = carts[hid].nm
            else:
                carts[hid].final_colour_group = "CommunityChoice"
                carts[hid].final_nm = select_community_nm(carts[hid])
        elif carts[hid].communitycolour == "grey2" and carts[hid].automatedcolour == ".":
            if carts[hid].refseqgene != "." or carts[hid].clinvar != "." or carts[hid].lrg != ".":
                carts[hid].final_colour_group = "Community"
                carts[hid].final_nm = carts[hid].community
            else:
                carts[hid].final_colour_group = "NoData"
                carts[hid].final_nm = "."
        elif carts[hid].communitycolour == "grey2" and carts[hid].automatedcolour == "green1":
            carts[hid].final_colour_group = "Algorithmic"
            carts[hid].final_nm = carts[hid].nm
        elif carts[hid].communitycolour == "grey2" and (
                        carts[hid].automatedcolour == "green2" or carts[hid].automatedcolour == "grey"):
            if (carts[hid].refseqgene == "." and carts[hid].clinvar == "." and carts[hid].lrg == "."):
                carts[hid].final_colour_group = "AlgorithmicFlag"
                carts[hid].final_nm = carts[hid].nm
            else:
                carts[hid].final_colour_group = "CommunityFlag"
                carts[hid].final_nm = carts[hid].community
        elif carts[hid].communitycolour == "grey3" and (
                            carts[hid].automatedcolour == "." or carts[hid].automatedcolour == "green2" or carts[
                    hid].automatedcolour == "grey"):
            if carts[hid].refseqgene != "." or carts[hid].clinvar != "." or carts[hid].lrg != ".":
                carts[hid].final_colour_group = "CommunityChoice"
                carts[hid].final_nm = select_community_nm(carts[hid])
            else:
                carts[hid].final_colour_group = "NoData"
                carts[hid].final_nm = "."
        elif carts[hid].communitycolour == "grey4" and (
                            carts[hid].automatedcolour == "green1" or carts[hid].automatedcolour == "green2" or carts[
                    hid].automatedcolour == "grey"):
            carts[hid].final_colour_group = "Algorithmic"
            carts[hid].final_nm = carts[hid].nm
        elif carts[hid].communitycolour == "grey4" and carts[hid].automatedcolour == ".":
            carts[hid].final_colour_group = "NoData"
            carts[hid].final_nm = "."
        else:
            print "missing %s and %s" % (carts[hid].automatedcolour, carts[hid].communitycolour)

        out_file.write(carts[hid].return_string() + "\n")
        colour_table[carts[hid].automatedcolour][carts[hid].communitycolour] += 1
        auto_total_table[carts[hid].automatedcolour] += 1
        comm_total_table[carts[hid].communitycolour] += 1
        if carts[hid].final_colour_group not in final_flag_total: final_flag_total[carts[hid].final_colour_group] = 0
        final_flag_total[carts[hid].final_colour_group] += 1

    out_file.close()

    counts_out_file = open('{}_counts.txt'.format(options.out_final), 'w')
    for c in tmp_colours1:
        counts_out_file.write("\t" + c)
    counts_out_file.write("\tTotal\n")

    for c2 in tmp_colours2:
        counts_out_file.write(c2)
        for c in tmp_colours1:
            counts_out_file.write("\t" + str(colour_table[c][c2]))
        counts_out_file.write("\t" + str(comm_total_table[c2]) + "\n")

    counts_out_file.write("Total")
    for c in tmp_colours1:
        counts_out_file.write("\t" + str(auto_total_table[c]))
    counts_out_file.write("\t18885\n")

    counts_out_file.write("\n\n")
    for f in final_flag_total.keys():
        counts_out_file.write("%s\t%s\n" % (f, final_flag_total[f]))
    counts_out_file.close()