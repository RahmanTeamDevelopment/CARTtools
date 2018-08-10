#!/usr/bin/env bash

green=`tput setaf 2`
reset=`tput sgr0`

msg () {
echo ""
echo "${green}>>> CART Pipeline: $1${reset}"
}

# Read in arguments and defining variables
for var in "$@"
do
	IFS=':' read -ra ARR <<< "$var"
    key=${ARR[0]}
	value=${ARR[1]}
	declare $key=$value
done

echo ""

# Absolute path to CARTtools
ROOTPATH=$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)


# Run RefSeqDB (NCBI interim RefSeq mapping)
#msg "Running RefSeqDB (NCBI interim RefSeq mapping)"
#$ROOTPATH/refseqdb --build GRCh37 --mapping ncbi --output refseqdb_GRCh37_ncbi

# Run RefSeqDB (UCSC mapping)
#msg "Running RefSeqDB (UCSC mapping)"
#$ROOTPATH/refseqdb --build GRCh37 --mapping ucsc --output refseqdb_GRCh37_ucsc


# Run RefSeqScan (NCBI mapping)
#msg "Running RefSeqScan (NCBI mapping)"
#$ROOTPATH/refseqscan --input refseqdb_GRCh37_ncbi.gz --reference $ref37 --output refseqscan_output_GRCh37_ncbi.txt

# Run RefSeqScan (UCSC mapping)
#msg "Running RefSeqScan (UCSC mapping)"
#$ROOTPATH/refseqscan --input refseqdb_GRCh37_ucsc.gz --reference $ref37 --output refseqscan_output_GRCh37_ucsc.txt


# Run SelectNMs
#msg "Running SelectNMs"
#$ROOTPATH/selectnms --input_genes $inputgenes --appr $apprisfile --refsdb refseqdb_GRCh37_ncbi.gz --refss refseqscan_output_GRCh37_ncbi.txt --genes $genesdict --build GRCh37 --series CART37A --out selected_nms_GRCh37


# Run MapNMs (GRCh37)
#msg "Running MapNMs (GRCh37)"
#$ROOTPATH/mapnms --input selected_nms_GRCh37_selected.txt --ncbi refseqdb_GRCh37_ncbi.gz --ucsc refseqdb_GRCh37_ucsc.gz --hgncid_to_symbol $hgncidtosymbol --output mapped_nms_GRCh37


# Run EnsemblDB (GRCh37)
#msg "Running EnsemblDB (GRCh37)"
#$ROOTPATH/ensembldb --release $ens37 --output ensembl_$ens37

# Run EnsemblDB (GRCh38)
#msg "Running EnsemblDB (GRCh38)"
#$ROOTPATH/ensembldb --release $ens38 --output ensembl_$ens38


# Run SelectENSTs (GRCh37)
msg "Running SelectENSTs (GRCh37)"
$ROOTPATH/selectensts --mapped_nms mapped_nms_GRCh37.gz --ensembl_data ensembl_$ens37.gz --gene_synonyms $genesynonyms --output selected_ensts_GRCh37


# Run FormatCARTs (GRCh37)
msg "Running FormatCARTs (GRCh37)"
$ROOTPATH/formatcarts --input selected_ensts_GRCh37.txt --ensembl ensembl_$ens37.gz --series $series --ref $ref37 --output $output