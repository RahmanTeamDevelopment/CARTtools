#!/usr/bin/env bash

green=`tput setaf 2`
cyan=`tput setaf 6`
reset=`tput sgr0`

info () {
d=$(date)
echo "${cyan}>>> $1 ($d)${reset}"
}

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

########################################################################################################################

echo ""
info "CART Pipeline started"
echo ""

# Start date/time
start=$(python -c "from datetime import datetime; print(datetime.now().strftime('%Y-%m-%d_%H:%M:%S'))")

# Absolute path to CARTtools
ROOTPATH=$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)

# Output directory for intermediate files
mkdir -p cart_pipeline_files
od="cart_pipeline_files"

# Run RefSeqDB (NCBI interim RefSeq mapping) for GRCh37 and GRCh38
msg "Running RefSeqDB (NCBI interim RefSeq mapping) for GRCh37"
$ROOTPATH/refseqdb --build GRCh37 --mapping ncbi --output ${od}/${prefix}_refseqdb_GRCh37_ncbi
msg "Running RefSeqDB (NCBI interim RefSeq mapping) for GRCh38"
$ROOTPATH/refseqdb --build GRCh38 --mapping ncbi --output ${od}/${prefix}_refseqdb_GRCh38_ncbi

# Run RefSeqDB (UCSC mapping) for GRCh37 and GRCh38
msg "Running RefSeqDB (UCSC mapping) for GRCh37"
$ROOTPATH/refseqdb --build GRCh37 --mapping ucsc --output ${od}/${prefix}_refseqdb_GRCh37_ucsc
msg "Running RefSeqDB (UCSC mapping) for GRCh38"
$ROOTPATH/refseqdb --build GRCh38 --mapping ucsc --output ${od}/${prefix}_refseqdb_GRCh38_ucsc

# Run RefSeqCheck (NCBI mapping) for GRCh37
msg "Running RefSeqCheck (NCBI mapping) for GRCh37"
$ROOTPATH/refseqcheck --input ${od}/${prefix}_refseqdb_GRCh37_ncbi.gz --output ${od}/${prefix}_refseqcheck_output_GRCh37_ncbi.txt \
--reference $ref37

# Run SelectNMs
msg "Running SelectNMs"
$ROOTPATH/selectnms --input_genes $inputgenes --appr $apprisfile --refsdb ${od}/${prefix}_refseqdb_GRCh37_ncbi.gz \
--refschk ${od}/${prefix}_refseqcheck_output_GRCh37_ncbi.txt --genes $genesdict --build GRCh37 \
--out_auto ${od}/${prefix}_nms_GRCh37 --refsdbinc ${od}/${prefix}_refseqdb_GRCh37_ncbi_included.txt \
--out ${od}/${prefix}_nms_GRCh37_final_selected

# Run MapNMs for GRCh37 and GRCh38
msg "Running MapNMs for GRCh37"
$ROOTPATH/mapnms --input ${od}/${prefix}_nms_GRCh37_final_selected.txt --ncbi ${od}/${prefix}_refseqdb_GRCh37_ncbi.gz \
--ucsc ${od}/${prefix}_refseqdb_GRCh37_ucsc.gz --hgncid_to_symbol $hgncidtosymbol --output ${od}/${prefix}_mapped_nms_GRCh37
msg "Running MapNMs for GRCh38"
$ROOTPATH/mapnms --input ${od}/${prefix}_nms_GRCh37_final_selected.txt --ncbi ${od}/${prefix}_refseqdb_GRCh38_ncbi.gz \
--ucsc ${od}/${prefix}_refseqdb_GRCh38_ucsc.gz --hgncid_to_symbol $hgncidtosymbol --output ${od}/${prefix}_mapped_nms_GRCh38

# Run EnsemblDB for GRCh37 and GRCh38
msg "Running EnsemblDB for GRCh37"
$ROOTPATH/ensembldb --release $ens37 --output ${od}/${prefix}_ensembl_$ens37
msg "Running EnsemblDB for GRCh38"
$ROOTPATH/ensembldb --release $ens38 --output ${od}/${prefix}_ensembl_$ens38

# Run SelectENSTs for GRCh37 and GRCh38
msg "Running SelectENSTs for GRCh37"
$ROOTPATH/selectensts --mapped_nms ${od}/${prefix}_mapped_nms_GRCh37.gz --ensembl_data ${od}/${prefix}_ensembl_$ens37.gz \
--gene_synonyms $genesynonyms --output ${od}/${prefix}_selected_ensts_GRCh37
msg "Running SelectENSTs for GRCh38"
$ROOTPATH/selectensts --mapped_nms ${od}/${prefix}_mapped_nms_GRCh38.gz --ensembl_data ${od}/${prefix}_ensembl_$ens38.gz \
--gene_synonyms $genesynonyms --output ${od}/${prefix}_selected_ensts_GRCh38

# Run FormatCARTs for GRCh37 and GRCh38
msg "Running FormatCARTs for GRCh37"
$ROOTPATH/formatcarts --input ${od}/${prefix}_selected_ensts_GRCh37.txt --ensembl ${od}/${prefix}_ensembl_$ens37.gz \
--series $series37 --ref $ref37 --output ${prefix}_GRCh37
msg "Running FormatCARTs for GRCh38"
$ROOTPATH/formatcarts --input ${od}/${prefix}_selected_ensts_GRCh38.txt --ensembl ${od}/${prefix}_ensembl_$ens38.gz \
--series $series38 --ref $ref38 --output ${prefix}_GRCh38 --prev_cava_db ${prefix}_GRCh37_cava.gz --prev_ref $ref37

# Run CompareENSTs
msg "Running CompareENSTs"
$ROOTPATH/compareensts --enstsx ${od}/${prefix}_selected_ensts_GRCh37.txt --enstsy ${od}/${prefix}_selected_ensts_GRCh38.txt \
--datax ${od}/${prefix}_ensembl_$ens37.gz --datay ${od}/${prefix}_ensembl_$ens38.gz --ref37  $ref37 --ref38  $ref38 \
--output ${od}/${prefix}_compared_ensts.txt --input $inputgenes

# End date/time
end=$(python -c "from datetime import datetime; print(datetime.now().strftime('%Y-%m-%d_%H:%M:%S'))")

# Run Summarize
msg "Running Summarize"
$ROOTPATH/summarize --prefix ${od}/${prefix} --output $prefix --start $start --end $end --inputfn $inputgenes \
--configfn $conffn --ens37 $ens37 --ens38 $ens38

echo ""
info "CART Pipeline finished"
echo ""