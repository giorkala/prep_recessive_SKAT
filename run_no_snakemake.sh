#!/bin/bash

## parameters ##
min_AC=5 # minimum number of carriers per gene
max_AF=0.05 # maximum allele frequency per variant
TAG="brava_s50.delta" # tag to identify the output files
CHR="chr$1" # chromosome index

# software requirements
module load /FIXTHIS/bcftools/1.0
call_chets="/FIXTHIS/call_chets/call_chets"
prep_genemap="python prepare_genemap.py"
prep_new_encoding="python prepare_new_rec_enc.py"
generate_vcf="python generate_vcf.py"

## General file paths and prefices ##
WORK_DIR="/FIXTHIS/PATH"
input_bcf="FIXTHIS/PATH/cohort.phased_all.$CHR.bcf"
input_annot_saige="FIXTHIS/PATH/cohort.brava_s50_nopopmax.$CHR.saige.txt"
# this is the output from BRaVa's variant annotation
input_maf="FIXTHIS/PATH/cohort.notrios.$CHR.snpinfo"
# this is obtained during the phasing process

input_samples="${WORK_DIR}/samples.bcf.txt"
# get a list of samples in the BCF file
if [ ! -f $input_samples ]; then
    bcftools query --list-samples $input_bcf > $input_samples
fi

output_prefix="${WORK_DIR}/output/${TAG}.$CHR"

mkdir -p $WORK_DIR/output $WORK_DIR/logs $WORK_DIR/sandbox

## 1. extract all genotypes from the phased data
genotypes="$WORK_DIR/sandbox/$CHR.PP90af05.txt.gz"
if [ ! -f $genotypes ]; then
    echo "Extracting biallelic genotypes for $CHR..."
    bcftools view ${input_bcf} --max-af ${max_AF} | bcftools query -i '(GT="1|1" | GT="0|1" | GT="1|0")' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT %PP\n]' | awk '(length($2)<50)' | awk '($4=="." || $4>0.90)' | gzip > $genotypes
else
    echo "Genotypes already extracted, moving on..."
fi

## 2. prepare the annotation file for the genotypes
biallelic_genotypes=${WORK_DIR}/sandbox/biallelic.${TAG}.merged.${CHR}.txt
var_annot=${WORK_DIR}/sandbox/annotation.${TAG}.${CHR}.txt

$prep_genemap -a $input_annot_saige \
    -c pLoF damaging_missense_or_protein_altering other_missense_or_protein_altering \
    -o $var_annot
$call_chets -g $genotypes -m $var_annot --show-variants | grep 'chet\|hom' > $biallelic_genotypes

# repeat for synonymous variants, adding to the existing file
annot_syn="${WORK_DIR}/sandbox/annotation.syn.${TAG}.${CHR}.txt"
$prep_genemap -a $input_annot_saige -c synonymous -o $annot_syn
$call_chets -g $genotypes -m $annot_syn --show-variants | grep 'chet\|hom' >> $biallelic_genotypes
cat $annot_syn >> $var_annot

## 3. generate the recessive encoding using the above as input
$prep_new_encoding -c $1 -t merged -g $biallelic_genotypes -m $input_maf -a $var_annot -o $output_prefix

## 4. generate a VCF files for the current chromosome
$generate_vcf -s $input_samples -c $1 -g $output_prefix.txt -o $output_prefix.vcf.gz

echo "All done for $CHR!"
rm $annot_syn $var_annot $biallelic_genotypes

# All done! 