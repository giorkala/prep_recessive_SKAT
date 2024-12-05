## Framework for Recessive Sequence Kernel Association Testing (ReSKAT) ##
Genotype preparation based on pseudo-variants for biallelic genotypes and association testing with SAIGE. Remember to change all required paths in the `Snakefile` and make sure all software is loaded appropriately.

### Overview of framework
This framework is developed as a Snakemake pipeline to streamline the preparation, with the following tasks:
* Step 1: Detect all biallelic genotypes, obtain variant annotation and allele frequencies. Make sure [call_chets](https://github.com/frhl/call_chets) is installed and can be referenced.
* Step 2: `prepare_new_rec_enc.py`
    1. For each individual, select the genotype with the worst consequence and assign an identifier
    2. Extract the frequency of the corresponding variant(s)
    3. Assign weights and generate annotation files for SAIGE/Regenie
* Step 3: `generate_vcf.py`; this is a simple way to create a VCF file based on the genotypes obtained previously, one per chromosome.
* Step 4: Concatenate all VCFs to one file, then convert to BGEN for higher efficiency.

Requirements:
* phased whole exome/genome sequencing data (post-QC)
* variant annotation
* allele frequencies
* ranking of paired consequences, e.g. "pLOF|damaging" > "damaging|damaging"
* collapsing criteria and weighting scheme for rare biallelic genotypes
    * any collapsing should follow the extraction of AF products, as we'll then need some sort of average
* SAIGE step1 files (for association testing)

### How to run
```
git clone --recurse-submodules git@github.com:giorkala/prep_recessive_SKAT.git
# then load module for snakemake
# then edit paths in Snakefile, then
`snakemake --cluster "bsub -M 4G -R 'select[mem>4G] rusage[mem=4G]' -n2 -G FIXTHIS -q normal -o run_all.stdout -e logs/run_all.stderr.%J" --jobs 22 --latency-wait 10 run_all`

```
### Merging all chromosomes
Sometimes it might be better to process all chromosomes at once (e.g. for small cohorts), which can be done with the following code
```
#!/bin/bash
final_prefix=output/brava_s50.delta.chrALL
chr_vcf_prefix=output/brava_s50.delta

files_to_concat=$final_prefix.files.txt
rm -f $files_to_concat
for c in {1..22}; do 
    echo $chr_vcf_prefix.chr$c.vcf.gz >> $files_to_concat
done
bcftools concat -f $files_to_concat -Oz -o $final_prefix.vcf.gz
plink2 --export bgen-1.2 bits=8 --vcf $final_prefix.vcf.gz dosage=DS --out $final_prefix
# fix the .sample file:
head -2 $final_prefix.sample > $final_prefix.sample1
awk 'NR>2{print 1,$2,0,0}' $final_prefix.sample >> $final_prefix.sample1
mv $final_prefix.sample1 $final_prefix.sample
# index the BGEN file
bgenix -index -g $final_prefix.bgen

```
### Association testing with SAIGE
When all the above steps are complete, you can proceed to association testing with SAIGE. For that you need to set the `model_prefix` (path to step-1 files), `grm_prefix`, `bgen`, `annotation`, and a `samples_file` (with the overlap between step-1 files and the bgen. Then you can invoke SAIGE as follows:
```
if [ ! -f $out_prefix.gz ]; then
    Rscript saige_step2_reskat.R \
        --model $model_prefix \
        --grm $grm_prefix \
        --bgen $bgen \
        --sample $samples_file \
        --annot $annotation \
        --out $out_prefix > $out_prefix.log

    gzip $out_prefix
else
    echo "Sum-stat file already exists, nothing to do here!"
fi
```

### Troubleshooting
* If the snakemake pipeline fails for some chromosomes (e.g. 1-3) but not others, try again after allocating more RAM.
* Note that you need to delete all `*.txt` output files before running `prepare_new_rec_enc`, as some are incrementally updated meaning that the output will be appended in an existing file. In that case, you'll get an error while generating the VCF in the next step.
