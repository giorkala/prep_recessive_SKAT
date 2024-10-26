## Development of Recessive SKAT ##
#### Genotype preparation based on pseudo-variants for biallelic genotypes. 

Requirements:
* biallelic genotypes (obtained from `call_chets`)
* variant annotation
* allele frequencies
* ranking of paired consequences, e.g. "pLOF|damaging" > "damaging|damaging"
* collapsing criteria for rare biallelic genotypes
    * any collapsing should follow the extraction of AF products, as we'll then need some sort of average

### Overview of framework
This framework is developed as a Snakemake pipeline to streamline the preparation, with the following tasks:
* Step 0: Detect all biallelic genotypes, obtain variant annotation and allele frequencies
* Step 1: `prepare_new_rec_enc.py`
    1. For each individual, select the genotype with the worst consequence and assign an identifier
    2. Extract the frequency of the corresponding variant(s)
    3. Assign weights and generate annotation files for SAIGE/Regenie
* Step 2: `generate_vcf.py`; this offers a quick way to create a VCF file based on the genotypes obtained previously, one per chromosome.
* Step 3: Concatenate all VCFs to one file, then convert to BGEN for higher efficiency. Repeat for other annotation files.

### How to run
```
git clone --recurse-submodules git@github.com:giorkala/prep_recessive_SKAT.git
# then load module for snakemake
# then edit paths in Snakefile, then
`snakemake --cluster "bsub -M 4G -R 'select[mem>4G] rusage[mem=4G]' -n2 -G FIXTHIS -q normal -o run_all.stdout -e logs/run_all.stderr.%J" --jobs 22 --latency-wait 10 run_all`

```

### Association testing with SAIGE
(work in progress)
