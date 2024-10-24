"""
Snakemake pipeline for preparing genotypes to enable Recessive SKAT analysis in SAIGE

#### Required input
1. BCF files with phased genotypes
2. annotation (SAIGE-like, output from https://github.com/BRaVa-genetics/variant-annotation)
3. variant consequences to work with (e.g. pLOF + damaging_missense + other_missense)
4. Parameters: MIN_AC, max_AF, WORK_DIR
5. a `samples.bcf.txt` file with the list of samples in the BCF files
Check `README.md` for more details.

#### Overview
1. [per chromosome] Use `call_chets` to generate a list of biallelic individuals, separately for syn and nonsyn variants.
2. [per chromosome] Encode all biallelic gentypes to the new format, after merging the two lists. Also make other necessary files for association testing.
3. [once for all] Concatenate all VCFs to one file, then convert to BGEN for higher efficiency. Repeat for other annotation files.

#### Notes
* We'll processs both nonsynonymous and synonymous genotypes, but you can use choose to skip the latter from testing.
* Please fix paths for `input.maf_file` and `input.annot` in rule get_recessive_encoding.

# for LSF:
`snakemake --cluster "bsub -M 4G -R 'select[mem>4G] rusage[mem=4G]' -n2 -G FIXTHIS -q normal -o run_all.stdout -e logs/run_all.stderr.%J" --jobs 50 --latency-wait 10 run_all`

GK - Oct 24th, 2024
"""

## parameters ##
min_AC=5 # minimum number of carriers per gene
max_AF=0.05 # maximum allele frequency per variant
TAG="brava_s50" # tag to identify the output files
cohort="GNH_39k_QC" # cohort identifier

## General file paths and prefices ##
parent_dir="/lustre/scratch123/hgi/mdt2/projects/gnh_industry/Genes_and_Health_2023_02_44k"
WORK_DIR="/FIXTHIS/PATH"
bcf_prefix="FIXTHIS/PATH/phased_genotypes/biobank_tag.notrios.phased_all"
WORK_DIR=f"{parent_dir}/phasing/dev_new_rec"
bcf_prefix=f"{parent_dir}/phasing/phased_genotypes/GNH_39k_QC.notrios.phased_all"
annot_saige_prefix=f"{parent_dir}/brava_vep_annotation/vep105_loftee/out/gnh_39k.brava_s50_nopopmax"
# this is the output from BRaVa's variant annotation; should end in ".chr{C}.saige.txt"
maf_prefix=f"{WORK_DIR}/input/GNH_39k_QC.notrios"
# this was obtained during the phasing process; should end in ".chr{C}.snpinfo"

final_prefix = f"{WORK_DIR}/output/{TAG}.reskat"

import os
for d in ["sandbox", "output", "logs"]:
    os.makedirs(f"{WORK_DIR}/{d}", exist_ok=True)

rule extract_genotypes:
    input:
        bcf = bcf_prefix+".chr{chrom}.bcf"
    output:
        WORK_DIR+"/sandbox/chr{chrom}.PP90af05.txt.gz"
    params:
        max_AF = max_AF
    shell:
        """
        source /etc/profile.d/modules.sh
        module load HGI/softpack/users/bh18/bcftools/1.0
        echo "Generating a file with all genotypes. This can be slow..."
        # select all heterozygous genotypes with PP>0.90 or missing (which implies a common variant)
        bcftools view {input.bcf} --max-af {params.max_AF} | bcftools query -i '(GT="1|1" | GT="0|1" | GT="1|0")' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT %PP\n]' | awk '(length($2)<50)' | awk '($4=="." || $4>0.90)' | gzip > {output}
        """

rule generate_biallelic_genotypes:
    params:
        prep_genemap = "GNH_BRaVa_dev_misc/recessive_encoding/prepare_genemap.py",
        call_chets = "/software/team281/gk18/call_chets/call_chets",
        wd = WORK_DIR,
        tag = TAG
    resources:
        mem_mb=4000,
        threads=1
    input:
        genotypes = rules.extract_genotypes.output,
        annot = annot_saige_prefix+".chr{chrom}.saige.txt",
        samples = f"{WORK_DIR}/samples.bcf.txt"
    output:
        genotypes = WORK_DIR+"/sandbox/biallelic.{TAG}.merged.chr{chrom}.txt",
        var_annot = WORK_DIR+"/sandbox/annotation.{TAG}.chr{chrom}.txt",
    shell:
        """
        # work for nonsynonymous variants
        python {params.prep_genemap} -a {input.annot} \
            -c pLoF damaging_missense_or_protein_altering other_missense_or_protein_altering \
            -o {output.var_annot} 
        {params.call_chets} -g {input.genotypes} -m {output.var_annot} --show-variants | grep 'chet\|hom' > {output}

        # repeat for synonymous variants, adding to the existing file
        annot_syn="{params.wd}/sandbox/annot.chr{wildcards.chrom}.syn.txt"
        python {params.prep_genemap} -a {input.annot} -c synonymous -o $annot_syn
        {params.call_chets} -g {input.genotypes} -m $annot_syn --show-variants | grep 'chet\|hom' >> {output}
        cat $annot_syn >> {output.var_annot}
        rm $annot_syn
        """

rule get_recessive_encoding:
    # Encode all biallelic gentypes to the new format, then make the corresponding annotation files
    input:
        genotypes = rules.generate_biallelic_genotypes.output.genotypes,
        annot = rules.generate_biallelic_genotypes.output.var_annot,
        maf_file = maf_prefix+".chr{chrom}.snpinfo"
    output:
        genotypes = WORK_DIR+"/output/{TAG}.chr{chrom}.txt",
        annot = WORK_DIR+"/output/{TAG}.chr{chrom}.annotation.txt",
        set_list = WORK_DIR+"/output/{TAG}.chr{chrom}.set_list.txt",
        markers = WORK_DIR+"/output/{TAG}.chr{chrom}.marker.info"
    params:
        C = lambda wildcards: wildcards.chrom,
        out_prefix = WORK_DIR+"/output/"+TAG+".chr{chrom}"
    shell:
        """
        # first merge the two lists with annotations from the previous step
        python prepare_new_rec_enc.py -c {params.C} -t merged -g {input.genotypes} -m {input.maf_file} -a {input.annot} -o {params.out_prefix}
        """

rule generate_vcf:
    input:
        rules.get_recessive_encoding.output.genotypes
    output:
        WORK_DIR+"/output/{TAG}.chr{chrom}.vcf.gz"
    params:
        samples = f"{WORK_DIR}/samples.bcf.txt",
        C = lambda wildcards: wildcards.chrom
    shell:
        """
        python generate_vcf.py -s {params.samples} -c {params.C} -g {input} -o {output}
        """

rule run_all:
    input:
        expand(WORK_DIR+"/sandbox/chr{chrom}.PP90af05.txt.gz", chrom=range(1,23), tag={TAG}),
        expand(WORK_DIR+"/sandbox/biallelic.{tag}.merged.chr{chrom}.txt", chrom=range(1,23), tag={TAG}),
        expand(WORK_DIR+"/output/{tag}.chr{chrom}.txt", chrom=range(20,23), tag={TAG}),
        expand(WORK_DIR+"/output/{tag}.chr{chrom}.vcf.gz", chrom=range(22,23), tag={TAG})

rule last_steps:
    input:
        expand(WORK_DIR+"/output/{tag}.chr{chrom}.vcf.gz", chrom=range(21,23), tag={TAG})
    output:
        bgen = final_prefix+".bgen",
        vcf = final_prefix+".vcf.gz",
        sample = final_prefix+".sample"
    params:
        prefix = final_prefix
    shell:
        """
        source /etc/profile.d/modules.sh
        module load HGI/softpack/users/bh18/bcftools/1.0

        files_to_concat={params.prefix}.files.txt
        rm -f $files_to_concat
        for c in $(seq 21 22); do echo {params.prefix}.chr$c.vcf.gz >> $files_to_concat; done
        bcftools concat -f $files_to_concat -Oz -o {output.vcf}

        plink2 --export bgen-1.2 bits=8 --vcf {output.vcf} dosage=DS --out {params.prefix}
        # fix the .sample file
        head -2 {output.sample} > {params.prefix}.sample1
        awk 'NR>2{{print 1,$2,0,0}}' {params.prefix}.sample >> {params.prefix}.sample1
        mv {params.prefix}.sample1 {output.sample}
        """