#!/bin/python
"""
## Genotype preparation for Recessive SKAT v2 ##
based on pseudo-variants for biallelic genotypes

For each gene,
* detect all biallelic genotypes/carriers
    * For each individual, select the genotype with the worst consequence and assign an identifier
    * Extract the frequency of the corresponding variants
* then make a variant-based index of the genotypes
* generate SAIGE or Regenie annotation files

Note: Regenie doesn't support custom weights, thus its better to proceed with SAIGE
 
###################
GK - Oct 1st, 2024
"""

import numpy as np
import pandas as pd
pd.options.mode.copy_on_write = True
import argparse
from utils_new_rec_enc import get_worst_consequence, get_freq_product, write_saige_file, set_class_weights

parser = argparse.ArgumentParser(description="Genotype preparation for Recessive SKAT v2.")
parser.add_argument("--anno", "-a", help="file with variant-gene annotation (as in Regenie)", required=True, type=str)
parser.add_argument("--chrom", "-c", help="chromosome index", required=True, type=int)
parser.add_argument("--maf", "-m", help="file to MAFs", required=True, type=str, default=None)
parser.add_argument("--maf_header", help="header available in the MAF file, e.g. if from PLINK", required=True, action='store_true', default=False)
parser.add_argument("--genotypes", "-g", help="list of genotypes (output from call_chets)", required=False, type=str)
parser.add_argument("--tag", "-t", help="marker identifier", required=True, type=str, default='nonsyn')
parser.add_argument("--out", "-o", help="prefix for any output", required=True, type=str)
args = parser.parse_args()

output_prefix=args.out
C=args.chrom
# the next will be used in the pseudo-variant identifiers
if args.tag == 'synonymous':
    snp_tag='syn'
else:
    snp_tag=''

### 1. Load and prepare input ###
df_annot = pd.read_csv( args.anno, sep='\s+', names=['SNP','Gene','Consq'])
# simplify terms
df_annot.replace({'other_missense_or_protein_altering':'other_missense',
                  'damaging_missense_or_protein_altering':'damaging_missense'}, inplace=True)

if args.maf_header:
   df_maf = pd.read_csv( args.maf, sep='\t' ).set_index('ID')
   assert 'MAF' in df_maf.columns, 'Mismatch with the MAF file header!'
   df_maf['MAF'] = df_maf.ALT_FREQS.values
else:
   df_maf = pd.read_csv( args.maf, sep='\t', header=None).set_index(0)
   # assuming MAF is the fourth column...
   df_maf['MAF'] = df_maf[3] / df_maf[4] / 2

df = pd.read_csv( args.genotypes, sep='\t', header=None)
df.columns= ['ID','CHR','Gene','GT','AC','Variants']

### 2. make a gene index ###
# POS might not always work due to overlapping variants. I'll first use POS to rank genes, then assign a unique index to each gene.
df_annot['POS'] = [int(x.split(':')[1]) for x in df_annot['SNP']]
gene_index = df_annot[['Gene','POS']].groupby('Gene').agg('min').reset_index()#.set_index('Gene')
gene_index = gene_index.set_index('Gene').sort_values(by='POS')
gene_index['geneID'] = [f'chr{C}:{i}' for i in range(1,gene_index.shape[0]+1)]

print(f'\nWill work for chrom-{C}-{args.tag} with {gene_index.shape[0]} genes and {df_annot.shape[0]} (total) variants.')

## check if any two genes have overlapping POS
if gene_index['geneID'].nunique() != gene_index.shape[0]:
    print('WARNING: overlapping gene IDs!')

### 3. go through all genotypes and assign a marker ID ###
# geno_samples_list = pd.DataFrame()
df_markers = pd.DataFrame(columns=['gtype_ID', 'worst_gtype', 'worst_consq', 'count', 'freq'])
geno_total = 0
duplicated = []
for gene in gene_index.index:

    df_tmp = df.query('Gene == @gene')
    df_tmp.reset_index(inplace=True, drop=True)

    gene_annot = df_annot.query('Gene==@gene').set_index('SNP')['Consq']
    # drop duplicated variants - TODO: check for potential issues
    duplicated.append(sum(gene_annot.index.duplicated()))
    gene_annot = gene_annot[~gene_annot.index.duplicated()]

    # get the worst consequence for each genotype and all combinations
    things = df_tmp.Variants.apply(lambda x: get_worst_consequence(x, gene_annot))
    df_tmp['worst_gtype'] = [x[0] for x in things]
    df_tmp['worst_consq'] = [x[1] for x in things]

    # make a table with info about each genotype
    df_sum = df_tmp[['worst_gtype','worst_consq']].value_counts().reset_index()
    df_sum.columns = ['worst_gtype','worst_consq','count']
    df_sum['gtype_ID'] = [f'{gene_index.loc[gene].geneID}:{snp_tag}{i+1}' for i in range(df_sum.shape[0])]
    df_sum['freq'] = df_sum.worst_gtype.apply(lambda x: get_freq_product(x, df_maf))
    df_sum['Gene'] = gene
    # save the marker info
    df_markers = pd.concat([df_markers, df_sum[['gtype_ID', 'worst_consq', 'Gene']] ])
    df_sum[['gtype_ID', 'worst_gtype', 'worst_consq', 'count', 'freq']].to_csv( output_prefix+'.marker_info.txt', 
                                                                            mode='a', header=False, sep='\t', index=False)

    # update the original list of genotypes and save to disk
    df_tmp['Marker'] = df_sum.set_index('worst_gtype').loc[df_tmp.worst_gtype, 'gtype_ID'].values
    df_tmp[['ID', 'CHR', 'Gene', 'GT', 'AC', 'Marker']].to_csv( output_prefix+'.txt', mode='a', 
                                                            header=False, sep='\t', index=False)
    
    geno_total += df_tmp.shape[0]

print("Done with genotype annotation. Genotypes loaded:", geno_total)
assert geno_total == df.shape[0], "Something went wrong with the genotype annotation."

### 4. make annotation files for Regenie ###

new_labels = {
    'pLoF': 'pLoF',
    'pLoF|pLoF': 'pLoF|pLoF',
    'pLoF|damaging_missense': 'pLoF|damaging_missense',
    'damaging_missense|pLoF': 'pLoF|damaging_missense',
    'pLoF|other_missense': 'pLoF|other_missense',
    'other_missense|pLoF': 'pLoF|other_missense',
    'damaging_missense': 'damaging_missense',
    'damaging_missense|damaging_missense': 'damaging_missense|damaging_missense',
    'damaging_missense|other_missense': 'damaging_missense|other_missense',
    'other_missense|damaging_missense': 'damaging_missense|other_missense',
    'other_missense': 'other_missense',
    'other_missense|other_missense': 'other_missense|other_missense',
    'pLoF|synonymous': 'NA',
    'synonymous|pLoF': 'NA',
    'damaging_missense|synonymous': 'NA',
    'synonymous|damaging_missense': 'NA',
    'other_missense|synonymous': 'NA',
    'synonymous|other_missense': 'NA',
    'synonymous': 'synonymous',
    'synonymous|synonymous': 'synonymous'
}
def annot_relabelling(x):
    return new_labels[x]

df_markers['consq'] = df_markers.worst_consq.apply(lambda x: annot_relabelling(x))
df_markers[['gtype_ID','Gene','consq']].to_csv( output_prefix+'.regenie_annotation.txt', sep='\t', index=False, mode='a', header=None)

df_markers['Weight'] = df_markers.worst_consq.apply(lambda x: set_class_weights(x))
df_markers = df_markers[['gtype_ID', 'Gene', 'consq', 'Weight']]
write_saige_file(df_markers, output_prefix + '.groupFile.txt', col_consq='consq', col_snp='gtype_ID')

# POS_all = []
# weird = []
# with open(f'{output_prefix}.set_list.txt', 'w') as fout:
#     for gene in df_markers.Gene.unique():
#         snps = df_markers.query('Gene == @gene').gtype_ID.values
#         chrom = snps[0].split(':')[0]
#         POS = min([int(x.split(':')[1]) for x in snps])
#         if POS in POS_all:
#             POS += 1
#             POS_all.append(POS) 
#             weird.append(gene)

#         fout.write(f'{gene}\t{chrom}\t{POS}\t' + ','.join(snps) + '\n')

# if len(weird) > 0:
#     print("WARNING: the following genes have overlapping POS:", weird)

print("All done! Proceed to the next step.")
    
