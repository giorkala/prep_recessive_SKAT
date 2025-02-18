#!/bin/python
"""
## Genotype preparation for Recessive SKAT v2 ##

Based on pseudo-variants for biallelic genotypes. Reads a list of biallelic genotypes and makes a VCF file, to use in downstream analyses. 
Could also be used with any list of dosages (though only {0,1} will be considered).

Note: might be slow as it's not optimised

==================
GK - Oct 3rd, 2024
"""

import numpy as np
import pandas as pd
pd.options.mode.copy_on_write = True
import argparse

parser = argparse.ArgumentParser(description="Generate a VCF file for Recessive-SKAT.")
parser.add_argument("--chrom", "-c", help="chromosome index", required=True, type=int)
parser.add_argument("--samples", "-s", help="list of samples to work with", required=True, type=str, default=None)
parser.add_argument("--genotypes", "-g", help="list of genotypes (output from call_chets)", required=False, type=str)
parser.add_argument("--out", "-o", help="prefix for any output", required=True, type=str)
args = parser.parse_args()

CHR=f'chr{args.chrom}'
geno_samples_list = pd.read_csv( args.genotypes, sep='\t', header=None)
geno_samples_list.columns = ['ID', 'CHR', 'Gene', 'GT', 'AC', 'Marker']
geno_samples_list = geno_samples_list[['ID','Marker']]

samples_all = pd.read_csv(args.samples, sep='\t', header=None) #[0].values
samples_all['index'] = range(samples_all.shape[0])
samples_all.set_index(0, inplace=True)
print(f"\nGenerating VCF for {len(samples_all)} samples")

geno_markers_dict = {}
for sampleID, markerID in geno_samples_list.to_numpy():
    if markerID not in geno_markers_dict:
        geno_markers_dict[markerID] = [sampleID]
    else:
        geno_markers_dict[markerID].append(sampleID)

# add a minimal header
header=[
"##fileformat=VCFv4.2",
"##FILTER=<ID=PASS,Description=\"All filters passed\">",
"##contig=<ID="+CHR+">",
"##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage\">",
"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + '\t'.join( samples_all.index )
]

# import gzip
# save to uncompressed format, then compress with bgzip

total_AC = 0
# with gzip.open( args.out, 'wb' ) as fout:
with open( args.out, 'w' ) as fout:
    fout.write('\n'.join(header).encode())
    
    for i, snp in enumerate(geno_markers_dict.keys()):

        genotypes = np.zeros(len(samples_all), dtype=int )
        genotypes[ samples_all.loc[geno_markers_dict[snp]] ] = 1

        # if i==10: 
        #     # print( f"Total AC for {gene}: {np.sum(genotypes)}" )
        #     break
        
        total_AC += np.sum(genotypes)

        # as a proxy for POS we'll use the index of the variant in the list of all variants
        row = '\n' + CHR + '\t' + str(i) + '\t' + snp + "\tA\tB\t.\t.\t.\tDS"

        for gt in genotypes:
            row += '\t'+str(gt)
        
        fout.write(row.encode())

assert total_AC == len(geno_samples_list), "Something went wrong with the genotype annotation."

print("Done, total AC=", total_AC)
# end-of-script