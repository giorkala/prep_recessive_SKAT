import numpy as np
# from scipy.stats import beta
    
def dummy_split(x):
    if type(x)==tuple:
        return ';'.join(x)
    else:
        return x
    
def write_saige_file(df, out_file, col_snp='SNP', col_gene='Gene', 
                     col_consq='Consq', col_weight='Weight'):
    # write the SAIGE group file to disk
    print(f"Saving the annotation to {out_file}")
    with open(out_file, 'w') as fout:
        for gene_id, variants in df.groupby(col_gene):
            annotated = ~variants[col_consq].isna()
            
            snps = ' '.join(
                [snp for snp in variants[col_snp][annotated].values]
            )
            annos = ' '.join(
                [anno for anno in variants[col_consq][annotated].values]
            )
            weights = ' '.join(
                [f'{w:.4f}' for w in variants[col_weight][annotated].values]
            )
            
            fout.write(f"{gene_id} var {snps}\n")
            fout.write(f"{gene_id} anno {annos}\n")
            fout.write(f"{gene_id} weight {weights}\n")


# ranking of consequences for pairs of variants (less is worse)
rank_consq = {
    'pLoF': 1,
    'pLoF|pLoF': 2,
    'pLoF|damaging_missense': 3,
    'damaging_missense|pLoF': 3,
    'damaging_missense': 4,
    'damaging_missense|damaging_missense': 4,
    'pLoF|other_missense': 5,
    'other_missense|pLoF': 5,
    'damaging_missense|other_missense': 7,
    'other_missense|damaging_missense': 7,
    'other_missense': 8,
    'other_missense|other_missense': 9,
    'pLoF|synonymous': 10,
    'synonymous|pLoF': 10,
    'damaging_missense|synonymous': 10,
    'synonymous|damaging_missense': 10,
    'other_missense|synonymous': 10,
    'synonymous|other_missense': 10,
    'synonymous': 10,
    'synonymous|synonymous': 10
 }

def get_worst_consequence(x, gene_annot):
    """
    this should go through all possible combinations of genotypes and return the worst consequence.
    """
    hap0 = x.split('|')[0]
    hap1 = x.split('|')[1]
    all_cases = {}
    # print('test')
    for A in hap0.split(';'):
        for B in hap1.split(';'):
            if A == B:
                all_cases[A] = gene_annot[A]
            else:
                # sort wrt POS
                if int(A.split(':')[1]) < int(B.split(':')[1]):
                    all_cases[A,B] = gene_annot[A] +'|'+ gene_annot[B]
                else:
                    all_cases[B,A] = gene_annot[B] +'|'+ gene_annot[A]

    # now select the worst one; k is variant(s), v is consequence(s)
    all_cases = {rank_consq[v]:(k,v) for k, v in all_cases.items()}
    m = min(all_cases.keys())
    return all_cases[m] 

def get_freq_product( markers, df_maf ):
    if len(markers) == 2:
        a = df_maf.loc[markers[0]].MAF
        b = df_maf.loc[markers[1]].MAF
        return 2 * a * b
    else:
        return df_maf.loc[markers].MAF ** 2

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
    'other_missense|other_missense': 'other_missense|other_missense',
    'other_missense': 'other_missense',
    'other_missense|other_missense': 'other_missense|other_missense',
    'pLoF|synonymous': 'NA',
    'pLoF|synonymous': 'NA',
    'damaging_missense|synonymous': 'NA',
    'synonymous|damaging_missense': 'NA',
    'other_missense|synonymous': 'NA',
    'synonymous|other_missense': 'NA',
    'synonymous': 'synonymous',
    'synonymous|synonymous': 'synonymous'
}

new_weights = {
    'synonymous|synonymous': 1,
    'synonymous': 1,
    'other_missense': 2,
    'other_missense|other_missense': 2,
    'other_missense|damaging_missense': 2,
    'damaging_missense|other_missense': 2,
    'other_missense|pLoF': 3,
    'pLoF|other_missense': 3,
    'damaging_missense': 4,
    'damaging_missense|damaging_missense': 4,
    'pLoF|damaging_missense': 5,
    'damaging_missense|pLoF': 5,
    'pLoF|pLoF': 5,
    'pLoF': 5,
}

def annot_relabelling(x):
    return new_labels[x]

# weighting scemes - suggested: set_weights
# to be used like df['weights'] = df.worst_consq.apply(lambda x: set_weights(x))
def set_class_weights(x, new_weights=new_weights):
    return new_weights[x]
    
### withdrawn/failed attempts ###
# def set_uniform_weights(x):
#     return beta.pdf(1, 1, 25)

# def set_af_weights(x):
#     return beta.pdf(x, 1, 25)

# def set_random_weights(x):
#     return beta.pdf(np.random.rand()/20, 1, 25)
