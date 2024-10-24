# ranking of consequences for pairs of variants (less is worse)
rank_consq = {
    'pLoF': 1,
    'pLoF|pLoF': 2,
    'pLoF|damaging_missense': 3,
    'damaging_missense|pLoF': 3,
    'pLoF|other_missense': 4,
    'other_missense|pLoF': 4,
    'damaging_missense': 5,
    'damaging_missense|damaging_missense': 6,
    'damaging_missense|other_missense': 7,
    'other_missense|damaging_missense': 7,
    'other_missense': 8,
    'other_missense|other_missense': 9,
    'pLoF|synonymous': 0,
    'synonymous|pLoF': 0,
    'damaging_missense|synonymous': 0,
    'synonymous|damaging_missense': 0,
    'other_missense|synonymous': 0,
    'synonymous|other_missense': 0,
    'synonymous': 0,
    'synonymous|synonymous': 0
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