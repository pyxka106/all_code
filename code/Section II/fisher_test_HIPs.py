import pandas as pd
import os
from scipy.stats import fisher_exact

def myfunc(a,b,c,d):
    odds_ratio, p_value = fisher_exact([[a,b], [c,d]], alternative='greater')
    return p_value

hla_df = pd.read_csv('/top_hlas.csv', sep='\t')
hla_list = hla_df.HLA.values.tolist()

hla_pheno = pd.read_csv('/meta_data/hips_hla_typing.csv', sep='\t')

seq_df = pd.read_csv('/hips_presence_table.csv', sep='\t')
for hla in hla_list:
    print(hla)
    hla_file = hla.replace('*', '')
    #os.mkdir(f'/hla_specific_seqs/{hla_file}')

    allele_df = hla_pheno[['subject_id', f'{hla}']]
    hla_pos = allele_df[allele_df[f'{hla}'] == 1]
    hla_pos_list = hla_pos.subject_id.values.tolist()
    hla_neg = allele_df[allele_df[f'{hla}'] == 0]
    hla_neg_list = hla_neg.subject_id.values.tolist()
    
    fisher_df = pd.DataFrame()
    fisher_df[['v_gene', 'CDR3', 'j_gene']] = seq_df[['v_resolved', 'amino_acid', 'j_resolved']]
    fisher_df['hla_pos_with'] = seq_df[hla_pos_list].sum(axis=1)
    fisher_df['hla_pos_miss'] = fisher_df.hla_pos_with.apply(lambda x: len(hla_pos_list) - x)

    fisher_df['hla_neg_with'] = seq_df[hla_neg_list].sum(axis=1)
    fisher_df['hla_neg_miss'] = fisher_df.hla_neg_with.apply(lambda x: len(hla_neg_list) - x)
    
    fisher_df['p_value'] = fisher_df.apply(lambda x: myfunc(x.hla_pos_with,  x.hla_pos_miss, x.hla_neg_with, x.hla_neg_miss), axis=1)
    fisher_df = fisher_df[fisher_df['p_value'] < 0.1]
    # print(fisher_df)

    fisher_df.to_csv(f'/{hla_file}.csv', sep='\t', index=False)
