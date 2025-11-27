from scipy.stats import fisher_exact
import numpy as np
import os
import glob
import pandas as pd


def func(a, b, c, d):
    contingency_table = [[a, b], [c, d]]
    odds_ratio, p_value = fisher_exact(contingency_table, alternative='greater')
    return p_value


os.chdir('/presence_tables/IFN_GLA_untreated/')
files = glob.glob('*.csv')

top_hla = pd.read_csv('/top_hlas.csv', sep='\t')
hla_list = top_hla.HLA.values.tolist()

therapy_df = pd.read_csv('/meta_data/twins_therapy.csv', sep='\t')
hla_df = pd.read_csv('/meta_data/all_hla_typing.csv', sep='\t')

subject_df = therapy_df[(therapy_df['therapy'] == 'none') | (therapy_df['therapy'] == 'IFN') | (therapy_df['therapy'] == 'GLA')]
subject_list = subject_df.subject_id.values.tolist()
hla_df = hla_df[hla_df['subject_id'].isin(subject_list)]

for file in files:
    print(file) ## file id single v gene segment
    df1 = pd.read_csv(file, sep='\t')

    for hla in hla_list:
        print(hla)
        hla1 = hla.replace('*', '')
        sub_df = hla_df[hla_df[hla] == 1]
        sub = sub_df.subject_id.values.tolist()

        df = df1[['v_resolved', 'amino_acid', 'j_resolved'] + sub]

        neg_df = df.loc[:, df.columns.str.contains('HD')]
        pos_df = df.loc[:, df.columns.str.contains('MS')]
        neg_list = neg_df.columns.tolist()
        pos_list = pos_df.columns.tolist()
        test_df = pd.DataFrame()
        test_df[['v_resolved', 'amino_acid', 'j_resolved']] = df[['v_resolved', 'amino_acid', 'j_resolved']]

        # present or absent for TCR
        test_df['MS_with'] = df[pos_list].sum(axis=1)
        test_df['HD_with'] = df[neg_list].sum(axis=1)
        test_df['MS_miss'] = len(pos_list) - df[pos_list].sum(axis=1)
        test_df['HD_miss'] = len(neg_list) - df[neg_list].sum(axis=1)
        test_df.iloc[:, 3:] = test_df.iloc[:, 3:].astype(np.int64)

        test_df['p_value'] = test_df.apply(lambda row: func(row['MS_with'],
                                                            row['HD_with'],
                                                            row['MS_miss'],
                                                            row['HD_miss']),
                                        axis=1)
        result_df = test_df[test_df['p_value'] < 0.05]
        # print(result_df)
        result_df.to_csv(f'/hla_ms_specific_seqs/{hla1}/{file}', sep='\t', index=False)
