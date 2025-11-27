import pandas as pd
import os

import numpy as np

seq_df = pd.read_csv('/twins_fisher_seqs_0.05.csv', sep='\t')
seq_df = seq_df.iloc[:, 0:3]
seq_df.columns=['v_resolved', 'amino_acid', 'j_resolved']
# v_list1 = [x.replace('*', '-') for x in v_list]
files = [x[2] for x in os.walk('/hips_normalized/')][0]

for file in files: ## file is single HIPs subject
    subject = file.split('.')[0]
    print(subject)

    df = pd.read_csv(f'/Users/jingyun/Documents/hips_normalized/{file}', sep='\t')
    df = df[df['v_gene'].str.contains('N') == False]
    df = df[df['j_gene'].str.contains('N') == False]
    df = df[df['v_gene'].str.contains('OR') == False]
    df = df[df['amino_acid'].str.contains('\*') == False]
    df = df[df['frame_type'].str.contains('In') == True]
    unique_df = df[['v_resolved', 'amino_acid', 'j_resolved']].drop_duplicates()


    all_df = pd.merge(seq_df, unique_df, on=['v_resolved', 'amino_acid', 'j_resolved'], how='left',
                    indicator='exists').drop_duplicates()

    # add column to show if each row in first DataFrame exists in second
    all_df['exists'] = np.where(all_df.exists == 'left_only', False, True)

    # final_df = all_df[all_df['exists'] == True]
    all_df[['exists']] = all_df[['exists']].replace({True: 1, False: 0})

    all_df = all_df.rename(columns={"exists": f'{subject}'})
    seq_df = seq_df.merge(all_df, how='right')
    seq_df.iloc[:, 3:] = seq_df.iloc[:, 3:].fillna(0).astype(int)
    # print(seq_df)
seq_df.to_csv(f'/Users/jingyun/Documents/genetics/hla_ms_specific_seqs_ifn_gla/single_hla_dpb104/hips_presence_table.csv', sep='\t', index=False)
