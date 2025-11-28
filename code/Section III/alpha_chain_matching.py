import os
import pandas as pd

folders = [x[1] for x in os.walk('/large_files/scrna/')][0]

rslt_df = pd.DataFrame()
for folder in folders:
    
    files = [x[2] for x in os.walk(f'/large_files/scrna/{folder}/')][0]
    for file in files:

        df = pd.read_csv(f'/large_files/scrna/{folder}/{file}', sep=',')
        subject_id = file.split('.')[0]

        if df.v_gene.isna().all():
            continue
        
        if df.j_gene.isna().all():
                
            df1 = df[(df['chain'] == 'TRB') & 
                     (df['v_gene'].str.contains('TRBV19')) & 
#                      (df['j_gene'].str.contains('TRBJ2-1')) & 
                     (df['cdr3'].str.contains('DR'))]['barcode'].unique()
            final_df = df[df['barcode'].isin(df1)] # .drop(columns='barcode')
            final_df.insert(7, 'subject_id', subject_id)
            if len(final_df) < 2:
                continue

            rslt_df = pd.concat([rslt_df, final_df], ignore_index=True, sort=False)
            rslt_df = rslt_df[['barcode', 'chain', 'v_gene', 'cdr3', 'j_gene', 'subject_id']]
            
        else: 
            df1 = df[(df['chain'] == 'TRB') & 
                     (df['v_gene'].str.contains('TRBV19')) & 
                     (df['j_gene'].str.contains('TRBJ2-1')) & 
                     (df['cdr3'].str.contains('DR'))]['barcode'].unique()
            final_df = df[df['barcode'].isin(df1)] # .drop(columns='barcode')
            final_df.insert(7, 'subject_id', subject_id)
            if len(final_df) < 2:
                continue

            rslt_df = pd.concat([rslt_df, final_df], ignore_index=True, sort=False)
            rslt_df = rslt_df[['barcode', 'chain', 'v_gene', 'cdr3', 'j_gene', 'subject_id']]
#         print(rslt_df)
    alpha = rslt_df[rslt_df['chain'] == 'TRA'].reset_index(drop=True)
    beta  = rslt_df[rslt_df['chain'] == 'TRB'].reset_index(drop=True)

    # Rename columns for each
    alpha = alpha.rename(columns={
        'barcode': 'barcode',
        'cdr3': 'cdr3a',
        'v_gene': 'TRAV',
        'j_gene': 'TRAJ'
    }).drop(columns=['chain'])

    beta = beta.rename(columns={
        'barcode': 'barcode',
        'cdr3': 'cdr3b',
        'v_gene': 'TRBV',
        'j_gene': 'TRBJ'
    }).drop(columns=['chain'])

    result_df = alpha.merge(beta, on=['barcode', 'subject_id'])

#     result_df = pd.concat([merged, rslt_df1], ignore_index=True, sort=False)

cols_at_end = ['subject_id']
result_df = result_df[[c for c in result_df if c not in cols_at_end] 
        + [c for c in cols_at_end if c in result_df]]
result_df = result_df.dropna()
