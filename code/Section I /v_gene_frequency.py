mport pandas as pd

subject_df = pd.read_csv('/meta_data/twins_hla_typing.csv', sep='\t')
subject_list = subject_df.subject_id.values.tolist()

final_df = pd.DataFrame()
final_df1 = pd.DataFrame()
for subject in subject_list:
    df = pd.read_csv(f'/PBMC_TWINS/TWIN_{subject}.csv', sep='\t')
    df1 = df[~df['v_resolved'].str.contains('OR', na=False)]
    df1 = df1[~df1['v_resolved'].str.contains('A', na=False)]
    df1 = df1[df1['v_resolved'].notna()]
    df1 = df1[df1['amino_acid'].str.contains('\*') == False]
    df1 = df1[df1['frame_type'].str.contains('In') == True]
    unique_df = df1[['v_resolved', 'amino_acid', 'j_resolved']] # .drop_duplicates()
    unique_df['cdr3_length'] = unique_df.amino_acid.apply(lambda x: len(x))

    new_df = pd.DataFrame(unique_df.groupby('v_resolved')['cdr3_length'].mean()).reset_index()
    new_df.insert(0, 'subject_id', subject)
    
    final_df = pd.concat([final_df, new_df], ignore_index=True, sort=False)
    
    new_df1 = pd.DataFrame(unique_df['v_resolved'].value_counts(normalize=True)).reset_index()
    new_df1.columns=['v_gene', 'freq']
    new_df1.insert(0, 'subject_id', subject)
    final_df1 = pd.concat([final_df1, new_df1], ignore_index=True, sort=False)
