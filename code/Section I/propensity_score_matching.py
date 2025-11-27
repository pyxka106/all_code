import pandas as pd

from psmpy import PsmPy

hla = pd.read_csv('/meta_data/top_hlas.csv', sep='\t')
hla_list = hla.HLA.values.tolist()

ms_datasets = ['ms01', 'ms02', 'ms03']
hd_datasets = ['hd01', 'hd02', 'microhd']

new_df = pd.DataFrame()
sample_df = pd.DataFrame()
for ms_dataset in ms_datasets:
    for hd_dataset in hd_datasets:
        
        ms_hla = pd.read_csv(f'/meta_data/{ms_dataset}_hla_typing.csv', sep='\t')
        ms_hla = ms_hla[['subject_id'] + hla_list]
        ms_depth = pd.read_csv(f'/seqs_depth/{ms_dataset}.csv', sep='\t')
        ms_depth = ms_depth[ms_depth['seqs_depth'] > 50000]
        ms_df = ms_hla.merge(ms_depth, on='subject_id')
        ms_df.insert(1, 'MS_status', 'MS')
#         ms_df = ms_df[ms_df['HLA_DRB1*15'] == 1]

        hd_hla = pd.read_csv(f'/meta_data/{hd_dataset}_hla_typing.csv', sep='\t')
        hd_hla = hd_hla[['subject_id'] + hla_list]
        hd_depth = pd.read_csv(f'/seqs_depth/{hd_dataset}.csv', sep='\t')
        hd_depth = hd_depth[hd_depth['seqs_depth'] > 50000]
        hd_df = hd_hla.merge(hd_depth, on='subject_id')
        hd_df.insert(1, 'MS_status', 'HD')
#         hd_df = hd_df[hd_df['HLA_DRB1*15'] == 1]
        
        df = pd.concat([ms_df, hd_df], ignore_index=True, sort=False)
        
        df['MS_status'] = df['MS_status'].map({'MS': 1, 'HD': 0})
        df = df.drop(columns=['seqs_depth'])
        
        psm = PsmPy(df, treatment='MS_status', indx='subject_id', exclude = [])

        # Estimate propensity scores
        psm.logistic_ps(balance=False)

        # Perform nearest neighbor matching
        psm.knn_matched(matcher='propensity_logit', replacement=False, caliper=0.2)

        matched_data = psm.df_matched
        test_df = pd.DataFrame(matched_data.groupby('MS_status').mean(numeric_only=True))
        test_df1 = test_df.T
        test_df1.columns=[0,1]
        test_df1.columns=[hd_dataset, ms_dataset]
        new_df = pd.concat([new_df, test_df1], axis=1)
#         print(ms_dataset, hd_dataset, test_df1)
    
        matched_df = matched_data[matched_data['MS_status'] == 0].drop_duplicates()
        control_df = matched_data[matched_data['MS_status'] == 1].drop_duplicates()
        print(len(matched_df), len(control_df), ms_dataset, hd_dataset)
        matched_data = matched_data[['subject_id']]

        model = ({
            'sample_number': [len(matched_df)]
            # 'group': ['positive', 'negative']
            })

        model_df = pd.DataFrame(model)
        sample_df = pd.concat([sample_df, model_df], ignore_index=True, sort=False)
new_df
