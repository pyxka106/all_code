import pandas as pd
from tcrdist.pgen import OlgaModel
from tcrdist.repertoire import TCRrep

df = pd.read_csv('/twins_hla_typing.csv', sep='\t')
# subject_df = df[df['HLA_DRB1*15'] == 1]
subject_list = df.subject_id.values.tolist()

for subject in subject_list:
    print(subject)
    df = pd.read_csv(f'/PBMC_TWINS/TWIN_{subject}.csv', sep='\t')
    df = df[df['v_gene'].str.contains('N') == False]
    df = df[df['j_gene'].str.contains('N') == False]
    df = df[df['v_gene'].str.contains('OR') == False]
    df = df[df['amino_acid'].str.contains('\*') == False]
    df = df[df['frame_type'].str.contains('In') == True]
    df1 = df.head(50000)
    df1 = df1[['v_resolved', 'amino_acid', 'j_resolved']]
    df1.columns = ['v_b_gene', 'cdr3_b_aa', 'j_b_gene']

    df1['v_b_gene'] = df1['v_b_gene'].str.replace('TRBV0', 'TRBV')
    df1['v_b_gene'] = df1['v_b_gene'].str.replace('-0', '-')
    df1['j_b_gene'] = df1['j_b_gene'].str.replace('TRBJ0', 'TRBJ')
    df1['j_b_gene'] = df1['j_b_gene'].str.replace('-0', '-')

    
    # df1 = df1[df1['v_b_gene'].isin(trbv_list)]
    # print(df1)
    
    tr = TCRrep(cell_df = df1, 
                organism = 'human', 
                chains = ['beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv', 
                store_all_cdr = False)

    olga_beta  = OlgaModel(chain_folder = "human_T_beta", recomb_type="VDJ")

    tr.clone_df['pgen_cdr3_b_aa'] = olga_beta.compute_aa_cdr3_pgens(
        CDR3_seq = tr.clone_df.cdr3_b_aa, V_usage_mask_in=tr.clone_df.v_b_gene, J_usage_mask_in=tr.clone_df.j_b_gene)  #, V_usage_mask_in=tr.clone_df.v_b_gene, J_usage_mask_in=tr.clone_df.j_b_gene

    result_df = tr.clone_df[['v_b_gene', 'cdr3_b_aa', 'j_b_gene', 'pgen_cdr3_b_aa']]

    result_df.to_csv(f'/genp_results/{subject}.csv', sep='\t', index=False)
