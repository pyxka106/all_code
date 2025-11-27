import pandas as pd
import itertools
import glob

import numpy as np

df = pd.read_csv('./user_output/targets.tsv', sep='\t')
pae_files = glob.glob('/TCRdock_mbp/*.npy')

df = df.iloc[: len(pae_files), :]

final_dfl = []
for counter, targetl in df.iterrows():
    query_chainseq = targetl.target_chainseq
    outl = targetl.copy()
    pae_file = pae_files[counter]
    paes = np.load(pae_file)
    
    cs = query_chainseq.split('/')
    # print(cs)
    chain_stops = list(itertools.accumulate(len(x) for x in cs))
    chain_starts = [0]+chain_stops[:-1]
    nres = chain_stops[-1]
    # print(nres)
    outl['model2_ptm_ft4_pae'] = np.mean(paes[:nres,:nres])
    for chain1,(start1,stop1) in enumerate(zip(chain_starts, chain_stops)):

        for chain2 in range(len(cs)):
            start2, stop2 = chain_starts[chain2], chain_stops[chain2]
            pae = np.mean(paes[start1:stop1,start2:stop2])
            outl[f'model2_ptm_ft4_pae_{chain1}_{chain2}'] = pae
    final_dfl.append(outl)

rslt_df = pd.DataFrame(final_dfl)

inter_paes = []
for _, l in rslt_df.iterrows():
    cs = l.target_chainseq.split('/')
    num_chains = len(cs)

    pmhc_chains = range(num_chains-2)
    tcr_chains = range(num_chains-2, num_chains)

    inter_pae = 0.
    for i in pmhc_chains:
        nres_i = len(cs[i])
        for j in tcr_chains:
            nres_j = len(cs[j])
            pae_ij = l[f'model2_ptm_ft4_pae_{i}_{j}']
            pae_ji = l[f'model2_ptm_ft4_pae_{j}_{i}']
            inter_pae += nres_i * nres_j * (pae_ij + pae_ji)
    nres_pmhc = sum(len(cs[x]) for x in pmhc_chains)
    nres_tcr = sum(len(cs[x]) for x in tcr_chains)
    inter_pae /= 2*nres_pmhc*nres_tcr
    inter_paes.append(inter_pae)

rslt_df['pmhc_tcr_pae'] = inter_paes
rslt_df1 = rslt_df[['organism', 'mhc_class', 'mhc', 'peptide', 'va', 'ja', 'cdr3a', 'vb', 'jb', 'cdr3b', 'pmhc_tcr_pae']]
rslt_df1
