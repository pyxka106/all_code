import numpy as np
import pandas as pd

## trim amino acids

df = pd.read_csv('/ms_seqs.csv', sep='\t')
df = df[(df['v_resolved'] == 'TRBV19*01') & (df['j_resolved'] == 'TRBJ02-01*01')]

df['cdr3_length'] = df.amino_acid.apply(lambda x: len(x))
df['kmers'] = df.apply(lambda row: ''.join(list(row['amino_acid'])[:-4]), axis=1)
df['kmers'] = df.apply(lambda row: ''.join(list(row['kmers'])[4:]), axis=1)
df = df[df['kmers'] != '']

ms_df = df[df['cdr3_length'] == 15] #.drop_duplicates()
ms_df.insert(5, 'MS_status', 'MS')

df = pd.read_csv('/hd_seqs.csv', sep='\t')
df = df[(df['v_resolved'] == 'TRBV19*01') & (df['j_resolved'] == 'TRBJ02-01*01')]

df['cdr3_length'] = df.amino_acid.apply(lambda x: len(x))
df['kmers'] = df.apply(lambda row: ''.join(list(row['amino_acid'])[:-4]), axis=1)
df['kmers'] = df.apply(lambda row: ''.join(list(row['kmers'])[4:]), axis=1)
df = df[df['kmers'] != '']

hd_df = df[df['cdr3_length'] == 15] #.drop_duplicates()
hd_df.insert(5, 'MS_status', 'HD')

final_df = pd.concat([ms_df, hd_df], ignore_index=True, sort=False)

## atchley factors

AAMetric_Atchley = np.array([
    [-0.59145974, -1.30209266, -0.7330651,  1.5703918, -0.14550842],
    [-1.34267179,  0.46542300, -0.8620345, -1.0200786, -0.25516894],
    [1.05015062,  0.30242411, -3.6559147, -0.2590236, -3.24176791],
    [1.35733226, -1.45275578,  1.4766610,  0.1129444, -0.83715681],
    [-1.00610084, -0.59046634,  1.8909687, -0.3966186,  0.41194139],
    [-0.38387987,  1.65201497,  1.3301017,  1.0449765,  2.06385566],
    [0.33616543, -0.41662780, -1.6733690, -1.4738898, -0.07772917],
    [-1.23936304, -0.54652238,  2.1314349,  0.3931618,  0.81630366],
    [1.83146558, -0.56109831,  0.5332237, -0.2771101,  1.64762794],
    [-1.01895162, -0.98693471, -1.5046185,  1.2658296, -0.91181195],
    [-0.66312569, -1.52353917,  2.2194787, -1.0047207,  1.21181214],
    [0.94535614,  0.82846219,  1.2991286, -0.1688162,  0.93339498],
    [0.18862522,  2.08084151, -1.6283286,  0.4207004, -1.39177378],
    [0.93056541, -0.17926549, -3.0048731, -0.5025910, -1.85303476],
    [1.53754853, -0.05472897,  1.5021086,  0.4403185,  2.89744417],
    [-0.22788299,  1.39869991, -4.7596375,  0.6701745, -2.64747356],
    [-0.03181782,  0.32571153,  2.2134612,  0.9078985,  1.31337035],
    [-1.33661279, -0.27854634, -0.5440132,  1.2419935, -1.26225362],
    [-0.59533918,  0.00907760,  0.6719274, -2.1275244, -0.18358096],
    [0.25999617,  0.82992312,  3.0973596, -0.8380164,  1.51150958]
])

# Labels
amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
atchley_factors = ["pah", "pss", "ms", "cc", "ec"]

# Create mapping dictionary
atchley_df = pd.DataFrame(AAMetric_Atchley, index=amino_acids, columns=atchley_factors)

# Function to convert one peptide string to Atchley vector
def ptuple2atchley(ptuple):
    return np.concatenate([atchley_df.loc[aa].values for aa in ptuple])

# Function to apply Atchley conversion for a list of peptides
def tcr2atchley(seq_list):
    return np.array([ptuple2atchley(seq) for seq in seq_list])

def col_name(seq_list):
    return [f"{factor}{pos+1}" for seq in seq_list[0:1] for pos in range(len(seq)) for factor in atchley_factors]

test_df = final_df.copy()

kmers = test_df['kmers'].tolist()
atchley_matrix = tcr2atchley(kmers)
columns_named = col_name(kmers)

atchley_result = pd.DataFrame(atchley_matrix)
atchley_result.columns=columns_named
atchley_result['MS_status'] = test_df['MS_status'].values.tolist()
atchley_result['amino_acid'] = test_df['amino_acid'].values.tolist()

### plot

import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import mannwhitneyu

from statsmodels.stats.multitest import multipletests

def convert_pvalue_to_asterisks(pvalue):
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
    return 'ns'

fig_df = atchley_result.copy()

fig = plt.figure(figsize=(20,20))
ax1 = plt.subplot2grid(shape=(4,6), loc=(0,0), colspan=2)
ax2 = plt.subplot2grid((4,6), (0,2), colspan=2)
ax3 = plt.subplot2grid((4,6), (0,4), colspan=2)
ax4 = plt.subplot2grid((4,6), (1,1), colspan=2)
ax5 = plt.subplot2grid((4,6), (1,3), colspan=2)
fig.suptitle(f"TRBV19*01 15AA TRBJ02-01*01 Atchley Factors", fontsize=18, y=0.91)
cmap=[sns.color_palette('Set1')[index] for index in [2,0]]

axs = [ax1, ax2, ax3, ax4 ,ax5]
factors_name = ['I Hydrophobicity', 'II Secondary Structure', 'III Size/Mass', 
                'IV Codon Degeneracy', 'V Electric Charge']

for factor, ax, title in zip(atchley_factors, axs, factors_name):
    
    group = [col for col in columns_named if factor in col]
    sub_df = fig_df[['MS_status'] + group].reset_index()
    df1 = sub_df.melt(id_vars=['index', 'MS_status'], var_name="position", value_name="AF_values").set_index('index')
    
    pvalus = []
    for p in group:

        min_df = df1[df1['position'] == p]
        ms_df = min_df[min_df['MS_status'] == 'MS']
        hd_df = min_df[min_df['MS_status'] == 'HD']

        statis, p_value = mannwhitneyu(hd_df['AF_values'], ms_df['AF_values'])
        
        pvalus.append(p_value)
    rejected, pvals_corrected, _, _  = multipletests(pvalus, alpha=0.05, method='bonferroni')
    
    sns.violinplot(x='position', y='AF_values', data=df1, hue='MS_status',
                     legend=True, hue_order=['HD', 'MS'],
                     palette=cmap, ax=ax) # .set_title('I Hydrophobicity')
    i=0
    texts = [t.get_text()  for t in ax.get_xticklabels()]
    for f, padj in zip(texts, pvals_corrected):
        min_df = df1[df1['position'] == p]
        value = min_df.AF_values.max()
        if padj < 0.05:
            star = convert_pvalue_to_asterisks(padj)
            ax.text(i, value, f'{round(padj, 2)}\n{star}',ha="center")
        i += 1
    ax.get_legend().remove()
    ax.set_title(f'{title}')
    
handles, labels = ax1.get_legend_handles_labels()

# Shared legend outside the figure
fig.legend(handles, labels, loc='upper left', ncol=1, bbox_to_anchor=(0.78, 0.68))
plt.show()
plt.close()

## clustering 
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
from numpy import unique

import seaborn as sns

fig_df1 = fig_df.drop_duplicates()
x1 = fig_df1.iloc[:, :-2].to_numpy()
pca = PCA(n_components = 2)
z1 = pca.fit_transform(x1)

dfx1 = pd.DataFrame()
# dfx1[['kmers']] = fig_df[['kmers']]
dfx1["comp-1"] = z1[:, 0]
dfx1["comp-2"] = z1[:, 1]
dfx1['MS_status'] = fig_df1.MS_status.values

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(18, 7))
fig.suptitle("Clustering with Atchley factors -- TRBV19*01 & TRBJ02-01*01 & 15AA", fontsize=18, y=0.92)

sns.scatterplot(x="comp-1", y="comp-2", alpha=0.5, edgecolors='none',data=dfx1, 
                      hue='MS_status', hue_order=['HD', 'MS'], palette=cmap,
                      s=30, legend=True, ax=axs[0])


x1 = dfx1.iloc[:, 0:2].to_numpy() 
model = DBSCAN(eps=0.3, min_samples=15)

yhat = model.fit_predict(x1)
dfx1.insert(3, 'clusters', yhat)

clusters = unique(yhat)

sns.scatterplot(x="comp-1", y="comp-2", alpha=0.5, edgecolors='none',
                      palette=sns.color_palette('tab10'),
                      data=dfx1, hue='clusters', s=30, legend=True, ax=axs[1])


plt.show()
plt.close()
