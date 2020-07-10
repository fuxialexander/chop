#!/usr/bin/env python
# %%
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600
# plt.rcParams['figure.figsize'] = 3, 3
# %%
exprs = pd.read_csv("exprs.ann.txt", sep="\t")
# %%
exprs_filt = exprs.loc[np.sum(exprs.iloc[:, 4:].values > np.log2(300), 1
                              ) > 0, :].iloc[:, 2:].dropna()
exprs_filt.columns = exprs_filt.columns.str.replace('.CEL.gz', '').str.replace(
    'GSM[0-9]+_[0-9]+_', '')
# %%
replicates_map = dict(
    zip(exprs_filt.columns[2:].values,
        exprs_filt.columns[2:].str.replace('_rep.', '').values))

# %%
exprs_repl_mean = pd.concat(
    (exprs_filt.iloc[:, 0:1], exprs_filt.iloc[:, 2:].groupby(
        by=replicates_map, axis=1).mean()), 1).set_index('SYMBOL')

exprs_repl_mean.index.names = ['Genes']
exprs_repl_mean = exprs_repl_mean.groupby('Genes').max()
# %%
new_mul_index = pd.MultiIndex.from_tuples(
    [tuple(i.split('_'))[::-1] for i in exprs_repl_mean.columns.values])

exprs_repl_mean.columns = new_mul_index
# %%
col_order = [
    ('WT', 'PZ'),
    ('WT', 'pHZ'),
    ('WT', 'Upper'),
    ('WT', 'Lower'),
    ('13del', 'PZ'),
    ('13del', 'pHZ'),
    ('13del', 'Upper'),
    ('13del', 'Middle'),
    ('13del', 'Lower'),
    ('13delChopKO', 'PZ'),
    ('13delChopKO', 'PHZ'),
    ('13delChopKO', 'Upper'),
    ('13delChopKO', 'Middle'),
    ('13delChopKO', 'Lower')
]

# %%
exprs_repl_mean = exprs_repl_mean.reindex(col_order, axis=1)
exprs_repl_mean.columns.names = ['Treatment', 'Segment']
# %%
# exprs_repl_mean.to_csv("exprs.formated.csv")

# %%

### Call significant DEG
from scipy.stats import ranksums

p_val_13del = np.array([ranksums(i[1].values[0:4], i[1].values[4:9])[1] for i in exprs_repl_mean.iloc[:, 0:9].iterrows()])

p_val_chop = np.array([ranksums(i[1].values[0:5], i[1].values[5:10])[1] for i in exprs_repl_mean.iloc[:, 4:].iterrows()])

p_val_WT_chop = np.array([ranksums(i[1].values[0:4], i[1].values[4:9])[1] for i in exprs_repl_mean.iloc[:, [0,1,2,3,9,10,11,12,13]].iterrows()])

# exprs_sig = exprs_repl_mean[(p_val_13del<0.05) | (p_val_chop<0.05) | (p_val_13_delchop<0.05)]

wt_13del = exprs_repl_mean[(p_val_13del<0.05)].copy()
exprs_repl_mean['P_WT_13del'] = p_val_13del
exprs_repl_mean['P_13del_chop'] = p_val_chop
exprs_repl_mean['P_WT_chop'] = p_val_WT_chop
exprs_sig = exprs_repl_mean[(p_val_13del<0.05)]

exprs_sig.to_csv("deg.wt_13del.csv")

# exprs_sig
# %%
import seaborn as sns

# %%
# sns.clustermap(exprs_sig.iloc[:,0:9], cmap="vlag", col_cluster=False, z_score=0, method='ward', metric='euclidean')

# %%
upr = np.where(wt_13del.iloc[:,0:4].mean(1).values < wt_13del.iloc[:,4:9].mean(1).values)[0]
# print(upr)
wt_13del.iloc[upr, :].to_csv("deg.wt_13del.13del_upregulated.csv")
genes = wt_13del.iloc[upr, :]
chop_invariant = exprs_repl_mean[(p_val_13del<0.05) & (p_val_chop>0.05)].index.values
chop_variant = exprs_repl_mean[(p_val_13del<0.05) & (p_val_chop<0.05)].index.values
chop_invariant = genes[genes.index.isin(chop_invariant)]
chop_variant = genes[genes.index.isin(chop_variant)]
# %%

from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()

# %%
# sns.clustermap(scaler.fit_transform(chop_variant.T).T, cmap="vlag", col_cluster=False, method='ward',z_score=0, metric='euclidean')
# %%
chop_invariant_kmeans = KMeans(n_clusters=2, random_state=0).fit(scaler.fit_transform(chop_invariant.T).T)
np.save('chop_invariant_kmeans.npy', chop_invariant_kmeans.labels_)

# %%
chop_invariant_0 = chop_invariant.loc[chop_invariant_kmeans.labels_==0,:].sort_values(('13del',  'Upper'), )
chop_invariant_1 = chop_invariant.loc[chop_invariant_kmeans.labels_==1,:].sort_values(('13delChopKO',  'Upper'), )

p = sns.clustermap(chop_invariant_0, cmap="vlag", col_cluster=False, z_score=0, method='ward', metric='euclidean')
p.savefig("deg.13del_upregulated.chop_invariant_0.pdf")

p = sns.clustermap(chop_invariant_1, cmap="vlag", col_cluster=False, z_score=0, method='ward', metric='euclidean')
p.savefig("deg.13del_upregulated.chop_invariant_1.pdf")

chop_invariant_0.to_csv("deg.13del_upregulated.chop_invariant_0.csv")
chop_invariant_1.to_csv("deg.13del_upregulated.chop_invariant_1.csv")
# %%
chop_variant_kmeans = KMeans(n_clusters=3, random_state=0).fit(scaler.fit_transform(chop_variant.T).T)
np.save('chop_variant_kmeans.npy', chop_variant_kmeans.labels_)
# %%
chop_variant_0 = chop_variant.loc[chop_variant_kmeans.labels_==0,:].sort_values(('13del',  'Upper'), )
chop_variant_1 = chop_variant.loc[chop_variant_kmeans.labels_==1,:].sort_values(('13delChopKO',  'Upper'), )
chop_variant_2 = chop_variant.loc[chop_variant_kmeans.labels_==2,:].sort_values(('13delChopKO',  'Upper'), )
# sns.clustermap(chop_variant_2, cmap="vlag", col_cluster=False, z_score=0, method='ward', metric='euclidean')
print(chop_variant_0.columns)
p = sns.clustermap(chop_variant_0, cmap="vlag", col_cluster=False, z_score=0, method='ward', metric='euclidean')
p.savefig("deg.13del_upregulated.chop_variant_0.pdf")

p = sns.clustermap(chop_variant_1, cmap="vlag", col_cluster=False, z_score=0, method='ward', metric='euclidean')
p.savefig("deg.13del_upregulated.chop_variant_1.pdf")

p = sns.clustermap(chop_variant_2, cmap="vlag", col_cluster=False, z_score=0, method='ward', metric='euclidean')
p.savefig("deg.13del_upregulated.chop_variant_2.pdf")

chop_variant_0.to_csv("deg.13del_upregulated.chop_variant_0.csv")
chop_variant_1.to_csv("deg.13del_upregulated.chop_variant_1.csv")
chop_variant_2.to_csv("deg.13del_upregulated.chop_variant_2.csv")
# # %%
# wt_13del.loc[:,:] = scaler.fit_transform(wt_13del.T).T
p = sns.clustermap(genes.sort_values(('13del',  'Upper'), ), cmap="vlag", col_cluster=False, z_score=0, method='ward', metric='euclidean')
p.savefig("deg.13del_upregulated.pdf")



# # %%
# er = exprs_sig.loc[kmeans.labels_==0,:].copy()
# # %%
# er_chop = er.iloc[:, 4:].copy()
# er_chop.loc[:,:] = scaler.fit_transform(er_chop.T).T
# # %%
# sns.clustermap(er_chop, cmap="vlag", col_cluster=False, z_score=0, method='ward', metric='euclidean')


# # %%
# er_chop_kmeans = KMeans(n_clusters=4, random_state=0).fit(er_chop)
# np.save('er_chop_kmeans_label.npy', er_chop_kmeans.labels_)


# # %%
# sns.clustermap(er_chop.loc[er_chop_kmeans.labels_==3,:], cmap="vlag", col_cluster=False, z_score=0, method='ward', metric='euclidean')

# # %%
# sns.clustermap(exprs_sig, cmap="vlag", col_cluster=False, method='ward',z_score=0, metric='euclidean')


# # %%
# from scipy.stats import ranksums

# chop_invariant = exprs_repl_mean[(p_val_13del<0.05) & (p_val_chop>0.05)]

# # %%
# sns.clustermap(chop_invariant, cmap="vlag", col_cluster=False, method='ward',z_score=0, metric='euclidean')


# # %%
# with open("./chop_invariant.txt", 'w') as f:
#     print('\n'.join(chop_invariant.index.values), file=f)  

# # %%
# chop_variant = exprs_repl_mean[(p_val_13del<0.05) & (p_val_chop<0.05)]

# # %%
# with open("./chop_variant.txt", 'w') as f:
#     print('\n'.join(chop_variant.index.values), file=f)  

# # %%
# sns.clustermap(chop_variant, cmap="vlag", col_cluster=False, method='ward',z_score=0, metric='euclidean')


# # %%


# %%
