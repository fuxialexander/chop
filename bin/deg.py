#!/usr/bin/env python
# %%
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# %%
import seaborn as sns
#%%
from scipy.stats import ranksums
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600

# %%
exprs = pd.read_csv("exprs.ann.txt", sep="\t")
# %%
exprs_filt = exprs.loc[np.sum(exprs.iloc[:, 4:].values > np.log2(300), 1
                              ) > 0, :].iloc[:, 2:].dropna()
exprs_filt.columns = exprs_filt.columns.str.replace('.CEL.gz', '').str.replace(
    'GSM[0-9]+_[0-9]+_', '')
exprs_filt = exprs_filt.groupby('SYMBOL').max()
p_val_13del = np.array([
    ranksums(i[1].values[8:14], i[1].values[20:29])[1]
    for i in exprs_filt.iterrows()
])

upr = (exprs_filt.iloc[:, 8:14].mean(1).values <
       exprs_filt.iloc[:, 20:29].mean(1).values)

exprs_filt_sig_upr = exprs_filt.iloc[upr & (p_val_13del < 0.05), :]

p_val_chop = np.array([
    ranksums(i[1].values[33:39], i[1].values[20:29])[1]
    for i in exprs_filt_sig_upr.iterrows()
])

chop_invariant = exprs_filt_sig_upr[p_val_chop > 0.05].index.values
chop_variant = exprs_filt_sig_upr[p_val_chop < 0.05].index.values

# %%
replicates_map = dict(
    zip(exprs_filt_sig_upr.columns[2:].values,
        exprs_filt_sig_upr.columns[2:].str.replace('_rep.', '').values))

# %%
exprs_repl_mean = exprs_filt_sig_upr.iloc[:, 2:].groupby(by=replicates_map,
                                                         axis=1).mean().copy()

exprs_repl_mean.index.names = ['Genes']
# %%
new_mul_index = pd.MultiIndex.from_tuples(
    [tuple(i.split('_'))[::-1] for i in exprs_repl_mean.columns.values])

exprs_repl_mean.columns = new_mul_index
# %%
col_order = [('WT', 'PZ'), ('WT', 'pHZ'), ('WT', 'Upper'), ('WT', 'Lower'),
             ('13del', 'PZ'), ('13del', 'pHZ'), ('13del', 'Upper'),
             ('13del', 'Middle'), ('13del', 'Lower'), ('13delChopKO', 'PZ'),
             ('13delChopKO', 'PHZ'), ('13delChopKO', 'Upper'),
             ('13delChopKO', 'Middle'), ('13delChopKO', 'Lower')]

# %%
exprs_repl_mean = exprs_repl_mean.reindex(col_order, axis=1)
exprs_repl_mean.columns.names = ['Treatment', 'Segment']
# %%

chop_invariant = exprs_repl_mean[exprs_repl_mean.index.isin(chop_invariant)]
# %%
chop_variant = exprs_repl_mean[exprs_repl_mean.index.isin(chop_variant)]
# exprs_repl_mean.to_csv("exprs.formated.csv")

# %%
exprs_repl_mean.to_csv("deg.wt_13del.13del_upregulated.csv")


# %%
# sns.clustermap(exprs_sig.iloc[:,0:9], cmap="vlag", col_cluster=False, z_score=0, method='ward', metric='euclidean')
# %%

scaler = StandardScaler()

# %%
# import seaborn as sns
# sns.clustermap(scaler.fit_transform(chop_invariant.T).T, cmap="vlag", col_cluster=False, method='ward',z_score=0, metric='euclidean')
# %%
chop_invariant_kmeans = KMeans(n_clusters=2, random_state=0).fit(
    scaler.fit_transform(chop_invariant.T).T)
# np.save('chop_invariant_kmeans.npy', chop_invariant_kmeans.labels_)

# %%

# def get_new_gene_label(df):
#     heatmap_genes = pd.read_csv('../heatmap_genes.txt', header=None).values
#     return [i if i in heatmap_genes else '                 ' for i in df.index.values]
    
def plot_clustermap(df, fname, xticklabels=False, cbar=False, total=1178):
    p = sns.clustermap(df,
                   cmap="vlag",
                   col_cluster=False,
                   z_score=0,
                   method='ward',
                   metric='euclidean',
                #    yticklabels=get_new_gene_label(df), 
                   xticklabels=xticklabels,
                   figsize=(5,8*df.shape[0]/1178.))
                   
    p.ax_heatmap.tick_params(axis='both', which='both', length=0)
    p.ax_heatmap.set_yticklabels(p.ax_heatmap.get_yticklabels(), rotation=0)
    p.ax_heatmap.set_ylabel('')
    p.ax_heatmap.set_xlabel('')
    p.cax.set_visible(cbar)
    print(p)
    p.savefig(fname+'.pdf', )
    df.to_csv(fname+'.csv')
    return 1
# %%

chop_invariant_0 = chop_invariant.loc[chop_invariant_kmeans.labels_ ==
                                      0, :].sort_values(('13del', 'Upper'), )
chop_invariant_1 = chop_invariant.loc[chop_invariant_kmeans.labels_ ==
                                      1, :].sort_values(
                                          ('13delChopKO', 'Upper'), )

plot_clustermap(chop_invariant_0, "deg.13del_upregulated.chop_invariant_0")
plot_clustermap(chop_invariant_1, "deg.13del_upregulated.chop_invariant_1")

# %%
chop_variant_kmeans = KMeans(n_clusters=3, random_state=0).fit(
    scaler.fit_transform(chop_variant.T).T)
# np.save('chop_variant_kmeans.npy', chop_variant_kmeans.labels_)
# %%
chop_variant_0 = chop_variant.loc[chop_variant_kmeans.labels_ ==
                                  0, :].sort_values(('13del', 'Upper'), )
chop_variant_1 = chop_variant.loc[chop_variant_kmeans.labels_ ==
                                  1, :].sort_values(('13delChopKO', 'Upper'), )
chop_variant_2 = chop_variant.loc[chop_variant_kmeans.labels_ ==
                                  2, :].sort_values(('13delChopKO', 'Upper'), )
# sns.clustermap(chop_variant_2, cmap="vlag", col_cluster=False, z_score=0, method='ward', metric='euclidean')


plot_clustermap(chop_variant_0, "deg.13del_upregulated.chop_variant_0")
plot_clustermap(chop_variant_1, "deg.13del_upregulated.chop_variant_1")
plot_clustermap(chop_variant_2, "deg.13del_upregulated.chop_variant_2")
# %%
plot_clustermap(exprs_repl_mean.sort_values(('13del', 'Upper')), "deg.13del_upregulated")


# %%
