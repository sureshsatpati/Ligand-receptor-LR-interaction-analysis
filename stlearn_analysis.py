import pandas as pd
import numpy as np
import stlearn as st
import scanpy as sc
import warnings
import pickle
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
from scipy.sparse import csr_array
adata = st.ReadXenium(feature_cell_matrix_file="/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250116__221632__SCC_BCC_JHJG99_01162025/output-XETG00074__0024835__813_SCC__20250116__221656/cell_feature_matrix.h5", cell_summary_file="/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Xenium/SCC_BCC_FF_IS_IC/20250116__221632__SCC_BCC_JHJG99_01162025/output-XETG00074__0024835__813_SCC__20250116__221656/cells.csv.gz")
adata.X

metadata = pd.read_csv('Final_IS_813_sample.csv')

metadata = metadata.set_index('Cell_IDs')
metadata.index = metadata.index.str.strip()

common_cells = adata.obs_names.intersection(metadata.index)
metadata = metadata.loc[common_cells]

adata.obs.loc[common_cells, "Detailed_Cell_Annotations"] = metadata["Detailed_Cell_Annotations"]

adata = adata[~adata.obs["Detailed_Cell_Annotations"].isna()].copy()

adata.obs["Detailed_Cell_Annotations"] = (adata.obs["Detailed_Cell_Annotations"].astype(float).astype(int).astype(str).astype("category"))

st.em.run_pca(adata,n_comps=50,random_state=0)
st.pp.neighbors(adata,n_neighbors=25,use_rep='X_pca',random_state=0)
st.tl.clustering.louvain(adata,random_state=0)

#fig = plt.figure()
#fig, axes = plt.subplots(ncols=2, figsize=(25,8),dpi = 600)
#st.pl.cluster_plot(adata, use_label='louvain',size=10, ax=axes[0], show_plot=False)
#st.pl.cluster_plot(adata, use_label='Detailed_Cell_Annotations',size=10, ax=axes[1], show_plot=False)


##axes[0].set_title(f'Cell louvain clustering')
#axes[1].set_title(f'Cell combined clustering')
##axes[1].set_title(f'Detailed_Cell_Annotations')
##plt.show()
#fig.savefig('scc_' + target_column_name + '_combined_clustering_gridding.png', dpi=fig.dpi)
##fig.savefig('813_scc_Detailed_Cell_Annotations_gridding.png', dpi=fig.dpi)

fig = plt.figure()
fig, axes = plt.subplots(ncols=2, figsize=(25,8),dpi = 600)
st.pl.cluster_plot(adata, use_label='louvain',size=10, ax=axes[0], show_plot=False)
st.pl.cluster_plot(adata, use_label='Detailed_Cell_Annotations',size=10, ax=axes[1], show_plot=False)

axes[0].set_title(f'Cell louvain clustering')
axes[1].set_title(f'Detailed_Cell_Annotations clustering')
plt.show()
fig.savefig('813_scc_Detailed_Cell_Annotations_clustering_gridding.png', dpi=fig.dpi)


### C2C Communication
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')
print(len(lrs))
# Running the analysis #
st.tl.cci.run(adata,
  lrs,
  min_spots = 20, #Filter out any LR pairs with no scores for less than min_spots
  distance=None, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
  n_pairs=1000, # Number of random pairs to generate; low as example, recommend ~10,000
  n_cpus=None, # Number of CPUs for parallel. If None, detects & use all available.
)

lr_info = adata.uns['lr_summary'] # A dataframe detailing the LR pairs ranked by number of significant spots.
print(lr_info.shape)
print(lr_info)

lr_info.to_csv('813_LR_pairs_ranking.csv')
# Showing the rankings of the LR from a global and local perspective.
# Ranking based on number of significant hotspots.
st.pl.lr_summary(adata, n_top=500)
st.pl.lr_summary(adata, n_top=20, figsize=(10,3))

### Can adjust significance thresholds.
st.tl.cci.adj_pvals(adata,
correct_axis='spot',
pval_adj_cutoff=0.05,
adj_method='fdr_bh')

#st.tl.cci.run_cci(adata, 'combined_clustering', # Spot cell information either in data.obs or data.uns
st.tl.cci.run_cci(adata, 'Detailed_Cell_Annotations', # Spot cell information either in data.obs or data.uns
  min_spots=2, # Minimum number of spots for LR to be tested.
  spot_mixtures=True, # If True will use the deconvolution data,
  # so spots can have multiple cell types if score>cell_prop_cutoff
  cell_prop_cutoff=0.1, # Spot considered to have cell type if score>0.1
  sig_spots=True, # Only consider neighbourhoods of spots which had significant LR scores.
  n_perms=100, # Permutations of cell information to get background, recommend ~1000
  n_cpus=None,
)

#fig, ax, ax2 = st.pl.cci_check(adata, 'combined_clustering', figsize=(16,5),show = False)
fig, ax, ax2 = st.pl.cci_check(adata, 'Detailed_Cell_Annotations', figsize=(16,5),show = False)
fig.savefig('813_scc_CCI_LR.png', dpi=600)

y_data = [line.get_ydata() for line in ax2.get_lines()]
y_data = np.concatenate(y_data)
# [array([ 19,   2,  44,  23,  29,  47,  61,  70,  60,  36,  16, 120,  46])]
text_labels = [text.get_text() for text in ax.texts]
# ['B_Cells', 'Melanocyte', 'Exhaustive_CD8+', 'Tregs', 'CD4+', 'Monocytes', 'CD8+', 'Endothelial', 'DCs', 'Fibroblasts', 'Skin', 'All_Macrophages', 'Tumor']
df = pd.DataFrame({'Text Label': text_labels, 'Y Data': y_data})

# Save DataFrame as CSV
df.to_csv('813_scc_cci_lr_numbers.csv', index=False)

#fig, ax = st.pl.lr_chord_plot(adata, 'combined_clustering',show = False)
fig, ax = st.pl.lr_chord_plot(adata, 'Detailed_Cell_Annotations',show = False)
fig.savefig('813_scc_CCI_LR_chord_plot_all.png', dpi=600)

lrs = adata.uns['lr_summary'].index.values[0:10]
#for lr in lrs[0:10]:
#fig, ax = st.pl.lr_chord_plot(adata, 'combined_clustering', lr,show = False,min_ints=0)
#fig.savefig(target_column_name + '_CCI_LR_' + lr + '_chord_plot_all.png', dpi=600)


with open('813_scc_uns.pkl', 'wb') as fp:
    pickle.dump(adata.uns, fp)
print('dictionary saved successfully to file')

del adata.uns

adata.write_h5ad("813_scc_.h5ad",compression='gzip')
