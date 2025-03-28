import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


adata = sc.read("/home/amomtaz/Beng199_Projects/scGPT/MacaqueStudyData/NeOrgAtlas.h5ad")  # Replace with your file
#duplicates = adata.var["gene_name"][adata.var["gene_name"].duplicated(keep=False)]
#print(adata.obs)

print(adata.rows.unique())

'''# Convert to DataFrame
df = pd.DataFrame(adata.X.toarray(), columns=adata.var["gene_name"], index=adata.obs_names)
print(adata.shape)

# Sum duplicate gene names
df_grouped = df.groupby(axis=1, level=0).sum()
print(df_grouped.shape)

# Assign back to adata
adata = adata[:, adata.var["gene_name"].drop_duplicates(keep="first").index]  # Keep only one instance per gene
adata.X = sp.csr_matrix(df_grouped.values)  # Assign summed values
adata.var = adata.var.loc[df_grouped.columns].copy()  # Update var to reflect merged genes
'''

'''
# Open an adata
adata = sc.read("/home/amomtaz/Beng199_Projects/scGPT/MacaqueStudyData/macaque_adata.h5ad")

# Normalize the data 
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.X = adata.X

# Group cell types with less than 600 cells as "Other"
celltype_counts = adata.obs["PredCellType"].value_counts()
rareTypes = celltype_counts[celltype_counts < 600].index
adata.obs["PredCellType"] = adata.obs["PredCellType"].replace(rareTypes, "Other")
print(adata.obs["PredCellType"].value_counts())

# Create a PCA, then UMAP, then save the UMAP plot as a new image file
sc.pp.pca(adata, n_comps=30)  # Reduce to 30 dimensions
sc.pp.neighbors(adata)
sc.tl.umap(adata,)
sc.pl.umap(adata,color="PredCellType", save="_umap.png")

# Write the adata as a new file
adata.write("clus_mac_adata.h5ad")

# Select highly variable genees
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
adata_hvg = adata[:, adata.var["highly_variable"]]

# Print the dimensions of the PCA and UMAP
print(adata.obsm["X_pca"].shape)
print(adata.obsm["X_umap"].shape)

Notes:
adata.obs holds cell barcode, other info (e.g. species), and celltype
  - This is basically the CO_4Mo_TrueCellTypes.tsv
Need deal with the duplicated gene names so that I can concatenate the adatas
  duplicates = adata.var["gene_name"][adata.var["gene_name"].duplicated(keep=False)]
  print(duplicates)
'''

 