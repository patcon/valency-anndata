# ============================================================
# ESS Human Values — Colab Exploration
# Copy each cell section into a separate Colab cell.
# ============================================================

# %% [markdown]
# ## 1. Install & setup

# %%
# !pip install -q git+https://github.com/patcon/valency-anndata.git@feat/ess-export
# !curl -sO https://raw.githubusercontent.com/patcon/valency-anndata/feat/ess-export/scripts/ess_to_adata.py


# %% [markdown]
# ## 2. Generate h5ad

# %%
ESS_USER_ID = "691c9862-aab9-4e05-9ed1-ca47fb4f7b70"

# !python ess_to_adata.py --user-id {ESS_USER_ID}
# Writes: exports/ess_round11_human_values.h5ad


# %% [markdown]
# ## 3. Load

# %%
import anndata as ad
import valency_anndata as val

adata = ad.read_h5ad("exports/ess_round11_human_values.h5ad")
print(adata)
print(adata.var[["content"]].to_string())


# %% [markdown]
# ## 4. Pipeline: Impute → PCA → k-means

# %%
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

val.pp.impute(adata, strategy="mean", target_layer="X_imputed")

X = adata.layers["X_imputed"]
pca = PCA(n_components=10, random_state=42)
adata.obsm["X_pca"] = pca.fit_transform(X)
print("Variance explained:", [f"PC{i+1}={r:.1%}" for i, r in enumerate(pca.explained_variance_ratio_[:5])])

for k in range(2, 6):
    km = KMeans(n_clusters=k, random_state=42, n_init=10)
    adata.obs[f"kmeans_{k}"] = km.fit_predict(adata.obsm["X_pca"]).astype(str)


# %% [markdown]
# ## 5. PaCMAP + LocalMAP

# %%
val.tl.pacmap(adata, layer="X_imputed")
val.tl.localmap(adata, layer="X_imputed")


# %% [markdown]
# ## 6. Scatter Plots

# %%
import matplotlib.pyplot as plt
import scanpy as sc

for basis, title in [("X_pacmap", "PaCMAP"), ("X_localmap", "LocalMAP"), ("X_pca", "PCA")]:
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(title, fontsize=14)
    for ax, color in zip(axes, ["kmeans_2", "kmeans_4", "country"]):
        sc.pl.embedding(adata, basis=basis, color=color, ax=ax, show=False,
                        s=3, title=color)
    plt.tight_layout()
    plt.show()


# %% [markdown]
# ## 7. Cluster Profiles
# Mean vote per variable per cluster — shows what each cluster cares about.

# %%
import seaborn as sns
import pandas as pd

K = 4  # change to match whichever k looked most meaningful above

cluster_means = (
    pd.DataFrame(adata.X, columns=adata.var_names)
    .assign(cluster=adata.obs[f"kmeans_{K}"].values)
    .groupby("cluster")
    .mean()
    .T
)
cluster_means.index = [adata.var["content"][v][:45] for v in cluster_means.index]

fig, ax = plt.subplots(figsize=(8, 7))
sns.heatmap(
    cluster_means,
    center=0, vmin=-1, vmax=1,
    cmap="RdYlGn",
    linewidths=0.5,
    annot=True, fmt=".2f",
    ax=ax,
)
ax.set_title(f"Mean vote per variable — k={K} clusters")
ax.set_xlabel("Cluster")
plt.tight_layout()
plt.show()


# %% [markdown]
# ## 8. Color by individual values (optional)
# Shows how agreement with each specific value is distributed across the embedding.

# %%
basis = "X_pacmap"  # or "X_localmap"

n_cols = 3
n_rows = -(-len(adata.var_names) // n_cols)  # ceiling division
fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
axes = axes.flatten()

for ax, var_name in zip(axes, adata.var_names):
    sc.pl.embedding(adata, basis=basis, color=var_name, ax=ax, show=False,
                    s=2, title=adata.var["content"][var_name][:40],
                    vmin=-1, vmax=1, cmap="RdYlGn")

for ax in axes[len(adata.var_names):]:
    ax.set_visible(False)

plt.tight_layout()
plt.show()
