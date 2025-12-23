<p align="center"><img src="https://imgur.com/UfnMqH0.png" alt="Logo" width=300></p>

# valency-anndata

Experimental tooling to support notebook analysis of polislike data.

:sparkles: Inspired by [scanpy][] and the [scverse][] ecosystem! :heart:

## Installation

```
pip install git+https://github.com/patcon/valency-anndata
```

## Usage

### Loading Polis Data

```py
import valency_anndata as val

adata = val.datasets.polis.load("https://pol.is/report/r29kkytnipymd3exbynkd")
val.viz.schematic_diagram(adata, diff_from=None)
```
<img height="200" alt="Screenshot 2025-12-23 at 12 00 16 AM" src="https://github.com/user-attachments/assets/af10a11a-0146-401e-afff-6567255cf51a" />

### Running Polis Pipelines

```py
with val.viz.schematic_diagram(diff_from=adata):
    val.tools.recipe_polis(adata, key_added_pca="X_pca_polis")
```
<img height="200" alt="Screenshot 2025-12-23 at 12 57 42 AM" src="https://github.com/user-attachments/assets/9341387b-0358-4f1e-bb27-3d2cd4beef7d" />


```py
val.viz.embedding(adata, basis="pca_polis", colors="kmeans_polis")
```
<img height="200" alt="Screenshot 2025-12-23 at 12 00 59 AM" src="https://github.com/user-attachments/assets/7cfe76d5-a03f-4024-bfe1-d152747845e4" />


### Exploring Polis Pipelines

```py
val.viz.schematic_diagram(diff_from=adata):
    val.preprocessing.calculate_qc_metrics(pacmap_adata, inplace=True)
```
<img height="200" alt="Screenshot 2025-12-23 at 12 58 18 AM" src="https://github.com/user-attachments/assets/7ef005b2-f6c3-4c7f-bebb-6c37c9a14290" />


```py
val.viz.embedding(adata, basis="pca_polis",
    colors=["kmeans_polis", "pct_seen", "pct_agree", "pct_pass"],
)
```
<img height="200" alt="Screenshot 2025-12-23 at 12 58 50 AM" src="https://github.com/user-attachments/assets/18133b45-cd2b-41a9-a7c5-2101426ba1f9" />


### Running & Exploring Alternative Pipelines

```py
from valency_anndata.tools._polis import _zero_mask, _cluster_mask

with val.viz.schematic_diagram(diff_from=adata):
    _zero_mask(adata)
    val.preprocessing.impute(
        adata,
        strategy="mean",
        source_layer="X_masked",
        target_layer="X_masked_imputed_mean",
    )
    val.tools.pacmap(
        adata,
        key_added="X_pacmap",
        layer="X_masked_imputed_mean",
    )
    _cluster_mask(adata)
    val.tools.kmeans(
        adata,
        k_bounds=(2, 9),
        use_rep="X_pacmap",
        mask_obs="cluster_mask",
        key_added="kmeans_pacmap",
    )
```
<img height="200" alt="Screenshot 2025-12-23 at 12 58 59 AM" src="https://github.com/user-attachments/assets/1fbed051-9a66-4e0e-82a2-f58c839d1a06" />


```py
val.scanpy.pl.embedding(adata, basis="pacmap",
    color=["kmeans_pacmap", "pct_seen", "pct_agree", "pct_pass"],
)
```
<img height="200" alt="Screenshot 2025-12-23 at 12 59 09 AM" src="https://github.com/user-attachments/assets/824ec281-d7fa-4f22-9b80-785119cf1529" />


For full examples and planned features, see: [`example-usage.ipynb`](./example-usage.ipynb)

## Contributing

We are maintaining a custom [`CONTRIBUTING.md`](./CONTRIBUTING.md) with specific links and a compiled list of entry tasks!

<!-- Links -->
   [scanpy]: https://scanpy.readthedocs.io/

   [scverse]: https://scverse.org/
