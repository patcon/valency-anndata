# Usage Principles

Import `valency-anndata` as:

```py
import valency_anndata as val
```

## Workflow

The typical workflow involves:

1. Loading raw vote data via methods in `val.datasets`.
2. Repeated rounds of...
    1. Processing data via methods in `val.preprocessing` and `val.tools`.
    2. Visualizing the processed data via methods in `val.viz`.

So we could import data from a Polis conversation's report URL, for e.g.,

```py
adata = val.datasets.polis.load("<some report url>", **import_params)
val.tools.recipe_polis(adata, **tool_params)
```

where `adata` is an [`AnnData`][anndata.AnnData] object (see below).

- `val.dataset` import-related methods populate the _valence matrix_ (aka _vote matrix_) `.X`, as well as add basic annotations to either participant rows in `.obs` or statement columns in `.var`.
    - These basic annotations might be the author ID or text content of a statement in `.var`, or the xid or sign-up timestamp of a participant in `.obs`.
- `val.tools` methods add additional calculated annotations to the valence matrix `.X`.
    - Simple annotations (strings or numbers) often get stored in `.obs` or `.var`, but more complex annotations have other places they can be stored.
    - New representations (of the high dimensional participant data in `.X`) when projected into lower-dimensional space tend to be stored in `.obsm`. For e.g., `.obsm["X_pca"]` might hold the representation of the first 50 components of a PCA projection of `.X`, or `.obsm["X_umap"]` might hold the representation of 2 or 3 components of a UMAP projection.

```py
val.viz.embedding(adata, **plotting_params)
```

To facilitate writing memory-efficient pipelines, by default, `valency-anndata` tools operate _inplace_ on `adata` and return `None`. If you want to return a copy of the [`AnnData`][anndata.AnnData] object and leave the passed `adata` unchanged, pass `inplace=False` (or, less often, `copy=True`).

## Underlying Tools

`valency-anndata` is building on the conventions and some methods of [Scanpy](https://scanpy.readthedocs.io/), which in turn builds on data structures of [`anndata`](https://anndata.readthedocs.io/)

### Scanpy

_To be written._

### AnnData

_To be written._

## Terminology

