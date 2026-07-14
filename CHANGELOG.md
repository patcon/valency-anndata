# Changelog

## [Unreleased][] (YYYY-MM-DD)

### Fixes
- `val.viz.schematic_diagram()` — fix display in marimo notebooks. It previously only detected Jupyter/IPython, so in marimo it fell through to `webbrowser.open()`, which raised `webbrowser.Error` when no GUI browser was reachable. Now detects marimo via `marimo.running_in_notebook()` and displays inline via `marimo.output.replace(marimo.Html(...))`.

## [0.4.0][] (2026-07-08)

### Added
- `val.viz.heatmap()` — plot a vote-matrix heatmap with Polis-friendly defaults. Wraps :func:`scanpy.pl.heatmap` with an optional `discrete=True` flag for a labelled segmented colorbar (`"disagree (-1)"`, `"pass (0)"`, `"agree (+1)"`), an optional `groupby` (defaults to index order when omitted), and a built-in `"RdYlGnBright"` colormap ([#92][]).
- `val.datasets.polis.load()` — new `include_precomputed_groups=True` flag stores Polis server's native group assignments in `adata.obs["kmeans_polis_precomputed"]` (nullable `Int64`) and the raw math dict in `adata.uns["polis_math"]`. Enables easy comparison of pipeline clustering vs. Polis-native grouping ([#93][]).
- `val.datasets.polis.load()` — new `trim_rule` parameter (default `1.0`, keeping all votes) passed through to `val.pp.rebuild_vote_matrix()`. Supports fractional, percent-integer, absolute-timestamp, and statistical-string trim rules. All reference dataset loaders (`aufstehen`, `vtaiwan`, `american_assembly`, `bg2050`, `chile_protest`, `cuba_protest`, `japanchoice`, `klimarat`) now accept `**kwargs` and forward them to `polis.load()`.

### Fixes
- `val.pp.rebuild_vote_matrix()` — fix `ValueError` raised when `trim_rule` reduced the matrix shape and the original `layers` or `obsm` entries no longer matched the new dimensions. Mismatched entries are now silently dropped.
- Pin `anndata<0.13.0` — anndata 0.13.0 introduces breaking changes that are incompatible with this package.

[#92]: https://github.com/patcon/valency-anndata/issues/92
[#93]: https://github.com/patcon/valency-anndata/issues/93

## [0.3.0][] (2026-03-04)

### Added
- `val.pp.filter_participants()` — filter participants (rows) by minimum number of statements voted on. Counts non-`NaN` entries (real votes), correctly treating `-1`, `0`, and `+1` as votes.
- `val.pp.filter_statements()` — filter statements (columns) by minimum number of participants who voted. Counts non-`NaN` entries (real votes), correctly treating `-1`, `0`, and `+1` as votes.
- `val.datasets.vtaiwan()` — load any of four Polis conversations from Taiwan's vTaiwan civic policymaking process, selectable by `topic=` keyword (`"uber"`, `"airbnb"`, `"online_alcohol"`, `"caning"`). Topic parameter uses `Literal` type hint for IDE/notebook autocomplete.
- `val.datasets.american_assembly()` — load Polis conversations run by the American Assembly in Kentucky cities, selectable by `city=` keyword (`"bowling_green"`, `"louisville"`). City parameter uses `Literal` type hint for IDE/notebook autocomplete.
- `val.datasets.bg2050()` — load the BG 2050 community visioning conversation from Bowling Green and Warren County, Kentucky (~7,900 participants).
- `val.datasets.cuba_protest()` — load any of three Polis conversations run around Cuba's planned 15N march (November 2021), selectable by `period=` keyword (`"before_1"`, `"before_2"`, `"after"`). Period parameter uses `Literal` type hint for IDE/notebook autocomplete.
- `val.datasets.japanchoice()` — load any of eight Polis conversations from Japan Choice (four policy topics × two election years: 2025 and 2026), selectable by positional `topic` argument (`"2025_foreign_affairs_security"`, `"2025_diversity_human_rights"`, `"2025_education_children_old_age"`, `"2025_economy_taxation_employment"`, and 2026 equivalents). Topic parameter uses `Literal` type hint for IDE/notebook autocomplete.
- `val.datasets.klimarat()` — load any of the five Polis conversations from Austria's Citizens' Climate Council (Klimarat), selectable by `topic=` keyword (`"food_land"`, `"mobility"`, `"energy"`, `"housing"`, `"production"`). Topic parameter uses `Literal` type hint for IDE/notebook autocomplete.
- Five Klimarat datasets added to the docs overview table with fingerprints: Food & Land Use, Mobility, Energy, Housing, Production & Consumption.
- "Data Overview" promoted to a top-level section on the datasets API page; `val.datasets.klimarat` added as a documented subsection alongside `aufstehen` and `chile_protest`.
- `scripts/generate_fingerprint_heatmap.py` — generates a square RdYlGn vote-matrix heatmap from any Polis report URL ([#89][]). Supports `--participant-vote-threshold`, `--statement-vote-threshold`, `--open`, and custom output path.
- `docs/api/datasets.yml` — machine-readable registry of reference datasets, rendered into an overview table in `docs/api/datasets.md` via mkdocs-macros ([#89][]).
- Two reference datasets added to the docs table: **Aufstehen** and **Chile Protests** ([#89][]).
- `mkdocs-glightbox` — clicking a fingerprint thumbnail opens the full-size image in a lightbox popup ([#89][]).
- New **Labs** page in docs — a curated showcase of experimental notebooks and apps built on valency-anndata, including Oval, Polis Report Processor, Perspective Map Explorer v2, Semantic Statement Map Creator, and more ([#84][]).
- `make strip-notebook-widgets` — strips ipywidget metadata from notebooks so they render correctly on GitHub ([#31][]).

### Changed
- `val.preprocessing.impute()` now uses `sklearn.impute.SimpleImputer` for `"zero"`, `"mean"`, and `"median"` strategies instead of a custom numpy implementation ([#71][]). Adds `strategy="knn"` backed by `sklearn.impute.KNNImputer`, with `n_neighbors` and `weights` parameters.
- `val.preprocessing.highly_variable_statements()` defaults changed: `variance_mode` is now `"valence"` (was `"overall"`), `bin_by` is now `"p_engaged"` (was `"coverage"`), and `n_bins` is now `10` (was `1`).

### Fixes
- `val.preprocessing.highly_variable_statements()` no longer emits `RuntimeWarning: Degrees of freedom <= 0 for slice` when a statement column has fewer than 2 non-NaN votes ([#86][]).
  - Variance is now computed only on columns with ≥ 2 observations; under-observed columns return `NaN`.
  - Added `TestVarianceNumerics` tests to verify computed values match `np.nanvar(ddof=1)` directly.
  - Added `TestNoRuntimeWarnings*` and `TestNoAnyWarnings*` regression tests covering all public API methods (`val.pp`, `val.tl`, `val.viz`) on real fixture data.
- Bugfix: scaling factors in `recipe_polis` were dividing instead of multiplying!

## [0.2.0][] (2026-02-16)

### Added
- `hf:` and `huggingface:` source prefixes for `val.datasets.polis.load()` — load any HuggingFace-hosted Polis export as a one-liner, e.g. `load("hf:patcon/polis-aufstehen-2018")` ([#81][]).
- `CLAUDE.md` guidance file for Claude Code contributors ([#58][]).
- Pytest infrastructure and test suite for `datasets.polis.load` ([#59][]).
  - 29 unit + local-fixture tests; 4 opt-in live network tests (`make test-live`).
  - Synthetic and real CSV fixtures checked in under `tests/fixtures/`.
  - `make test` and `make test-live` targets added to Makefile.
- Unit and integration tests for `tools.kmeans` ([#63][]).
  - 22 mocked unit tests + 1 real-clustering integration test.
  - 3 `k-means++` smoke tests.
- `val.tl.recipe_polis2_statements()` — embeds and clusters statements (var axis) via polismath ([#44][]).
  - New `polis2` optional-dependency group (`pip install valency-anndata[polis2]`).
  - 13 unit tests with all polismath helpers mocked.
  - Noise/unassigned cluster labels (`-1`) in `evoc_polis2_top` are stored as `NA` so scanpy renders them as lightgray by default.
  - `show_progress=False` (the default) now fully silences HF download progress bars and mlx model-load stdout.
  - "Polis 2.0 Pipeline" tutorial added to docs nav.
- `val.preprocessing.highly_variable_statements()` — identify highly variable statements in vote matrices ([#52][]).
  - Analogous to scanpy's highly_variable_genes for single-cell data.
  - Supports multiple variance modes (overall, valence, engagement) and binning strategies.
  - `key_added` parameter allows running multiple times with different settings.
  - `val.viz.highly_variable_statements()` plotting function for visualizing dispersion metrics.
  - `mask_var` parameter added to `val.tools.recipe_polis()`, `val.tools.pacmap()`, and `val.tools.localmap()` for filtering statements before dimensionality reduction.
- `val.write()` — export AnnData to h5ad with automatic sanitization for webapp compatibility ([#57][]).
  - `include` parameter for selective export using glob-style `"namespace/key"` paths (e.g. `"obsm/X_*"`).
- `make lint` and `make fmt` targets for ruff.
- Claude Code skill for guided Polis conversation exploration ([#42][]).
  - Interactive prompts for projection selection (PaCMAP, LocalMAP, UMAP, t-SNE) and QC annotation selection.
  - Fixed CLI plotting to support multi-color `val.viz.embedding()` calls.
- Cache downloaded Polis report files locally for 24 hours using `platformdirs` ([#70][]).
  - `skip_cache` parameter on `val.datasets.polis.load()` to bypass the cache.
  - Smart cache revalidation using `last_vote_timestamp` from the Polis math endpoint — stale cache is reused without re-fetching when no new votes have been cast ([#78][]).
- `mask_obs` parameter on `val.tools.kmeans()` for clustering a subset of participants ([#77][]).
- `val.datasets.polis.export_csv()` — export an AnnData object to Polis CSV format (`votes.csv` + `comments.csv`).
- `include_huggingface_metadata` parameter on `val.datasets.polis.export_csv()` — opt-in generation of a HuggingFace dataset card (`README.md` with YAML frontmatter) alongside the CSV export.
- `show_progress` parameter on `val.datasets.polis.load()` — displays a tqdm progress bar when fetching votes per-participant from the API; auto-detects notebooks vs terminal ([#79][]).

### Fixes
- Fixed `uns["statements"]` having `comment-id` as both index and column, which prevented h5ad serialization ([#57][]).
- Fixed API vote sign inversion — the Polis API returns inverted vote signs vs the CSV export convention; votes are now negated on ingest so `+1` = agree and `-1` = disagree everywhere.
- Replaced deprecated `use_highly_variable=False` with `mask_var=None` in `recipe_polis` PCA call to eliminate `FutureWarning` from scanpy ([#82][]).

[#82]: https://github.com/patcon/valency-anndata/issues/82
[#81]: https://github.com/patcon/valency-anndata/issues/81
[#42]: https://github.com/patcon/valency-anndata/issues/42
[#57]: https://github.com/patcon/valency-anndata/issues/57
[#44]: https://github.com/patcon/valency-anndata/issues/44
[#52]: https://github.com/patcon/valency-anndata/issues/52
[#58]: https://github.com/patcon/valency-anndata/pull/58
[#59]: https://github.com/patcon/valency-anndata/issues/59
[#63]: https://github.com/patcon/valency-anndata/pull/63
[#70]: https://github.com/patcon/valency-anndata/issues/70
[#77]: https://github.com/patcon/valency-anndata/pull/77
[#79]: https://github.com/patcon/valency-anndata/issues/79
[#78]: https://github.com/patcon/valency-anndata/pull/78
[#89]: https://github.com/patcon/valency-anndata/pull/89
[#86]: https://github.com/patcon/valency-anndata/pull/86
[#84]: https://github.com/patcon/valency-anndata/pull/84
[#31]: https://github.com/patcon/valency-anndata/issues/31
[#71]: https://github.com/patcon/valency-anndata/issues/71

## [0.1.1][] (2026-01-20)

### Fixes
- Fixed the README image so not broken on PyPI.

## [0.1.0][] (2026-01-19)

Initial release includes:

- `val.viz.schematic_diagram()` helper to showing data structure and visual diffs.
- Helper methods for PaCMAP and LocalMAP dimensional reduction.
- Langevitour visualisation for exploring high dimensional space.
- Basic Jupyter Scatter support for up to 1M participants, including animations.
- Import of Polis conversation data.
- Basic Polis v1 pipeline support.
- Added `val.tools.kmeans()`.
- Large reference datasets for Aufstehen political party consultation (33k participants) and the #ChileDesperto protest (3k).
- `val.viz.voter_vignette_widget()` for exploring data stories of random individuals.
- Comprehensive documentation website.
- Wrappers for various scanpy methods.

<!-- Links -->

[Unreleased]: https://github.com/patcon/valency-anndata/compare/v0.4.0...main
[0.4.0]: https://github.com/patcon/valency-anndata/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/patcon/valency-anndata/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/patcon/valency-anndata/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/patcon/valency-anndata/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/patcon/valency-anndata/releases/tag/v0.1.0
