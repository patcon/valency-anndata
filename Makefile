install-dev: ## Install dev dependencies
	uv sync --extra dev

install-all: ## Install all dependencies
	uv sync --all-extras

strip-nb-widgets: ## Strip ipywidget metadata from notebooks so they render on GitHub
	uv run python scripts/strip_widget_metadata.py docs/notebooks/*.ipynb

notebook-docs: ## Build the markdown versions of notebook files for docs site
	IS_GENERATING_DOCS=true uv run jupyter nbconvert docs/notebooks/*.ipynb \
		--config jupyter_nbconvert_config.py \
		--NbConvertApp.output_files_dir notebook-assets \
		--to markdown

notebook-docs-debug:
	IS_GENERATING_DOCS=true uv run jupyter nbconvert docs/notebooks/*.ipynb \
		--config jupyter_nbconvert_config.py \
		--NbConvertApp.output_files_dir notebook-assets \
		--log-level=DEBUG \
		--to markdown

serve: install-dev ## Serve documentation website for local development
	uv run mkdocs serve

docs: install-dev ## Build documentation website directory
	uv run mkdocs build

lint: ## Lint source with ruff
	uv run ruff check src/

fmt: ## Format source with ruff
	uv run ruff format src/

test: ## Run unit and local-fixture tests (excludes live network tests)
	uv run pytest

test-live: ## Run only the live network tests (requires pol.is access)
	uv run pytest -m live

csv-export: ## Export Polis CSV files (usage: make csv-export URL=https://pol.is/2fzs33rmbr)
	@test -n "$(URL)" || (echo "Usage: make csv-export URL=<polis-url-or-id>" && exit 1)
	$(eval SLUG := $(lastword $(subst /, ,$(URL))))
	uv run python scripts/export_csv.py "$(URL)" "exports/$(SLUG)"

h5ad: ## Export a Polis conversation to h5ad (writes to exports/)
	uv run python scripts/export_h5ad.py

release: ## Print no-op documentation to guide the release process
	uv run python scripts/release.py

js: ## Build the TypeScript widget bundles into src/valency_anndata/viz/static/
	cd js && pnpm install --frozen-lockfile && pnpm build

js-watch: ## Rebuild widget bundles on change (pair with ANYWIDGET_HMR=1 in Jupyter)
	cd js && pnpm watch

build: ## Build wheel package for publishing to PyPI
	rm -rf dist/
	uv build

publish: ## Publish built package to PyPI
	uv publish

# These make tasks allow the default help text to work properly.
%:
	@true

.PHONY: help js js-watch strip-nb-widgets notebook-docs notebook-docs-debug serve docs lint fmt test test-live csv-export h5ad release build publish

help:
	@echo 'Usage: make <command>'
	@echo
	@echo 'where <command> is one of the following:'
	@echo
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.DEFAULT_GOAL := help
