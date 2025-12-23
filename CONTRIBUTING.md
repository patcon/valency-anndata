# How To Contribute

There are perhaps two different paths, based on your personal preferences:
- A. Coding: Get into the code and get to a small pull request ASAP, or
- B. Learning: Spend some time learning the context and ecosystem of the project.

Both are great approaches, depending on whether you like working inside-to-out (actions, then context) or outside-to-in (context, then action).

## A. Action: Code First

Small things we need done: (these will move to pull requests soon)
- [ ] some simple tests of existing functionality
- [ ] CI running on GitHub actions, running a basic install on a few versions of Python
- [ ] a simple documentation website using mkdocs (e.g., front-page from README + API page)
- [ ] write a dataset processor for importing representative likert survey data into valency-anndata
    - See "Data Loading" section of [`example-usage.ipynb`](./example-usage.ipynb)
- [x] add [PaCMAP algorithm](https://github.com/YingfanWang/PaCMAP) support (model off of `scanpy.tl.umap()`)
- [ ] add [LocalMAP algorithm](https://github.com/williamsyy/LocalMAP) support (model off of `scanpy.tl.umap()`)
    - despite the claims of PaCMAP's README, LocalMAP is already merged into PaCMAP's codebase
- [ ] val.viz.schematic_diagram: remove misleading grid lines from all non-X blocks
- [ ] val.datasets.polis.translate_statements: make `is_translated` more clever (to know when things AREN'T translated)
- [ ] val.datasets.polis: summarize changes during .load()
- [x] val.viz.schematic_diagram: when no diff_from arg provided, assume diff from empty AnnData ([#4](https://github.com/patcon/valency-anndata/pull/4))
- [ ] val.datasets.polis: when above n participants, require allow_large_scrape=True to use convo url
- [ ] val.datasets.polis: rename adata.uns objects to have `raw_` prefix
- [ ] val.tools.recipe_polis: document polis pipeline, e.g. flowchart
- [ ] val.tools.impute: convert to use scikit-learn's SimpleImputer
- [ ] val.tools.polis: Add helper method to extract metadata statements into obs
    - clustering doesn't happen on these, as these votes get zeroed out in `zero_mask`
- [ ] val.tools.kmeans: get kmeans++ init strategy working (in red-dwarf)


See [`TODO.md`](./TODO.md) or ask patcon for other ideas!

## B. Learning: Context First
- [ ] Watch a video about the internals of Pol.is
    - Video: [Pol.is Explainer](https://www.youtube.com/watch?v=FrIin_omVn4) (40 min)
- [ ] Watch this video about polislike tools from Hackers on Planet Earth (HOPE) 2024 conference
    - Video: [A Revolution in Representation: Computation Comes to Democracy's Aid](https://www.youtube.com/watch?v=_2DXwPDcJ1U) (55 min)
- [ ] Review a flowchart of the Polis data processing pipeline
    - https://github.com/polis-community/red-dwarf/issues/53#issuecomment-2988756562
- [ ] Perform the basic tutorial for scanpy, the tool that inspires this project's approach.
    - Scanpy tutorial: [Preprocessing and clustering](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html)
    - Use [Google CoLab](https://colab.research.google.com/)'s Notebook environment to run through it in the browser
- [ ] Open and play around with valency-anndata's own Python Notebook
    - https://colab.research.google.com/github/patcon/valency-anndata/blob/main/example-usage.ipynb
- [ ] Review this notebook that migrates large-scale representative likert survey data into polislike form, for side-by-side analysis
    - https://github.com/nishio/UTAS-UMAP/
- [ ] Learn about matrix completion, which we'll eventually need to explore to help show the statements to people in the most helpful order for placing them in the perspective landscape
    - Wikipedia: [matrix completion](https://en.wikipedia.org/wiki/Matrix_completion)
    - 
