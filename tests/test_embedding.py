import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest
from anndata import AnnData
import scanpy as sc

import valency_anndata as val
from valency_anndata.viz._embedding import (
    _expand_color_specs,
    _parse_color_spec,
)


def _embedding_adata():
    adata = AnnData(np.zeros((4, 2)))
    adata.obs_names = [f"cell_{i}" for i in range(4)]
    adata.obs["cluster"] = ["a", "b", "a", "b"]
    adata.obs["kmeans_pacmap"] = ["0", "1", "0", "1"]
    adata.obsm["X_pca"] = np.array(
        [
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 10.0, 11.0, 12.0],
            [13.0, 14.0, 15.0, 16.0],
        ]
    )
    adata.obsm["evoc_polis2"] = np.array(
        [
            [10.0, 20.0],
            [30.0, 40.0],
            [50.0, 60.0],
            [70.0, 80.0],
        ]
    )
    adata.obsm["X_pacmap"] = np.array(
        [
            [0.0, 0.0],
            [1.0, 1.0],
            [2.0, 2.0],
            [3.0, 3.0],
        ]
    )
    return adata


class TestParseColorSpec:
    def test_single_index_parse(self):
        assert _parse_color_spec("X_pca[2]") == ("X_pca", 2, None)

    def test_non_x_key_parse(self):
        assert _parse_color_spec("evoc_polis2[1]") == ("evoc_polis2", 1, None)

    def test_plain_key_returns_none(self):
        assert _parse_color_spec("cluster") is None


class TestExpandColorSpecs:
    def test_bounded_slice_expansion(self):
        assert _expand_color_specs("X_pca[2:6]") == [
            "X_pca[2]",
            "X_pca[3]",
            "X_pca[4]",
            "X_pca[5]",
        ]

    def test_prefix_slice_expansion(self):
        assert _expand_color_specs("X_pca[:2]") == [
            "X_pca[0]",
            "X_pca[1]",
        ]

    def test_non_x_single_index_stays_single(self):
        assert _expand_color_specs("evoc_polis[0]") == "evoc_polis[0]"

    def test_plain_key_stays_plain(self):
        assert _expand_color_specs("cluster") == "cluster"

    def test_non_x_slice_expansion(self):
        assert _expand_color_specs("evoc_polis2[1:3]") == [
            "evoc_polis2[1]",
            "evoc_polis2[2]",
        ]

    def test_mixed_list_preserves_plain_obs_keys(self):
        assert _expand_color_specs(["kmeans_pacmap", "X_pca[2:4]", "pct_seen"]) == [
            "kmeans_pacmap",
            "X_pca[2]",
            "X_pca[3]",
            "pct_seen",
        ]

    def test_mixed_list_with_non_x_key(self):
        assert _expand_color_specs(["cluster", "evoc_polis2[0]"]) == [
            "cluster",
            "evoc_polis2[0]",
        ]

    def test_malformed_syntax_rejected(self):
        with pytest.raises(ValueError, match="Invalid embedding color spec"):
            _expand_color_specs("X_pca[nope]")

    @pytest.mark.parametrize(
        "color",
        [
            "X_pca[2:]",
            "X_pca[:]",
            "X_pca[2:2]",
            "X_pca[6:2]",
        ],
    )
    def test_invalid_slices_rejected(self, color):
        with pytest.raises(ValueError, match="Invalid embedding color slice"):
            _expand_color_specs(color)


class TestEmbeddingWrapper:
    def test_normal_color_key_forwards_unchanged(self, monkeypatch):
        adata = _embedding_adata()
        captured = {}

        def fake_embedding(adata, *args, color=None, **kwargs):
            captured["adata"] = adata
            captured["args"] = args
            captured["color"] = color
            captured["kwargs"] = kwargs
            return "ok"

        monkeypatch.setattr(sc.pl, "embedding", fake_embedding)

        result = val.viz.embedding(adata, basis="pacmap", color="cluster", show=False)

        assert result == "ok"
        assert captured["adata"] is adata
        assert captured["color"] == "cluster"
        assert captured["kwargs"]["basis"] == "pacmap"

    def test_single_obsm_color_gets_rewritten(self, monkeypatch):
        adata = _embedding_adata()
        captured = {}

        def fake_embedding(adata, *args, color=None, **kwargs):
            captured["adata"] = adata
            captured["color"] = color
            return "ok"

        monkeypatch.setattr(sc.pl, "embedding", fake_embedding)

        val.viz.embedding(adata, basis="pacmap", color="X_pca[2]", show=False)

        assert captured["adata"] is not adata
        assert captured["color"] == "X_pca[2]"
        np.testing.assert_array_equal(
            captured["adata"].obs["X_pca[2]"].to_numpy(),
            adata.obsm["X_pca"][:, 2],
        )

    def test_slice_color_gets_expanded_and_rewritten(self, monkeypatch):
        adata = _embedding_adata()
        captured = {}

        def fake_embedding(adata, *args, color=None, **kwargs):
            captured["adata"] = adata
            captured["color"] = color
            return "ok"

        monkeypatch.setattr(sc.pl, "embedding", fake_embedding)

        val.viz.embedding(adata, basis="pacmap", color="X_pca[2:4]", show=False)

        assert captured["color"] == ["X_pca[2]", "X_pca[3]"]
        np.testing.assert_array_equal(
            captured["adata"].obs["X_pca[2]"].to_numpy(),
            adata.obsm["X_pca"][:, 2],
        )
        np.testing.assert_array_equal(
            captured["adata"].obs["X_pca[3]"].to_numpy(),
            adata.obsm["X_pca"][:, 3],
        )

    def test_mixed_plain_and_obsm_colors_work(self, monkeypatch):
        adata = _embedding_adata()
        captured = {}

        def fake_embedding(adata, *args, color=None, **kwargs):
            captured["adata"] = adata
            captured["color"] = color
            return "ok"

        monkeypatch.setattr(sc.pl, "embedding", fake_embedding)

        val.viz.embedding(adata, basis="pacmap", color=["kmeans_pacmap", "X_pca[2]"], show=False)

        assert captured["color"] == ["kmeans_pacmap", "X_pca[2]"]
        np.testing.assert_array_equal(
            captured["adata"].obs["X_pca[2]"].to_numpy(),
            adata.obsm["X_pca"][:, 2],
        )

    def test_mixed_plain_and_non_x_obsm_colors_work(self, monkeypatch):
        adata = _embedding_adata()
        captured = {}

        def fake_embedding(adata, *args, color=None, **kwargs):
            captured["adata"] = adata
            captured["color"] = color
            return "ok"

        monkeypatch.setattr(sc.pl, "embedding", fake_embedding)

        val.viz.embedding(adata, basis="pacmap", color=["cluster", "evoc_polis2[1]"], show=False)

        assert captured["color"] == ["cluster", "evoc_polis2[1]"]
        np.testing.assert_array_equal(
            captured["adata"].obs["evoc_polis2[1]"].to_numpy(),
            adata.obsm["evoc_polis2"][:, 1],
        )

    def test_original_obs_not_mutated(self, monkeypatch):
        adata = _embedding_adata()
        original_columns = list(adata.obs.columns)

        def fake_embedding(adata, *args, color=None, **kwargs):
            return "ok"

        monkeypatch.setattr(sc.pl, "embedding", fake_embedding)

        val.viz.embedding(adata, basis="pacmap", color="X_pca[2]", show=False)

        assert list(adata.obs.columns) == original_columns
        assert "X_pca[2]" not in adata.obs.columns

    def test_missing_obsm_key_raises_clear_error(self):
        adata = _embedding_adata()

        with pytest.raises(KeyError, match="references missing adata.obsm\\['missing'\\]"):
            val.viz.embedding(adata, basis="pacmap", color="missing[0]", show=False)

    def test_out_of_range_index_raises_clear_error(self):
        adata = _embedding_adata()

        with pytest.raises(IndexError, match="out of bounds"):
            val.viz.embedding(adata, basis="pacmap", color="X_pca[10]", show=False)


class TestEmbeddingSmoke:
    def test_real_plotting_path_runs(self):
        import matplotlib.pyplot as plt

        adata = _embedding_adata()
        result = val.viz.embedding(
            adata,
            basis="pacmap",
            color=["cluster", "X_pca[1]"],
            show=False,
        )

        assert result is not None
        plt.close("all")
