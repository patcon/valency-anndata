from typing import Literal, Optional
import valency_anndata as val


JapanChoiceTopic = Literal[
    "2025_foreign_affairs_security",
    "2025_diversity_human_rights",
    "2025_education_children_old_age",
    "2025_economy_taxation_employment",
    "2026_foreign_affairs_security",
    "2026_diversity_human_rights",
    "2026_education_children_old_age",
    "2026_economy_taxation_employment",
]

_TOPIC_SOURCES = {
    "2025_foreign_affairs_security": "huggingface:patcon/japanchoice-2025-foreign-policy-and-security",
    "2025_economy_taxation_employment": "huggingface:patcon/japanchoice-2025-economy-taxation-employment",
    "2026_foreign_affairs_security": "huggingface:patcon/japanchoice-2026-foreign-policy-and-security",
    "2026_diversity_human_rights": "huggingface:patcon/japanchoice-2026-diversity-and-human-rights",
    "2026_education_children_old_age": "huggingface:patcon/japanchoice-2026-education-children-old-age",
    "2026_economy_taxation_employment": "huggingface:patcon/japanchoice-2026-economy-taxation-employment",
}

# Topics that were once available at pol.is but have since been taken down,
# with no HuggingFace archive to fall back on.
_UNAVAILABLE_TOPICS = {
    "2025_diversity_human_rights",
    "2025_education_children_old_age",
}


def japanchoice(
    topic: JapanChoiceTopic | str,
    translate_to: Optional[str] = None,
    **kwargs,
):
    """
    Polis conversations from Japan Choice, a Japanese civic engagement platform.

    Japan Choice runs Polis conversations on key policy topics ahead of Japanese
    elections, allowing citizens to share and compare their views on national issues.
    Conversations are in Japanese.

    Data is archived as CSV exports on HuggingFace, under
    <https://huggingface.co/patcon>.

    See: <https://japanchoice.jp/polis>

    Parameters
    ----------
    topic : str
        The policy topic and year to load. One of:

        - ``"2025_foreign_affairs_security"`` — Foreign Affairs & Security (2025)
        - ``"2025_diversity_human_rights"`` — Diversity & Human Rights (2025) — **no longer available**
        - ``"2025_education_children_old_age"`` — Education, Children & Old Age Care (2025) — **no longer available**
        - ``"2025_economy_taxation_employment"`` — Economy, Taxation & Employment (2025)
        - ``"2026_foreign_affairs_security"`` — Foreign Affairs & Security (2026)
        - ``"2026_diversity_human_rights"`` — Diversity & Human Rights (2026)
        - ``"2026_education_children_old_age"`` — Education, Children & Old Age Care (2026)
        - ``"2026_economy_taxation_employment"`` — Economy, Taxation & Employment (2026)

        The topics marked "no longer available" were originally hosted as
        live pol.is conversations but have since been taken down, and no
        HuggingFace archive currently exists for them. Requesting one of
        these raises a ``ValueError``.

    translate_to : str or None, optional
        Target language code (e.g., ``"en"``, ``"fr"``) for translating
        statement text. Defaults to None (no translation).

    Returns
    -------
    adata : anndata.AnnData
        AnnData object containing the loaded Polis conversation.

    Examples
    --------
    Load the 2025 Economy, Taxation & Employment conversation:

    ```py
    adata = val.datasets.japanchoice("2025_economy_taxation_employment")
    ```

    Load the 2026 Foreign Affairs & Security conversation translated to English:

    ```py
    adata = val.datasets.japanchoice("2026_foreign_affairs_security", translate_to="en")
    ```

    Attribution
    -----------

    Data was gathered using the Polis software (see:
    <https://compdemocracy.org/polis> and
    <https://github.com/compdemocracy/polis>) and is sub-licensed under CC BY
    4.0 with Attribution to The Computational Democracy Project.
    """
    if topic in _UNAVAILABLE_TOPICS:
        raise ValueError(
            f"The {topic!r} dataset is no longer available. It was "
            "originally hosted as a live pol.is conversation, which has "
            "since been taken down, and no HuggingFace archive currently "
            "exists for it."
        )

    if topic not in _TOPIC_SOURCES:
        raise ValueError(
            f"Unknown topic {topic!r}. Must be one of: {list(JapanChoiceTopic.__args__)}"
        )
    source = _TOPIC_SOURCES[topic]

    adata = val.datasets.polis.load(source, translate_to=translate_to, **kwargs)

    return adata
