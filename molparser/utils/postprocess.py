"""Post-processing for MolParser model outputs."""

from __future__ import annotations

import logging
import re
from typing import Any, Dict, Sequence

try:
    from .translator import Translator
except ImportError:  # Support running from package directory as working directory.
    from translator import Translator


logger = logging.getLogger(__name__)


_CONCRETE_REPEAT_PATTERN = re.compile(
    r"\?[1-9]\d*(?=</(?:a|r)>)|\|Sg:[1-9]\d*\|"
)


def _physicalize_concrete_repeats(caption: str, error_msg: bool) -> str:
    """Expand one deterministic fixed repeat without choosing an isomer."""
    if _CONCRETE_REPEAT_PATTERN.search(caption) is None:
        return caption
    try:
        expanded = Translator.substitute_markush(
            caption,
            {},
            error_msg=error_msg,
            repeat_policy="best_effort",
        )
    except ValueError as exc:
        if error_msg:
            logger.warning("Concrete repeat expansion was preserved: %s", exc)
        return caption
    return expanded if isinstance(expanded, str) else caption


def _has_top_level_sru(caption: str) -> bool:
    if "<sep>" not in caption:
        return False
    trailing = caption.split("<sep>", 1)[1]
    _, ext = Translator.parse_trailing(trailing)
    return bool(
        re.fullmatch(
            r"Sg:[A-Za-z0-9]+(?:-[A-Za-z0-9]+)?",
            Translator.parse_extension(ext),
        )
    )


def postprocess_caption(caption: str, error_msg: bool = False) -> Dict[str, object]:
    """Refactor a raw caption into normalized SMILES, E-SMILES, and CXSMILES."""
    raw_caption = str(caption).strip()
    source_sru = _has_top_level_sru(raw_caption)
    effective_caption = _physicalize_concrete_repeats(raw_caption, error_msg)
    result = Translator.refactor(effective_caption, error_msg=error_msg)
    normalized_caption = result.esmi if result is not None else effective_caption
    cxsmiles = Translator.esmiles_to_cxsmiles(
        normalized_caption,
        error_msg=error_msg,
    )
    if result is None:
        raw_smi = effective_caption.split("<sep>", 1)[0]
        raw_groups = (
            effective_caption.split("<sep>", 1)[1]
            if "<sep>" in effective_caption
            else ""
        )
        return {
            "caption": raw_caption,
            "smi": raw_smi,
            "esmi": (
                effective_caption
                if "<sep>" in effective_caption
                else f"{raw_smi}<sep>"
            ),
            "cxsmiles": cxsmiles,
            "markush": "<sep>" in raw_caption and raw_groups != "",
            "sru": source_sru,
            "groups": raw_groups,
        }
    return {
        "caption": raw_caption,
        "smi": result.smi,
        "esmi": result.esmi,
        "cxsmiles": cxsmiles,
        "markush": result.markush,
        "sru": source_sru or result.sru,
        "groups": result.groups,
    }


def extract_confidence(sequence: Sequence[int], scores: Sequence[Any]) -> float:
    """Take the lowest token-probability over the generated sequence."""
    if not scores:
        return 0.0
    generated_length = len(scores)
    chosen_tokens = sequence[-generated_length:]
    min_prob = 1.0
    for step_scores, token_id in zip(scores, chosen_tokens):
        probs = step_scores.softmax(dim=-1)
        prob = probs[int(token_id)].item()
        min_prob = min(min_prob, prob)
    return float(min_prob)


__all__ = ["postprocess_caption", "extract_confidence"]
