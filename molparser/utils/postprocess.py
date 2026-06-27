"""Post-processing for MolParser model outputs."""

from __future__ import annotations

from typing import Any, Dict, Sequence

try:
    from .translator import Translator
except ImportError:  # Support running from package directory as working directory.
    from translator import Translator


def postprocess_caption(caption: str, error_msg: bool = False) -> Dict[str, object]:
    """Refactor a raw caption into normalized SMILES, E-SMILES, and CXSMILES."""
    raw_caption = str(caption).strip()
    result = Translator.refactor(raw_caption, error_msg=error_msg)
    cxsmiles = Translator.esmiles_to_cxsmiles(raw_caption, error_msg=error_msg)
    if result is None:
        raw_smi = raw_caption.split("<sep>", 1)[0]
        raw_groups = raw_caption.split("<sep>", 1)[1] if "<sep>" in raw_caption else ""
        return {
            "caption": raw_caption,
            "smi": raw_smi,
            "esmi": raw_caption if "<sep>" in raw_caption else f"{raw_smi}<sep>",
            "cxsmiles": cxsmiles,
            "markush": "<sep>" in raw_caption and raw_groups != "",
            "sru": False,
            "groups": raw_groups,
        }
    return {
        "caption": result.caption,
        "smi": result.smi,
        "esmi": result.esmi,
        "cxsmiles": cxsmiles,
        "markush": result.markush,
        "sru": result.sru,
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
