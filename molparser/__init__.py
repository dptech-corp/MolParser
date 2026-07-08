from __future__ import annotations

from . import utils

__version__ = "0.1.0"

__all__ = ["MolParser", "MolParserConfig", "MolParserResult", "__version__", "utils"]


def __getattr__(name: str):
    if name in {"MolParser", "MolParserConfig", "MolParserResult"}:
        from .models import MolParser, MolParserConfig, MolParserResult

        return {"MolParser": MolParser, "MolParserConfig": MolParserConfig, "MolParserResult": MolParserResult}[name]
    raise AttributeError(name)
