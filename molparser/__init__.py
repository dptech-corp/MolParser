"""Public package namespace for MolParser."""

from __future__ import annotations

import sys
from importlib import import_module

__version__ = "0.1.0"

utils = import_module("utils")
sys.modules[__name__ + ".utils"] = utils

__all__ = ["__version__", "utils"]
