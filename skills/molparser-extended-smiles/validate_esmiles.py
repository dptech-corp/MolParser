#!/usr/bin/env python3
"""Lightweight validator for MolParser E-SMILES (current project scope).

Supported extension records:
  - <a>[ATOM_INDEX]:[GROUP_NAME]</a>
  - <r>[RING_INDEX]:[GROUP_NAME]</r>
  - <c>[ATOM_INDEX]:[RING_LABEL]</c>
  - |Sg:n| (structural repeating unit marker)

Notes:
  - This script validates notation shape and token structure, not full chemistry.
  - The base SMILES before <sep> must be parseable by RDKit.
  - It does not perform valence, aromaticity, stereochemical, or reaction-mechanism checks.
  - Warnings are aligned with the current utils translator/drawer parsing scope.
"""

from __future__ import annotations

import argparse
import importlib
import re
import sys
from dataclasses import dataclass


RECORD_RE = re.compile(r"<(?P<tag>a|r|c)>(?P<body>.*?)</(?P=tag)>", re.DOTALL)
SG_RE = re.compile(r"\|Sg:(?P<count>[^|]+)\|")
INDEX_VALUE_RE = re.compile(r"^(?P<index>\d+):(?P<value>.+)$", re.DOTALL)
RING_VIRTUAL_C_RE = re.compile(r"^<c>(?P<index>\d+):(?P<value>.+)$", re.DOTALL)
REPEAT_SUFFIX_RE = re.compile(r"\?(?:[a-z]|\d|\d-\d)$")


@dataclass
class Message:
    level: str
    text: str


def add(messages: list[Message], level: str, text: str) -> None:
    messages.append(Message(level=level, text=text))


def _is_rdkit_smiles(smiles: str) -> bool:
    chem = importlib.import_module("rdkit.Chem")
    rd_logger = importlib.import_module("rdkit.RDLogger")

    rd_logger.DisableLog("rdApp.*")
    try:
        return chem.MolFromSmiles(smiles) is not None
    finally:
        rd_logger.EnableLog("rdApp.*")


def _validate_record(tag: str, body: str, messages: list[Message]) -> None:
    stripped = body.strip()
    match = INDEX_VALUE_RE.match(stripped)
    if not match and tag == "r":
        # Advanced source-level form accepted by current utils parsing:
        # <r><c>[INDEX]:[VALUE]</r>
        match = RING_VIRTUAL_C_RE.match(stripped)
        if match:
            add(messages, "warning", "<r><c>... virtual-ring form detected; treated as advanced reference syntax")
    if not match:
        add(messages, "error", f"<{tag}> should use [INDEX]:[VALUE], got: {body!r}")
        return

    index = match.group("index")
    value = match.group("value").strip()
    if not index.isdigit():
        add(messages, "error", f"<{tag}> index is not a non-negative integer: {index!r}")
    if not value:
        add(messages, "error", f"<{tag}> value is empty")
        return

    # Current utils parsing does not preserve group names containing spaces.
    if re.search(r"\s", value):
        add(messages, "error", f"<{tag}> value contains whitespace and may be dropped by parser: {value!r}")

    if "<sep>" in value:
        add(messages, "warning", f"<{tag}> value contains <sep>; check for accidental nesting")

    if tag == "a" and value == "<dum>":
        return

    if "? " in value:
        add(messages, "warning", f"<{tag}> repeat suffix contains whitespace: {value!r}")

    if "?" in value:
        parts = value.rsplit("?", 1)
        if len(parts) == 2 and parts[1]:
            suffix = "?" + parts[1]
            if not REPEAT_SUFFIX_RE.match(suffix):
                add(
                    messages,
                    "warning",
                    f"<{tag}> repeat suffix may not be parsed by current utils; use ?n / ?1-3 / ?3: {suffix!r}",
                )


def validate(esmiles: str, strict: bool = False) -> list[Message]:
    messages: list[Message] = []
    text = str(esmiles).strip()
    if not text:
        add(messages, "error", "input is empty")
        return messages

    sep_count = text.count("<sep>")
    if sep_count != 1:
        add(messages, "error", f"top-level <sep> count should be 1, got {sep_count}")
        return messages

    base, extension = text.split("<sep>", 1)
    if not base.strip():
        add(messages, "error", "base SMILES is empty")
    elif not _is_rdkit_smiles(base.strip()):
        add(messages, "error", f"base SMILES is not parseable by RDKit: {base!r}")

    if extension:
        for tag in ("a", "r", "c"):
            opens = len(re.findall(fr"<{tag}>", extension))
            if tag == "c":
                # In <r><c>[INDEX]:[VALUE]</r>, <c> is an inline target marker
                # and is not expected to have a closing </c>.
                virtual_c_opens = len(re.findall(r"<r>\s*<c>\s*\d+:", extension))
                opens -= virtual_c_opens
                opens = max(opens, 0)
            closes = len(re.findall(fr"</{tag}>", extension))
            if opens != closes:
                add(messages, "error", f"<{tag}> opening/closing tag mismatch: {opens} != {closes}")

    spans: list[tuple[int, int]] = []
    for match in RECORD_RE.finditer(extension):
        spans.append(match.span())
        _validate_record(match.group("tag"), match.group("body"), messages)

    for sg_match in SG_RE.finditer(extension):
        count = sg_match.group("count").strip()
        if not count:
            add(messages, "error", "|Sg:n| count is empty")
        elif count != "n":
            add(
                messages,
                "warning",
                f"|Sg:{count}| is syntactically accepted here, but current utils only flags |Sg:n| as SRU",
            )

    leftovers = extension
    for start, end in sorted(spans, reverse=True):
        leftovers = leftovers[:start] + leftovers[end:]
    leftovers = SG_RE.sub("", leftovers)
    leftovers = leftovers.strip()
    if leftovers:
        level = "error" if strict else "warning"
        add(messages, level, f"unparsed extension text remains: {leftovers!r}")

    if not any(item.level == "error" for item in messages):
        add(messages, "ok", "basic E-SMILES syntax checks passed")
    return messages


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Validate MolParser E-SMILES basic syntax.")
    parser.add_argument("esmiles", nargs="?", help="E-SMILES string. If omitted, read from stdin.")
    parser.add_argument("--strict", action="store_true", help="Treat unparsed extension text as error.")
    args = parser.parse_args(argv)

    text = args.esmiles if args.esmiles is not None else sys.stdin.read()
    messages = validate(text, strict=args.strict)
    for message in messages:
        print(f"[{message.level}] {message.text}")
    return 1 if any(msg.level == "error" for msg in messages) else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
