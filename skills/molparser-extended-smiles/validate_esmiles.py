#!/usr/bin/env python3
"""Lightweight validator for MolParser E-SMILES (current project scope).

Supported extension records:
  - <a>[ATOM_INDEX]:[GROUP_NAME]</a>
  - <d>[ATOM_INDEX]:<dum></d>
  - <r>[RING_INDEX]:[GROUP_NAME]</r>
  - <r><v>[VIRTUALARC_INDEX]:[GROUP_NAME]</r>
  - <c>[ATOM_INDEX]:[RING_LABEL]</c>
  - <s>[SUBSTRUCTURE_ESMILES]</s>
  - <g>[INNER_PORT:OUTER_PORT]:...:|Sg:n|</g>
  - <v>[VIRTUALARC_INDEX]:[VIRTUALARC_NAME]:[FROM_ATOM:TO_ATOM]</v>
  - |Sg:n| (structural repeating unit marker)

Notes:
  - This script validates notation shape and token structure, not full chemistry.
  - The base SMILES before <sep> must be parseable by RDKit.
  - It does not perform valence, aromaticity, stereochemical, or reaction-mechanism checks.
  - Warnings are aligned with the current molparser.utils translator/drawer parsing scope.
"""

from __future__ import annotations

import argparse
import importlib
import re
import sys
from dataclasses import dataclass


RECORD_RE = re.compile(r"<(?P<tag>a|d|r|c)>(?P<body>.*?)</(?P=tag)>", re.DOTALL)
SUBSTRUCT_RE = re.compile(r"<s>(?P<body>.*?)</s>", re.DOTALL)
SGROUP_RE = re.compile(r"<g>(?P<body>.*?)</g>", re.DOTALL)
VIRTUAL_RE = re.compile(r"<v>(?P<body>.*?)</v>", re.DOTALL)
SG_RE = re.compile(r"\|Sg:(?P<count>[^|]+)\|")
INDEX_VALUE_RE = re.compile(r"^(?P<index>\d+):(?P<value>.+)$", re.DOTALL)
RING_VIRTUAL_C_RE = re.compile(r"^<c>(?P<index>\d+):(?P<value>.+)$", re.DOTALL)
RING_VIRTUAL_V_RE = re.compile(r"^<v>(?P<index>\d+):(?P<value>.+)$", re.DOTALL)
VIRTUAL_ARC_RE = re.compile(
    r"^(?P<index>\d+):(?P<name>[^:]+):\[(?P<from>\d+):(?P<to>\d+)\]$",
    re.DOTALL,
)
PORT_RE = re.compile(r"\[(?P<inner>\d+):(?P<outer>\d+)\]")
SG_COUNT_RE = re.compile(r"^[A-Za-z0-9]+(?:-[A-Za-z0-9]+)?$")
REPEAT_SUFFIX_RE = re.compile(r"\?(?:[a-z]|\d+|\d+-\d+)$")
SPECIAL_ID_RE = re.compile(r"^<id>\[[^\]]+\]$")


@dataclass
class Message:
    level: str
    text: str


def add(messages: list[Message], level: str, text: str) -> None:
    messages.append(Message(level=level, text=text))


def _merged_spans(spans: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(spans):
        if not merged or start > merged[-1][1]:
            merged.append((start, end))
        else:
            prev_start, prev_end = merged[-1]
            merged[-1] = (prev_start, max(prev_end, end))
    return merged


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
        # Advanced source-level form accepted by current molparser.utils parsing:
        # <r><c>[INDEX]:[VALUE]</r>
        match = RING_VIRTUAL_C_RE.match(stripped)
        if match:
            add(messages, "warning", "<r><c>... virtual-ring form detected; treated as advanced reference syntax")
        else:
            match = RING_VIRTUAL_V_RE.match(stripped)
            if match:
                add(messages, "warning", "<r><v>... virtualArc substituent detected; currently preserved as pre-compatible syntax")
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

    # Current molparser.utils parsing does not preserve group names containing spaces.
    if re.search(r"\s", value):
        add(messages, "error", f"<{tag}> value contains whitespace and may be dropped by parser: {value!r}")

    if "<sep>" in value:
        add(messages, "warning", f"<{tag}> value contains <sep>; check for accidental nesting")

    if tag == "d":
        if value != "<dum>":
            add(messages, "error", f"<d> is reserved for dummy attachment points; expected <dum>, got: {value!r}")
        return

    if tag == "a" and value == "<dum>":
        return

    if value.startswith("<id>"):
        if not SPECIAL_ID_RE.match(value):
            add(messages, "error", f"<id> special Markush label should use <id>[NOTE], got: {value!r}")
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
                    f"<{tag}> repeat suffix may not be parsed by current molparser.utils; use ?n / ?1-3 / ?3 / ?10: {suffix!r}",
                )


def _validate_sgroup_count(count: str, messages: list[Message], context: str) -> None:
    clean = count.strip()
    if not clean:
        add(messages, "error", f"{context} count is empty")
    elif not SG_COUNT_RE.match(clean):
        add(messages, "warning", f"{context} count is syntactically unusual: {clean!r}")


def _validate_substruct(body: str, messages: list[Message]) -> None:
    if "<sep>" not in body:
        add(messages, "error", "<s> substructure should contain SMILES<sep>EXTENSION")
        return
    add(messages, "warning", "<s> substructure syntax detected; currently preserved without semantic expansion")


def _validate_sgroup(body: str, messages: list[Message]) -> None:
    sg_match = SG_RE.search(body)
    if not sg_match:
        add(messages, "error", "<g> should include a |Sg:...| repeat marker")
        return
    _validate_sgroup_count(sg_match.group("count"), messages, "<g> |Sg:...|")

    prefix = body[:sg_match.start()]
    suffix = body[sg_match.end():].strip()
    if suffix:
        add(messages, "warning", f"<g> has trailing content after |Sg:...|: {suffix!r}")

    ports = list(PORT_RE.finditer(prefix))
    if not ports:
        add(messages, "error", "<g> should include at least one [INNER_PORT:OUTER_PORT] pair")
        return
    leftover = PORT_RE.sub("", prefix).replace(":", "").strip()
    if leftover:
        add(messages, "error", f"<g> contains unparsed port text: {leftover!r}")
    add(messages, "warning", "<g> s-group repeat syntax detected; currently preserved without semantic expansion")


def _validate_virtual_arc(body: str, messages: list[Message]) -> None:
    match = VIRTUAL_ARC_RE.match(body.strip())
    if not match:
        add(messages, "error", f"<v> should use [INDEX]:[NAME]:[FROM_ATOM:TO_ATOM], got: {body!r}")
        return
    if not match.group("name").strip():
        add(messages, "error", "<v> virtualArc name is empty")
    add(messages, "warning", "<v> virtualArc syntax detected; endpoints are treated as unordered and currently preserved")


def validate(esmiles: str, strict: bool = False) -> list[Message]:
    messages: list[Message] = []
    text = str(esmiles).strip()
    if not text:
        add(messages, "error", "input is empty")
        return messages

    if "<sep>" not in text:
        add(messages, "error", "top-level <sep> count should be at least 1, got 0")
        return messages

    base, extension = text.split("<sep>", 1)
    if not base.strip():
        add(messages, "error", "base SMILES is empty")
    elif not _is_rdkit_smiles(base.strip()):
        add(messages, "error", f"base SMILES is not parseable by RDKit: {base!r}")

    if extension:
        for tag in ("a", "d", "r", "c", "s", "g", "v"):
            opens = len(re.findall(fr"<{tag}>", extension))
            if tag == "c":
                # In <r><c>[INDEX]:[VALUE]</r>, <c> is an inline target marker
                # and is not expected to have a closing </c>.
                virtual_c_opens = len(re.findall(r"<r>\s*<c>\s*\d+:", extension))
                opens -= virtual_c_opens
                opens = max(opens, 0)
            if tag == "v":
                # In <r><v>[INDEX]:[VALUE]</r>, <v> is an inline virtualArc marker
                # and is not expected to have a closing </v>.
                virtual_v_opens = len(re.findall(r"<r>\s*<v>\s*\d+:", extension))
                opens -= virtual_v_opens
                opens = max(opens, 0)
            closes = len(re.findall(fr"</{tag}>", extension))
            if opens != closes:
                add(messages, "error", f"<{tag}> opening/closing tag mismatch: {opens} != {closes}")

    spans: list[tuple[int, int]] = []
    for match in RECORD_RE.finditer(extension):
        spans.append(match.span())
        _validate_record(match.group("tag"), match.group("body"), messages)

    for match in SUBSTRUCT_RE.finditer(extension):
        spans.append(match.span())
        _validate_substruct(match.group("body"), messages)

    for match in SGROUP_RE.finditer(extension):
        spans.append(match.span())
        _validate_sgroup(match.group("body"), messages)

    for match in VIRTUAL_RE.finditer(extension):
        spans.append(match.span())
        _validate_virtual_arc(match.group("body"), messages)

    for sg_match in SG_RE.finditer(extension):
        count = sg_match.group("count").strip()
        _validate_sgroup_count(count, messages, "|Sg:...|")
        if count != "n":
            add(messages, "warning", f"|Sg:{count}| is accepted, but current molparser.utils only flags |Sg:n| as SRU")

    leftovers = extension
    for start, end in sorted(_merged_spans(spans), reverse=True):
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
