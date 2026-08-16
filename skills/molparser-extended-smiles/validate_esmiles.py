#!/usr/bin/env python3
"""Lightweight validator for MolParser E-SMILES (current project scope).

Supported extension records:
  - <a>[ATOM_INDEX]:[GROUP_NAME]</a>
  - <a>[ATOM_INDEX]:<id>[NOTE]</a>
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
RING_VIRTUAL_V_RE = re.compile(
    r"^<v>\s*(?P<index>\d+):(?P<value>.+)$",
    re.DOTALL,
)
VIRTUAL_ARC_RE = re.compile(
    r"^(?P<index>\d+):(?P<name>[^:]*):\[(?P<from>\d+):(?P<to>\d+)\]$",
    re.DOTALL,
)
PORT_RE = re.compile(r"\[(?P<inner>\d+):(?P<outer>\d+)\]")
PORT_PREFIX_RE = re.compile(r"^(?:\[\d+:\d+\]:)+$")
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


def _mask_spans(text: str, spans: list[tuple[int, int]]) -> str:
    """Blank spans while preserving offsets for top-level record matching."""
    chars = list(text)
    for start, end in _merged_spans(spans):
        chars[start:end] = " " * (end - start)
    return "".join(chars)


def _validate_record(
    tag: str,
    body: str,
    messages: list[Message],
    *,
    atom_count: int | None = None,
    ring_count: int | None = None,
) -> tuple[str, int] | None:
    stripped = body.strip()
    match = INDEX_VALUE_RE.match(stripped)
    namespace = "ring" if tag == "r" else "atom"
    if not match and tag == "r":
        # Advanced source-level form accepted by current molparser.utils parsing:
        # <r><c>[INDEX]:[VALUE]</r>
        match = RING_VIRTUAL_C_RE.match(stripped)
        if match:
            namespace = "virtual_ring"
            add(messages, "warning", "<r><c>... virtual-ring form detected; treated as advanced reference syntax")
        else:
            match = RING_VIRTUAL_V_RE.match(stripped)
            if match:
                namespace = "virtual_arc"
    if not match:
        add(messages, "error", f"<{tag}> should use [INDEX]:[VALUE], got: {body!r}")
        return None

    index = match.group("index")
    value = match.group("value").strip()
    if not index.isdigit():
        add(messages, "error", f"<{tag}> index is not a non-negative integer: {index!r}")
    if not value:
        add(messages, "error", f"<{tag}> value is empty")
        return None

    numeric_index = int(index)
    record_info = (
        (namespace, numeric_index) if namespace == "virtual_arc" else None
    )
    if namespace == "atom" and atom_count is not None and numeric_index >= atom_count:
        add(
            messages,
            "error",
            f"<{tag}> atom index {numeric_index} is outside the base molecule "
            f"(atom count {atom_count})",
        )
    elif namespace == "ring" and ring_count is not None and numeric_index >= ring_count:
        add(
            messages,
            "error",
            f"<r> ring index {numeric_index} is outside the base molecule "
            f"(ring count {ring_count})",
        )

    # Current molparser.utils parsing does not preserve group names containing spaces.
    if re.search(r"\s", value):
        add(messages, "error", f"<{tag}> value contains whitespace and may be dropped by parser: {value!r}")

    if "<sep>" in value:
        add(messages, "warning", f"<{tag}> value contains <sep>; check for accidental nesting")

    if tag == "d":
        if value != "<dum>":
            add(messages, "error", f"<d> is reserved for dummy attachment points; expected <dum>, got: {value!r}")
        return None

    if tag == "a" and value == "<dum>":
        return None

    if "<id>" in value:
        if tag != "a" or namespace != "atom":
            add(
                messages,
                "error",
                "<id>[NOTE] is supported only as an atom-indexed <a> payload",
            )
        if not SPECIAL_ID_RE.match(value):
            add(messages, "error", f"<id> special Markush label should use <id>[NOTE], got: {value!r}")
        return record_info

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
    return record_info


def _validate_sgroup_count(count: str, messages: list[Message], context: str) -> None:
    clean = count.strip()
    if not clean:
        add(messages, "error", f"{context} count is empty")
    elif not SG_COUNT_RE.match(clean):
        add(messages, "warning", f"{context} count is syntactically unusual: {clean!r}")
    elif clean.isdigit() and int(clean) < 1:
        add(messages, "error", f"{context} count must be positive")
    elif "-" in clean:
        start, end = clean.split("-", 1)
        if start.isdigit() and end.isdigit() and (
            int(start) < 1 or int(end) < int(start)
        ):
            add(messages, "error", f"{context} numeric range must be positive and ascending")


def _validate_substruct(body: str, messages: list[Message]) -> None:
    if "<sep>" not in body:
        add(messages, "error", "<s> substructure should contain SMILES<sep>EXTENSION")
        return
    add(
        messages,
        "warning",
        "<s> substructure syntax detected; labels and symbolic repeat counts can "
        "be substituted recursively, while the record remains a nested annotation",
    )


def _validate_sgroup(
    body: str,
    messages: list[Message],
    atom_count: int | None = None,
) -> None:
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
    if not PORT_PREFIX_RE.fullmatch(prefix.strip()):
        add(
            messages,
            "error",
            "<g> ports must use [INNER:OUTER]: records separated by colons",
        )
    if atom_count is not None:
        for port in ports:
            for role in ("inner", "outer"):
                index = int(port.group(role))
                if index >= atom_count:
                    add(
                        messages,
                        "error",
                        f"<g> {role} port index {index} is outside the base molecule "
                        f"(atom count {atom_count})",
                    )
    port_pairs = [
        (int(port.group("inner")), int(port.group("outer")))
        for port in ports
    ]
    if any(inner == outer for inner, outer in port_pairs):
        add(messages, "error", "<g> inner and outer indices in a port must differ")
    if len(set(port_pairs)) != len(port_pairs):
        add(messages, "error", "<g> boundary port pairs must be unique")


def _validate_virtual_arc(
    body: str,
    messages: list[Message],
    atom_count: int | None = None,
) -> re.Match[str] | None:
    match = VIRTUAL_ARC_RE.match(body.strip())
    if not match:
        add(messages, "error", f"<v> should use [INDEX]:[NAME]:[FROM_ATOM:TO_ATOM], got: {body!r}")
        return None
    if "<id>" in match.group("name"):
        add(
            messages,
            "error",
            "<id>[NOTE] is supported only as an atom-indexed <a> payload, not as a virtualArc name",
        )
    start = int(match.group("from"))
    end = int(match.group("to"))
    if start == end:
        add(messages, "error", "<v> virtualArc endpoints must be different atoms")
    if atom_count is not None and (start >= atom_count or end >= atom_count):
        add(
            messages,
            "error",
            f"<v> endpoint [{start}:{end}] is outside the base molecule "
            f"(atom count {atom_count})",
        )
    if start > end:
        add(
            messages,
            "warning",
            f"<v> endpoint order [{start}:{end}] is legacy-compatible; "
            f"new E-SMILES 2.0 data should use [{end}:{start}]",
        )
    return match


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
    mol = None
    if not base.strip():
        add(messages, "error", "base SMILES is empty")
    else:
        chem = importlib.import_module("rdkit.Chem")
        rd_logger = importlib.import_module("rdkit.RDLogger")
        rd_logger.DisableLog("rdApp.*")
        try:
            mol = chem.MolFromSmiles(base.strip())
        finally:
            rd_logger.EnableLog("rdApp.*")
        if mol is None:
            add(messages, "error", f"base SMILES is not parseable by RDKit: {base!r}")
    atom_count = mol.GetNumAtoms() if mol is not None else None
    ring_count = len(mol.GetRingInfo().AtomRings()) if mol is not None else None

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

    substruct_matches = list(SUBSTRUCT_RE.finditer(extension))
    substruct_spans = [match.span() for match in substruct_matches]
    spans: list[tuple[int, int]] = list(substruct_spans)
    top_extension = _mask_spans(extension, substruct_spans)

    for match in substruct_matches:
        body = match.group("body").strip()
        _validate_substruct(body, messages)
        if "<sep>" in body:
            for nested in validate(body, strict=strict):
                if nested.level != "ok":
                    add(messages, nested.level, f"<s> {nested.text}")

    virtual_ref_ids: list[int] = []
    for match in RECORD_RE.finditer(top_extension):
        spans.append(match.span())
        record_info = _validate_record(
            match.group("tag"),
            match.group("body"),
            messages,
            atom_count=atom_count,
            ring_count=ring_count,
        )
        if record_info is not None and record_info[0] == "virtual_arc":
            virtual_ref_ids.append(record_info[1])

    for match in SGROUP_RE.finditer(top_extension):
        spans.append(match.span())
        _validate_sgroup(match.group("body"), messages, atom_count=atom_count)

    virtual_arc_matches: list[re.Match[str]] = []
    for match in VIRTUAL_RE.finditer(top_extension):
        spans.append(match.span())
        parsed_arc = _validate_virtual_arc(
            match.group("body"),
            messages,
            atom_count=atom_count,
        )
        if parsed_arc is not None:
            virtual_arc_matches.append(parsed_arc)

    if virtual_arc_matches:
        arc_ids = [int(match.group("index")) for match in virtual_arc_matches]
        if len(set(arc_ids)) != len(arc_ids):
            add(messages, "error", "<v> virtualArc ids must be unique")
        if arc_ids != list(range(len(arc_ids))):
            add(
                messages,
                "warning",
                "new E-SMILES 2.0 data should order virtualArc records by "
                "consecutive ids 0..n-1",
            )
        pairs = [
            (int(match.group("from")), int(match.group("to")))
            for match in virtual_arc_matches
        ]
        canonical_pairs = [(min(start, end), max(start, end)) for start, end in pairs]
        if len(set(canonical_pairs)) != len(canonical_pairs):
            add(messages, "error", "<v> virtualArc endpoint pairs must be unique")
        if pairs != sorted(canonical_pairs):
            add(
                messages,
                "warning",
                "new E-SMILES 2.0 data should sort canonical virtualArc pairs "
                "lexicographically before assigning ids",
            )

    known_arc_ids = {int(match.group("index")) for match in virtual_arc_matches}
    for ref_id in virtual_ref_ids:
        if ref_id not in known_arc_ids:
            add(
                messages,
                "error",
                f"<r><v> references virtualArc id {ref_id}, but no matching <v> record exists",
            )

    top_level_remainder = _mask_spans(extension, spans)
    top_level_sgroups = list(SG_RE.finditer(top_level_remainder))
    if len(top_level_sgroups) > 1:
        add(messages, "error", "top-level extension may contain at most one |Sg:...| marker")
    for sg_match in top_level_sgroups:
        _validate_sgroup_count(
            sg_match.group("count"),
            messages,
            "top-level |Sg:...|",
        )

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
