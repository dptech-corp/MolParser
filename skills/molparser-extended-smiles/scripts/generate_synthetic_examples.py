#!/usr/bin/env python3
"""Validate and, after drawer freeze, render the small skill fixture set.

The default operation is read-only validation. Rendering is deliberately
fail-closed: the caller must provide the exact SHA-256 of the drawer source and
the drawer must expose every feature used by the fixture manifest.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
SKILL_DIR = SCRIPT_DIR.parent
REPO_ROOT = SKILL_DIR.parents[1]
DEFAULT_MANIFEST = SKILL_DIR / "references" / "synthetic-examples.json"
DRAWER_PATH = REPO_ROOT / "molparser" / "utils" / "drawer.py"

SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
ARC_RE = re.compile(
    r"<v>(?P<id>\d+):(?P<name>[^:]*):\[(?P<start>\d+):(?P<end>\d+)\]</v>"
)
ARC_REF_RE = re.compile(r"<r><v>(?P<id>\d+):(?P<label>.+?)</r>")
PORT_RE = re.compile(r"\[(?P<inner>\d+):(?P<outer>\d+)\]")
SGROUP_RE = re.compile(r"<g>(?P<body>.*?)</g>")
SG_COUNT_RE = re.compile(r"\|Sg:(?P<count>[A-Za-z0-9]+(?:-[A-Za-z0-9]+)?)\|")
CH2_RE = re.compile(r"<a>(?P<index>\d+):CH2\?(?P<count>[A-Za-z0-9]+)</a>")
DUMMY_RE = re.compile(r"<d>(?P<index>\d+):<dum></d>")
BALL_RE = re.compile(r"<a>(?P<index>\d+):<id>\[(?P<token>[a-z]+)\]</a>")


class FixtureError(ValueError):
    """Raised when a fixture violates the documentation contract."""


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _portable_manifest_name(path: Path) -> str:
    try:
        return path.resolve().relative_to(SKILL_DIR.resolve()).as_posix()
    except ValueError:
        return path.name


def _load_manifest(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise FixtureError(f"cannot read manifest {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise FixtureError("manifest root must be an object")
    return data


def _deep_merge(base: dict[str, Any], overlay: dict[str, Any]) -> dict[str, Any]:
    merged = copy.deepcopy(base)
    for key, value in overlay.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def _split_esmiles(esmiles: str) -> tuple[str, str]:
    if esmiles.count("<sep>") != 1:
        raise FixtureError("fixture E-SMILES must contain exactly one <sep>")
    base, extension = esmiles.split("<sep>", 1)
    if not base or not extension:
        raise FixtureError("fixture E-SMILES requires a non-empty base and extension")
    return base, extension


def _validate_rdkit_base(base: str, fixture_id: str) -> int:
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise FixtureError("RDKit is required for fixture validation") from exc
    mol = Chem.MolFromSmiles(base)
    if mol is None:
        raise FixtureError(f"{fixture_id}: base SMILES is not RDKit-parseable")
    return mol.GetNumAtoms()


def _validate_arc(example: dict[str, Any], extension: str) -> None:
    fixture_id = example["id"]
    arcs = [
        {
            "id": int(match.group("id")),
            "name": match.group("name"),
            "start": int(match.group("start")),
            "end": int(match.group("end")),
        }
        for match in ARC_RE.finditer(extension)
    ]
    if not arcs:
        raise FixtureError(f"{fixture_id}: virtual_arc fixture has no <v> record")
    ids = [arc["id"] for arc in arcs]
    if ids != list(range(len(arcs))):
        raise FixtureError(f"{fixture_id}: arc ids must be consecutive from zero")
    pairs = [(arc["start"], arc["end"]) for arc in arcs]
    if any(start >= end for start, end in pairs):
        raise FixtureError(f"{fixture_id}: new fixture arcs require START < END")
    if pairs != sorted(pairs):
        raise FixtureError(f"{fixture_id}: arc pairs must be lexicographically sorted")
    if len(set(pairs)) != len(pairs):
        raise FixtureError(f"{fixture_id}: duplicate arc pair")

    refs = [int(match.group("id")) for match in ARC_REF_RE.finditer(extension)]
    if any(ref not in ids for ref in refs):
        raise FixtureError(f"{fixture_id}: <r><v> refers to an unknown arc id")

    named = all(bool(arc["name"]) for arc in arcs)
    unnamed = all(not arc["name"] for arc in arcs)
    if not (named or unnamed):
        raise FixtureError(f"{fixture_id}: one fixture may not mix named and unnamed arcs")
    appearance = example.get("arc_appearance")
    expected = {
        "named_no_refs": (True, False),
        "named_with_refs": (True, True),
        "unnamed_with_refs": (False, True),
        "unnamed_no_refs": (False, False),
    }
    if appearance not in expected:
        raise FixtureError(f"{fixture_id}: unknown arc appearance {appearance!r}")
    if (named, bool(refs)) != expected[appearance]:
        raise FixtureError(f"{fixture_id}: arc appearance does not match its E-SMILES")


def _validate_example(example: dict[str, Any], source_sets: dict[str, Any]) -> None:
    fixture_id = example.get("id")
    if not isinstance(fixture_id, str) or not re.fullmatch(r"[a-z0-9-]+", fixture_id):
        raise FixtureError(f"invalid fixture id: {fixture_id!r}")
    esmiles = example.get("esmiles")
    if not isinstance(esmiles, str):
        raise FixtureError(f"{fixture_id}: esmiles must be a string")
    base, extension = _split_esmiles(esmiles)
    atom_count = _validate_rdkit_base(base, fixture_id)

    source = example.get("source")
    if not isinstance(source, dict) or source.get("set") not in source_sets:
        raise FixtureError(f"{fixture_id}: unknown source set")
    if not SHA256_RE.fullmatch(str(source.get("svg_sha256", ""))):
        raise FixtureError(f"{fixture_id}: invalid source SVG SHA-256")

    asset = example.get("asset")
    if not isinstance(asset, dict):
        raise FixtureError(f"{fixture_id}: missing asset object")
    asset_path = Path(str(asset.get("path", "")))
    if (
        asset_path.is_absolute()
        or ".." in asset_path.parts
        or asset_path.suffix.lower() != ".svg"
        or asset_path.parts[:2] != ("assets", "synthetic-examples")
    ):
        raise FixtureError(f"{fixture_id}: unsafe or unexpected asset path")
    asset_sha = asset.get("sha256")
    if asset_sha is not None and not SHA256_RE.fullmatch(str(asset_sha)):
        raise FixtureError(f"{fixture_id}: invalid asset SHA-256")

    compatibility = example.get("compatibility")
    if not isinstance(compatibility, dict) or set(compatibility) != {
        "esmiles_1_0",
        "esmiles_2_0",
    }:
        raise FixtureError(f"{fixture_id}: compatibility must explain both 1.0 and 2.0")

    category = example.get("category")
    if category == "whole_sru":
        dummy_indices = [int(match.group("index")) for match in DUMMY_RE.finditer(extension)]
        if len(dummy_indices) != 2 or any(index >= atom_count for index in dummy_indices):
            raise FixtureError(f"{fixture_id}: Whole SRU requires two valid <d> endpoints")
        if not SG_COUNT_RE.search(extension):
            raise FixtureError(f"{fixture_id}: Whole SRU requires |Sg:COUNT|")
    elif category == "sgroup":
        subtype = example.get("subtype")
        if subtype == "single_ch2":
            match = CH2_RE.fullmatch(extension)
            if match is None or int(match.group("index")) >= atom_count:
                raise FixtureError(f"{fixture_id}: invalid single-CH2 annotation")
        elif subtype in {"simple_short", "complex_local"}:
            match = SGROUP_RE.fullmatch(extension)
            if match is None:
                raise FixtureError(f"{fixture_id}: local repeat requires one <g> record")
            body = match.group("body")
            ports = [(int(m.group("inner")), int(m.group("outer"))) for m in PORT_RE.finditer(body)]
            if len(ports) != 2 or any(max(pair) >= atom_count for pair in ports):
                raise FixtureError(f"{fixture_id}: <g> requires two valid port pairs")
            if SG_COUNT_RE.search(body) is None:
                raise FixtureError(f"{fixture_id}: <g> requires |Sg:COUNT|")
            if example.get("bracket_style") not in {"round", "square"}:
                raise FixtureError(f"{fixture_id}: local repeat requires a bracket style")
        else:
            raise FixtureError(f"{fixture_id}: unknown sgroup subtype {subtype!r}")
    elif category == "virtual_arc":
        _validate_arc(example, extension)
    elif category == "endpoint_ball":
        match = BALL_RE.fullmatch(extension)
        if match is None or int(match.group("index")) >= atom_count:
            raise FixtureError(f"{fixture_id}: invalid endpoint-ball annotation")
        if match.group("token") != example.get("subtype"):
            raise FixtureError(f"{fixture_id}: endpoint-ball token/subtype mismatch")
    else:
        raise FixtureError(f"{fixture_id}: unknown category {category!r}")


def validate_manifest(manifest: dict[str, Any]) -> dict[str, Any]:
    if manifest.get("schema_version") != 1:
        raise FixtureError("unsupported manifest schema_version")
    source_sets = manifest.get("source_sets")
    examples = manifest.get("examples")
    if not isinstance(source_sets, dict) or not isinstance(examples, list):
        raise FixtureError("manifest requires source_sets and examples")
    if len(examples) != 10:
        raise FixtureError(f"expected exactly 10 fixtures, got {len(examples)}")

    for source_name, source in source_sets.items():
        if not isinstance(source, dict):
            raise FixtureError(f"source set {source_name!r} must be an object")
        for key, value in source.items():
            if key.endswith("sha256") and not SHA256_RE.fullmatch(str(value)):
                raise FixtureError(f"source set {source_name!r} has invalid {key}")

    for example in examples:
        if not isinstance(example, dict):
            raise FixtureError("every fixture must be an object")
        _validate_example(example, source_sets)

    ids = [example["id"] for example in examples]
    esmiles = [example["esmiles"] for example in examples]
    assets = [example["asset"]["path"] for example in examples]
    for label, values in (("id", ids), ("E-SMILES", esmiles), ("asset path", assets)):
        if len(set(values)) != len(values):
            raise FixtureError(f"duplicate fixture {label}")

    category_counts = Counter(example["category"] for example in examples)
    expected_categories = Counter(
        {"whole_sru": 1, "sgroup": 3, "virtual_arc": 4, "endpoint_ball": 2}
    )
    if category_counts != expected_categories:
        raise FixtureError(f"category coverage mismatch: {dict(category_counts)}")

    sgroup_subtypes = Counter(
        example["subtype"] for example in examples if example["category"] == "sgroup"
    )
    if sgroup_subtypes != Counter(
        {"single_ch2": 1, "simple_short": 1, "complex_local": 1}
    ):
        raise FixtureError(f"sgroup subtype coverage mismatch: {dict(sgroup_subtypes)}")

    arc_appearances = {
        example["arc_appearance"]
        for example in examples
        if example["category"] == "virtual_arc"
    }
    if arc_appearances != {
        "named_no_refs",
        "named_with_refs",
        "unnamed_with_refs",
        "unnamed_no_refs",
    }:
        raise FixtureError(f"arc appearance coverage mismatch: {sorted(arc_appearances)}")

    pending_assets = sum(example["asset"]["sha256"] is None for example in examples)
    return {
        "all_valid": True,
        "fixture_count": len(examples),
        "category_counts": dict(sorted(category_counts.items())),
        "sgroup_subtype_counts": dict(sorted(sgroup_subtypes.items())),
        "arc_appearances": sorted(arc_appearances),
        "pending_asset_count": pending_assets,
    }


def verify_asset_state(
    manifest: dict[str, Any], output_root: Path = SKILL_DIR
) -> dict[str, Any]:
    status = manifest.get("status")
    if status not in {"pending-renderer-freeze", "frozen"}:
        raise FixtureError(f"unknown manifest status: {status!r}")
    verified: list[dict[str, Any]] = []
    for example in manifest["examples"]:
        target = output_root / example["asset"]["path"]
        expected_sha = example["asset"]["sha256"]
        if status == "pending-renderer-freeze":
            if expected_sha is not None or target.exists():
                raise FixtureError(
                    f"{example['id']}: pending manifest must not contain a rendered asset"
                )
            continue
        if expected_sha is None:
            raise FixtureError(f"{example['id']}: frozen asset has no SHA-256")
        if not target.is_file():
            raise FixtureError(f"{example['id']}: frozen asset is missing")
        actual_sha = _sha256_file(target)
        if actual_sha != expected_sha:
            raise FixtureError(
                f"{example['id']}: asset SHA mismatch, expected {expected_sha}, got {actual_sha}"
            )
        svg = target.read_text(encoding="utf-8")
        metrics = _validate_svg(svg, example["id"])
        verified.append(
            {
                "id": example["id"],
                "asset": example["asset"]["path"],
                "sha256": actual_sha,
                "size": target.stat().st_size,
                **metrics,
            }
        )
    return {
        "asset_state": status,
        "verified_asset_count": len(verified),
        "verified_assets": verified,
    }


def _model_fields(model_class: Any) -> set[str]:
    fields = getattr(model_class, "model_fields", None)
    if fields is None:
        fields = getattr(model_class, "__fields__", {})
    return set(fields)


def _load_frozen_drawer(expected_sha256: str) -> Any:
    if not SHA256_RE.fullmatch(expected_sha256):
        raise FixtureError("--renderer-sha256 must be a lowercase SHA-256")
    actual_sha256 = _sha256_file(DRAWER_PATH)
    if actual_sha256 != expected_sha256:
        raise FixtureError(
            f"drawer SHA mismatch: expected {expected_sha256}, got {actual_sha256}"
        )
    sys.path.insert(0, str(REPO_ROOT))
    from molparser.utils import drawer

    required = {
        drawer.StylingConfig: {"sgroup_bracket_style"},
        drawer.FeaturesConfig: {
            "repeat_units",
            "academic_virtual_arcs",
            "endpoint_balls",
        },
    }
    missing: list[str] = []
    for model, field_names in required.items():
        missing.extend(
            f"{model.__name__}.{name}"
            for name in sorted(field_names - _model_fields(model))
        )
    if missing:
        raise FixtureError(
            "drawer is not feature-complete for these fixtures: " + ", ".join(missing)
        )
    return drawer


def _validate_svg(svg: str, fixture_id: str) -> dict[str, Any]:
    try:
        root = ET.fromstring(svg)
    except ET.ParseError as exc:
        raise FixtureError(f"{fixture_id}: rendered SVG is not XML: {exc}") from exc
    if root.tag.rsplit("}", 1)[-1] != "svg":
        raise FixtureError(f"{fixture_id}: rendered document root is not SVG")
    view_box = root.get("viewBox", "").split()
    if len(view_box) != 4:
        raise FixtureError(f"{fixture_id}: SVG has no four-number viewBox")
    try:
        width, height = float(view_box[2]), float(view_box[3])
    except ValueError as exc:
        raise FixtureError(f"{fixture_id}: invalid SVG viewBox") from exc
    if width <= 0 or height <= 0:
        raise FixtureError(f"{fixture_id}: non-positive SVG viewBox")
    for node in root.iter():
        if node.tag.rsplit("}", 1)[-1].lower() == "script":
            raise FixtureError(f"{fixture_id}: SVG contains a script")
        for attr_name, attr_value in node.attrib.items():
            if attr_name.rsplit("}", 1)[-1] == "href" and re.match(
                r"(?i)^(?:https?|file|javascript):", attr_value.strip()
            ):
                raise FixtureError(f"{fixture_id}: SVG contains an external reference")
    return {"viewBox": [float(value) for value in view_box]}


def render_examples(
    manifest: dict[str, Any],
    manifest_path: Path,
    output_root: Path,
    renderer_sha256: str,
    overwrite: bool,
) -> dict[str, Any]:
    drawer = _load_frozen_drawer(renderer_sha256)
    defaults = manifest["render_defaults"]["config"]
    rendered: list[dict[str, Any]] = []
    payloads: list[tuple[Path, bytes]] = []

    for example in manifest["examples"]:
        config = _deep_merge(defaults, example.get("render_config", {}))
        svg = drawer.draw(example["esmiles"], config=config, output_format="svg")
        if not isinstance(svg, str):
            raise FixtureError(f"{example['id']}: drawer did not return SVG text")
        metrics = _validate_svg(svg, example["id"])
        payload = svg.encode("utf-8")
        target = output_root / example["asset"]["path"]
        if target.exists() and not overwrite:
            raise FixtureError(f"refusing to overwrite existing asset: {target}")
        payloads.append((target, payload))
        rendered.append(
            {
                "id": example["id"],
                "asset": example["asset"]["path"],
                "sha256": _sha256_bytes(payload),
                "size": len(payload),
                **metrics,
            }
        )

    for target, payload in payloads:
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = target.with_suffix(target.suffix + ".tmp")
        temporary.write_bytes(payload)
        temporary.replace(target)

    return {
        "all_valid": True,
        "manifest": _portable_manifest_name(manifest_path),
        "manifest_sha256": _sha256_file(manifest_path),
        "drawer_sha256": renderer_sha256,
        "asset_count": len(rendered),
        "assets": rendered,
    }


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate or render the ten audited E-SMILES skill fixtures."
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Validate only (the default when --render is absent).",
    )
    parser.add_argument("--render", action="store_true", help="Render SVG assets.")
    parser.add_argument(
        "--renderer-sha256",
        help="Required with --render; must match molparser/utils/drawer.py.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=SKILL_DIR,
        help="Root under which manifest-relative asset paths are written.",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--report",
        type=Path,
        help="Optional JSON report path. Without it, the report is printed only.",
    )
    args = parser.parse_args(argv)
    if args.render and args.check_only:
        parser.error("choose --check-only or --render, not both")
    if args.render and not args.renderer_sha256:
        parser.error("--render requires --renderer-sha256")
    return args


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(sys.argv[1:] if argv is None else argv)
    try:
        manifest_path = args.manifest.resolve()
        manifest = _load_manifest(manifest_path)
        validation = validate_manifest(manifest)
        report: dict[str, Any] = {
            "mode": "render" if args.render else "check-only",
            "manifest": _portable_manifest_name(manifest_path),
            "manifest_sha256": _sha256_file(manifest_path),
            **validation,
        }
        if args.render:
            report.update(
                render_examples(
                    manifest,
                    manifest_path,
                    args.output_root.resolve(),
                    args.renderer_sha256,
                    args.overwrite,
                )
            )
        else:
            report.update(verify_asset_state(manifest))
    except FixtureError as exc:
        print(json.dumps({"all_valid": False, "error": str(exc)}, ensure_ascii=False))
        return 1

    rendered_report = json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.report.with_suffix(args.report.suffix + ".tmp")
        temporary.write_text(rendered_report, encoding="utf-8")
        temporary.replace(args.report)
    print(rendered_report, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
