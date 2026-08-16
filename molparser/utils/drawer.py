"""SVG renderer for E-SMILES captions."""

from __future__ import annotations

import logging
import math
import random
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple, Union

from pydantic import BaseModel, Field
from rdkit import Chem, RDLogger
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdchem import Mol
from rdkit.Geometry import Point2D

try:
    from .translator import (
        AtomIndex,
        GroupDesc,
        Patterns,
        RingIndex,
        TextType,
        Tokens,
    )
except ImportError:  # Support running from the package directory as working directory.
    from translator import (
        AtomIndex,
        GroupDesc,
        Patterns,
        RingIndex,
        TextType,
        Tokens,
    )


__all__ = ["DrawingConfig", "draw"]


logger = logging.getLogger(__name__)

_SVG_NAMESPACE = "http://www.w3.org/2000/svg"

# `<id>[...]` is otherwise an open-ended textual annotation.  Keep endpoint
# balls behind an explicit allow-list so existing patent labels continue to be
# rendered as text.  Colours use RGB values in RDKit's 0..1 range.
_ENDPOINT_BALL_STYLES: Dict[
    str,
    Tuple[Tuple[float, float, float], Tuple[float, float, float]],
] = {
    "ball": ((1.00, 1.00, 1.00), (0.08, 0.08, 0.08)),
    "grey": ((0.78, 0.78, 0.78), (0.28, 0.28, 0.28)),
    "black": ((0.00, 0.00, 0.00), (0.00, 0.00, 0.00)),
    "green": ((0.62, 0.89, 0.76), (0.07, 0.53, 0.31)),
    "blue": ((0.66, 0.73, 0.92), (0.18, 0.30, 0.60)),
    "yellow": ((0.97, 0.85, 0.42), (0.66, 0.47, 0.00)),
    "purple": ((0.78, 0.65, 0.91), (0.44, 0.25, 0.63)),
    "orange": ((0.97, 0.75, 0.51), (0.78, 0.42, 0.09)),
    "pink": ((0.90, 0.65, 0.76), (0.61, 0.19, 0.36)),
    "brown": ((0.73, 0.55, 0.42), (0.42, 0.25, 0.15)),
}
_ENDPOINT_BALL_RADIUS = 0.43


class VisualConfig(BaseModel):
    # Do not change these defaults: callers and tracked figures rely on the
    # exact SVG emitted by the original drawer.
    padding: float = Field(default=0.1, ge=0.0, le=1.0)
    additionalAtomLabelPadding: float = Field(default=0.05, ge=0.0, le=0.5)
    fixedFontSize: int = Field(default=10, ge=1, le=100)
    # RDKit accepts an explicit TTF/OTF path through MolDrawOptions.fontFile.
    # Keep this optional so existing callers retain the built-in font exactly.
    fontFile: Optional[str] = Field(default=None)
    # Optional fail-closed escape hatch for unusually wide font glyphs in
    # virtual-arc annotations.  None preserves the historical arc layout
    # byte-for-byte; controlled augmentation may request modest extra spacing.
    virtualArcLabelSpacingScale: Optional[float] = Field(
        default=None, ge=1.0, le=2.0
    )
    bondLineWidth: int = Field(default=1, ge=1, le=10)
    multipleBondOffset: float = Field(default=0.1, ge=0.01, le=1.0)
    dummiesAreAttachments: bool = Field(default=False)
    addAtomIndices: Union[bool, float] = Field(default=False)
    singleColourWedgeBonds: bool = Field(default=True)
    legendFontSize: int = Field(default=9, ge=1, le=50)


class StylingConfig(BaseModel):
    palette: Union[Literal["cdk", "bw"], float] = Field(default="cdk")
    use_modern_symbols: bool = Field(default=True)
    enhanced_contrast: bool = Field(default=True)
    highlight_atoms: Union[bool, float] = Field(default=False)
    highlight_bonds: Union[bool, float] = Field(default=False)
    custom_colors: Union[bool, float] = Field(default=False)
    bold_r_groups: bool = Field(default=True)
    ring_connector_style: Literal["dashed", "solid"] = Field(default="solid")
    ring_connector_extension: float = Field(default=0.8, ge=0.0, le=5.0)
    dummy_line_style: Literal["wavy", "dashed", "dotted", "solid"] = Field(default="wavy")
    show_default_sru_count: bool = Field(default=True)
    # Local E-SMILES 2.0 <g> repeats may use either publication-style
    # parentheses or square brackets.  The caller resolves the style for a
    # record; the drawer itself remains deterministic.  Whole-SRU rendering
    # deliberately ignores this option and stays round for compatibility.
    sgroup_bracket_style: Literal["round", "square"] = Field(default="round")


class FeaturesConfig(BaseModel):
    dummy_atoms: bool = Field(default=True)
    circled_r_groups: bool = Field(default=True)
    ring_annotations: bool = Field(default=True)
    # None selects syntax-aware compatibility mode: new constructs render
    # academically while already-supported 1.0 captions keep the legacy path.
    # Explicit False disables an overlay and True forces its new implementation.
    repeat_units: Optional[bool] = Field(default=None)
    academic_virtual_arcs: Optional[bool] = Field(default=None)
    endpoint_balls: Optional[bool] = Field(default=None)


class DrawingConfig(BaseModel):
    visual: VisualConfig = Field(default_factory=VisualConfig)
    styling: StylingConfig = Field(default_factory=StylingConfig)
    features: FeaturesConfig = Field(default_factory=FeaturesConfig)

    @classmethod
    def from_dict(cls, config_dict: Optional[Dict[str, Any]] = None) -> "DrawingConfig":
        if config_dict is None:
            return cls()
        merged_dict: Dict[str, Dict] = {"visual": {}, "styling": {}, "features": {}}
        for section in ("visual", "styling", "features"):
            if section in config_dict and isinstance(config_dict[section], dict):
                merged_dict[section] = config_dict[section]
        return cls(**merged_dict)


def _validate_config(config: Optional[Dict[str, Any]] = None) -> DrawingConfig:
    return DrawingConfig.from_dict(config)


def _resolve_probability(value: Union[bool, float]) -> bool:
    if isinstance(value, bool):
        return value
    prob = min(max(value, 0.0), 1.0)
    return random.random() < prob


def _random_color() -> Tuple[float, float, float]:
    return (
        random.random() / 2 + 0.5,
        random.random() / 2 + 0.5,
        random.random() / 2 + 0.5,
    )


def _update_bounds(bounds: Optional[List[float]], x: float, y: float) -> List[float]:
    if bounds is None:
        return [x, y, x, y]
    bounds[0] = min(bounds[0], x)
    bounds[1] = min(bounds[1], y)
    bounds[2] = max(bounds[2], x)
    bounds[3] = max(bounds[3], y)
    return bounds


def _estimate_text_bounds(
    label_text: str,
    x: float,
    y: float,
    font_px: float,
    margin: float,
    strip_markup: bool = True,
    left_aligned: bool = False,
) -> Tuple[float, float, float, float]:
    visible = re.sub(r"<[^>]+>", "", label_text) if strip_markup else label_text
    text_len = max(len(visible), 1)
    if left_aligned:
        # RDKit's DrawString(..., orient=1) places fallback annotation text to
        # the right of its anchor.  A one-sided estimate avoids doubling the
        # canvas width while retaining a conservative glyph-width allowance.
        text_w = font_px * 0.45 * text_len + 2 * margin
    else:
        text_w = font_px * 0.6 * text_len + 2 * margin
    text_h = font_px + 2 * margin
    if left_aligned:
        return x - margin, y - text_h, x + text_w, y + text_h
    return x - text_w, y - text_h, x + text_w, y + text_h


def _format_symbol_with_subscripts(symbol: str) -> str:
    if "<" in symbol and ">" in symbol:
        return symbol
    return re.sub(r"(\d+)", r"<sub>\1</sub>", symbol)


def _ring_annotation_screen_slots(count: int) -> List[Tuple[float, float]]:
    """Return stable screen-space slots for regio-uncertain ring labels."""

    if count <= 0:
        return []
    diagonal = math.sqrt(0.5)
    ordered = [
        (-1.0, 0.0),
        (1.0, 0.0),
        (0.0, -1.0),
        (0.0, 1.0),
        (-diagonal, -diagonal),
        (diagonal, -diagonal),
        (-diagonal, diagonal),
        (diagonal, diagonal),
    ]
    if count == 1:
        return [(0.0, -1.0)]
    if count <= len(ordered):
        return ordered[:count]
    return [
        (math.cos(2.0 * math.pi * idx / count), math.sin(2.0 * math.pi * idx / count))
        for idx in range(count)
    ]


def _format_group_label(desc: GroupDesc) -> str:
    if desc.symbol == Tokens.special_id and desc.script:
        return desc.script
    label = desc.symbol or ""
    if desc.script is not None:
        label += f"<sub>{desc.script}</sub>"
    if desc.prime is not None:
        label += desc.prime
    if desc.multiple is not None:
        label = "(" + label + ")" + f"<sub>{desc.multiple}</sub>"
    return label


def _endpoint_ball_token(desc: GroupDesc) -> Optional[str]:
    """Return the approved ball token encoded by an atom-group description."""
    if (
        desc.symbol != Tokens.special_id
        or desc.script not in _ENDPOINT_BALL_STYLES
        or desc.prime is not None
        or desc.multiple is not None
        or desc.is_circle
        or desc.is_dummy
    ):
        return None
    return desc.script


def _draw_endpoint_ball(
    drawer,
    center: Point2D,
    token: str,
    line_width: float,
) -> None:
    """Draw one deterministic flat ball from the explicit allow-list."""

    fill_color, outline_color = _ENDPOINT_BALL_STYLES[token]
    radius = _ENDPOINT_BALL_RADIUS
    lower_left = Point2D(center.x - radius, center.y - radius)
    upper_right = Point2D(center.x + radius, center.y + radius)

    # Painting after the molecule covers the bond underneath the disc, so the
    # visible bond naturally terminates at the ball boundary.
    drawer.SetFillPolys(True)
    drawer.SetColour(fill_color)
    drawer.DrawEllipse(lower_left, upper_right)
    drawer.SetFillPolys(False)
    drawer.SetColour(outline_color)
    drawer.SetLineWidth(line_width)
    drawer.DrawEllipse(lower_left, upper_right)


def _truncate_annotation(text: str, limit: int = 64) -> str:
    compact = re.sub(r"\s+", " ", text).strip()
    return compact if len(compact) <= limit else compact[: limit - 3] + "..."


def _precompat_annotations(
    groups: str,
    rendered_sgroups: Optional[set[str]] = None,
) -> List[str]:
    if rendered_sgroups is None:
        # Preserve the original fallback text path exactly when academic
        # repeat overlays are not requested.
        annotations: List[str] = []
        for match in re.finditer(r"<s>(.*?)</s>", groups, re.DOTALL):
            annotations.append(
                f"substructure: {_truncate_annotation(match.group(1))}"
            )
        for match in re.finditer(r"<g>(.*?)</g>", groups, re.DOTALL):
            annotations.append(f"s-group: {_truncate_annotation(match.group(1))}")
        return annotations

    annotations: List[str] = []
    for match in re.finditer(r"<s>(.*?)</s>", groups, re.DOTALL):
        annotations.append(f"substructure: {_truncate_annotation(match.group(1))}")
    rendered = rendered_sgroups or set()
    outer_groups = re.sub(r"<s>.*?</s>", "", groups, flags=re.DOTALL)
    for match in re.finditer(r"<g>(.*?)</g>", outer_groups, re.DOTALL):
        body = match.group(1)
        if body not in rendered:
            annotations.append(f"s-group: {_truncate_annotation(body)}")
    return annotations


def _sgroup_records(groups: str) -> List[Tuple[str, List[Tuple[int, int]], str]]:
    """Parse pre-compatible S-group repeat records for depiction.

    The first index in each pair is the atom inside the repeated fragment and
    the second is the atom outside it.  The parser is intentionally strict:
    malformed records remain available to the fallback annotation renderer.
    """

    records: List[Tuple[str, List[Tuple[int, int]], str]] = []
    outer_groups = re.sub(r"<s>.*?</s>", "", groups, flags=re.DOTALL)
    for match in re.finditer(r"<g>(.*?)</g>", outer_groups, re.DOTALL):
        body = match.group(1)
        repeat = re.search(r"\|Sg:([^|]+)\|", body)
        if repeat is None:
            continue
        prefix = body[: repeat.start()]
        suffix = body[repeat.end() :]
        ports = [
            (int(inner), int(outer))
            for inner, outer in re.findall(r"\[(\d+):(\d+)\]", prefix)
        ]
        leftover = re.sub(r"\[\d+:\d+\]", "", prefix).replace(":", "").strip()
        if (
            not ports
            or len(ports) != len(set(ports))
            or leftover
            or suffix.strip()
        ):
            continue
        count = repeat.group(1).strip()
        if not re.fullmatch(r"[A-Za-z0-9]+(?:-[A-Za-z0-9]+)?", count):
            continue
        records.append((body, ports, count))
    return records


def _single_ch2_repeat_records(
    mol: Mol,
    descriptions: List[GroupDesc],
) -> List[Tuple[int, List[Tuple[int, int]], str]]:
    """Return valid inline ``<a>ATOM:CH2?COUNT</a>`` repeat records.

    The source atom remains a degree-two dummy in canonical E-SMILES.  For
    depiction it acts as an unlabeled carbon-chain vertex, while the two
    incident single bonds define the repeat boundaries.  Keeping this
    recognition deliberately narrow prevents ordinary atom abbreviations from
    acquiring polymer parentheses.
    """

    records: List[Tuple[int, List[Tuple[int, int]], str]] = []
    seen: set[int] = set()
    atom_description_counts: Dict[int, int] = {}
    for desc in descriptions:
        if isinstance(desc.id, AtomIndex):
            atom_idx = int(desc.id)
            atom_description_counts[atom_idx] = (
                atom_description_counts.get(atom_idx, 0) + 1
            )
    for desc in descriptions:
        if (
            not isinstance(desc.id, AtomIndex)
            or desc.is_dummy
            or desc.multiple is None
            or desc.script is not None
            or desc.prime is not None
        ):
            continue
        atom_idx = int(desc.id)
        if atom_description_counts.get(atom_idx) != 1:
            continue
        visible_symbol = re.sub(r"<[^>]+>", "", desc.symbol or "")
        if visible_symbol != "CH2" or atom_idx in seen:
            continue
        if not 0 <= atom_idx < mol.GetNumAtoms():
            continue
        atom = mol.GetAtomWithIdx(atom_idx)
        bonds = list(atom.GetBonds())
        if (
            atom.GetSymbol() != "*"
            or atom.GetDegree() != 2
            or len(bonds) != 2
            or any(
                bond.GetBondType() != Chem.BondType.SINGLE
                or bond.GetIsAromatic()
                or bond.IsInRing()
                for bond in bonds
            )
        ):
            continue
        ports = sorted(
            (atom_idx, bond.GetOtherAtomIdx(atom_idx)) for bond in bonds
        )
        if len({outer for _inner, outer in ports}) != 2:
            continue
        seen.add(atom_idx)
        records.append((atom_idx, ports, desc.multiple))
    return records


def _sgroup_repeat_atoms(
    mol: Mol,
    ports: List[Tuple[int, int]],
) -> Optional[List[int]]:
    """Infer the repeated atom component after removing its boundary bonds."""

    if not ports:
        return None
    atom_count = mol.GetNumAtoms()
    if len(ports) != len(set(ports)):
        return None
    removed: set[Tuple[int, int]] = set()
    inner_atoms: set[int] = set()
    outer_atoms: set[int] = set()
    for inner, outer in ports:
        if (
            inner == outer
            or not 0 <= inner < atom_count
            or not 0 <= outer < atom_count
            or mol.GetBondBetweenAtoms(inner, outer) is None
        ):
            return None
        inner_atoms.add(inner)
        outer_atoms.add(outer)
        removed.add(tuple(sorted((inner, outer))))

    component: set[int] = set()
    pending = [next(iter(inner_atoms))]
    while pending:
        atom_idx = pending.pop()
        if atom_idx in component:
            continue
        component.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        for bond in atom.GetBonds():
            other_idx = bond.GetOtherAtomIdx(atom_idx)
            if tuple(sorted((atom_idx, other_idx))) in removed:
                continue
            if other_idx not in component:
                pending.append(other_idx)

    if not inner_atoms.issubset(component) or component & outer_atoms:
        return None
    if len(component) == atom_count:
        return None
    return sorted(component)


def _whole_sru_endpoints(mol: Mol, groups: str) -> Optional[List[int]]:
    """Return the two validated terminal dummy atoms of a whole SRU.

    E-SMILES 2.0 identifies the endpoints explicitly with two ``<d>``
    records.  The original MolParser figure set also contains a legacy SRU
    example with only two terminal ``*`` atoms and a top-level ``|Sg:n|``.
    Keep that unambiguous form drawable, but never guess when the molecule
    contains any other dummy atom or a partial/malformed ``<d>`` annotation.
    """

    top_level_groups = re.sub(
        r"<s>.*?</s>|<g>.*?</g>|<v>.*?</v>",
        "",
        groups,
        flags=re.DOTALL,
    )
    matches = re.findall(r"<d>(\d+):<dum></d>", top_level_groups)
    if matches:
        if len(matches) != 2 or len(set(matches)) != 2:
            return None
        dummy_indices = {int(value) for value in matches}
    else:
        # Legacy, pre-2.0 fallback.  Requiring exactly two dummy atoms in the
        # entire graph avoids interpreting ordinary Markush placeholders as
        # polymer endpoints.
        dummy_indices = {
            atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0
        }
        if len(dummy_indices) != 2:
            return None
    for atom_idx in dummy_indices:
        if not 0 <= atom_idx < mol.GetNumAtoms():
            return None
        atom = mol.GetAtomWithIdx(atom_idx)
        bonds = list(atom.GetBonds())
        if (
            atom.GetSymbol() != "*"
            or atom.GetDegree() != 1
            or len(bonds) != 1
            or bonds[0].GetBondType() != Chem.BondType.SINGLE
            or bonds[0].IsInRing()
        ):
            return None
    return sorted(dummy_indices)


def _whole_sru_atoms(mol: Mol, groups: str) -> Optional[List[int]]:
    """Return the connected core atoms for a well-formed two-ended whole SRU."""

    endpoint_indices = _whole_sru_endpoints(mol, groups)
    if endpoint_indices is None:
        return None
    dummy_indices = set(endpoint_indices)

    core = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetIdx() not in dummy_indices
    ]
    if not core:
        return None
    allowed = set(core)
    seen: set[int] = set()
    pending = [core[0]]
    while pending:
        atom_idx = pending.pop()
        if atom_idx in seen:
            continue
        seen.add(atom_idx)
        for bond in mol.GetAtomWithIdx(atom_idx).GetBonds():
            other_idx = bond.GetOtherAtomIdx(atom_idx)
            if other_idx in allowed and other_idx not in seen:
                pending.append(other_idx)
    return core if seen == allowed else None


def _whole_sru_ports(mol: Mol, groups: str) -> Optional[List[Tuple[int, int]]]:
    """Return inner-core/outer-dummy boundary pairs for a valid whole SRU."""

    endpoints = _whole_sru_endpoints(mol, groups)
    if endpoints is None or _whole_sru_atoms(mol, groups) is None:
        return None
    ports: List[Tuple[int, int]] = []
    for dummy_idx in endpoints:
        atom = mol.GetAtomWithIdx(dummy_idx)
        neighbor = next(iter(atom.GetNeighbors()))
        ports.append((neighbor.GetIdx(), dummy_idx))
    return ports


def _sru_count(extension: str) -> Optional[str]:
    match = re.fullmatch(
        r"\|Sg:([A-Za-z0-9]+(?:-[A-Za-z0-9]+)?)\|",
        extension,
    )
    if match is None:
        return None
    return match.group(1)


def _draw_repeat_count(
    drawer,
    count: str,
    label_point: Point2D,
    font_size: float,
    bounds: Optional[List[float]],
) -> Optional[List[float]]:
    """Draw a compact lower-right repeat subscript without leaking font state."""

    # Academic polymer subscripts are visibly smaller than atom labels.  A
    # slightly smaller scale is especially important for multi-digit counts.
    scale = 0.72 if len(count) > 1 else 0.75
    subscript_font = max(8.0, font_size * scale)
    previous_font = drawer.FontSize() if hasattr(drawer, "FontSize") else None
    if hasattr(drawer, "SetFontSize"):
        drawer.SetFontSize(subscript_font)
    try:
        drawer.DrawString(count, label_point, 1, True)
    finally:
        if previous_font is not None and hasattr(drawer, "SetFontSize"):
            drawer.SetFontSize(previous_font)
    label_bounds = _estimate_text_bounds(
        count,
        label_point.x,
        label_point.y,
        font_px=subscript_font,
        margin=4.0,
    )
    bounds = _update_bounds(bounds, label_bounds[0], label_bounds[1])
    bounds = _update_bounds(bounds, label_bounds[2], label_bounds[3])
    return bounds


def _parenthesis_curve(
    midpoint: Point2D,
    axis_x: float,
    axis_y: float,
    opening_x: float,
    opening_y: float,
    half_height: float,
) -> Tuple[Point2D, Point2D, Point2D, Point2D]:
    """Build the shared slim parenthesis used by whole and nested repeats."""

    endpoint_inset = min(max(half_height * 0.06, 0.03), 0.07)
    outward_bulge = min(max(half_height * 0.18, 0.07), 0.16)
    control_height = half_height * 0.43
    return (
        Point2D(
            midpoint.x + axis_x * half_height + opening_x * endpoint_inset,
            midpoint.y + axis_y * half_height + opening_y * endpoint_inset,
        ),
        Point2D(
            midpoint.x + axis_x * control_height - opening_x * outward_bulge,
            midpoint.y + axis_y * control_height - opening_y * outward_bulge,
        ),
        Point2D(
            midpoint.x - axis_x * control_height - opening_x * outward_bulge,
            midpoint.y - axis_y * control_height - opening_y * outward_bulge,
        ),
        Point2D(
            midpoint.x - axis_x * half_height + opening_x * endpoint_inset,
            midpoint.y - axis_y * half_height + opening_y * endpoint_inset,
        ),
    )


def _snap_paired_bracket_axes(
    axes: List[Tuple[float, float]],
    threshold_degrees: float = 45.0,
) -> List[Tuple[float, float]]:
    """Snap a more-vertical-than-horizontal pair to parallel vertical axes.

    Snapping is deliberately all-or-nothing: a single oblique boundary keeps
    both parentheses perpendicular to their own bonds.  This avoids mixing a
    typographic vertical bracket with a chemically aligned oblique partner.
    The threshold is inclusive; at the default 45 degrees, ``abs(axis_y) >=
    abs(axis_x)`` is the exact eligibility rule for each side.
    """

    if len(axes) != 2:
        return list(axes)
    if math.isclose(threshold_degrees, 45.0, abs_tol=1e-12):
        eligible = all(
            abs(axis_y) + 1e-12 >= abs(axis_x)
            for axis_x, axis_y in axes
        )
    else:
        eligible = all(
            math.degrees(math.atan2(abs(axis_x), abs(axis_y)))
            <= threshold_degrees + 1e-12
            for axis_x, axis_y in axes
        )
    if not eligible:
        return list(axes)
    return [
        (0.0, 1.0 if axis_y >= 0.0 else -1.0)
        for _axis_x, axis_y in axes
    ]


def _draw_sgroup_boundary_parentheses(
    drawer,
    mol: Mol,
    ports: List[Tuple[int, int]],
    count: str,
    line_width: float,
    font_size: float,
    bounds: Optional[List[float]],
    show_count: bool = True,
    bracket_style: Literal["round", "square"] = "round",
) -> Optional[List[float]]:
    """Draw an inward-opening bracket across each S-group boundary bond."""

    conf = mol.GetConformer()
    draw_points: List[Point2D] = []
    draw_curves: List[List[Point2D]] = []
    drawer.SetColour((0.05, 0.05, 0.05))
    drawer.SetLineWidth(line_width)
    layouts = []
    for inner_idx, outer_idx in ports:
        inner = conf.GetAtomPosition(inner_idx)
        outer = conf.GetAtomPosition(outer_idx)
        inner_dir_x = inner.x - outer.x
        inner_dir_y = inner.y - outer.y
        length = math.hypot(inner_dir_x, inner_dir_y)
        if length <= 1e-6:
            continue
        inner_dir_x /= length
        inner_dir_y /= length
        perp_x = -inner_dir_y
        perp_y = inner_dir_x
        midpoint = Point2D((inner.x + outer.x) / 2.0, (inner.y + outer.y) / 2.0)
        layouts.append(
            [midpoint, inner_dir_x, inner_dir_y, perp_x, perp_y]
        )
    snapped_axes = _snap_paired_bracket_axes(
        [(layout[3], layout[4]) for layout in layouts]
    )
    for layout, (perp_x, perp_y) in zip(layouts, snapped_axes):
        midpoint, inner_dir_x, inner_dir_y, _old_perp_x, _old_perp_y = layout
        if bracket_style == "square":
            half_height = 0.52
            cap_length = 0.19
            upper = Point2D(
                midpoint.x + perp_x * half_height,
                midpoint.y + perp_y * half_height,
            )
            lower = Point2D(
                midpoint.x - perp_x * half_height,
                midpoint.y - perp_y * half_height,
            )
            upper_inner = Point2D(
                upper.x + inner_dir_x * cap_length,
                upper.y + inner_dir_y * cap_length,
            )
            lower_inner = Point2D(
                lower.x + inner_dir_x * cap_length,
                lower.y + inner_dir_y * cap_length,
            )
            drawer.DrawLine(upper, lower)
            drawer.DrawLine(upper, upper_inner)
            drawer.DrawLine(lower, lower_inner)
            curve_points = [upper, lower, upper_inner, lower_inner]
        else:
            upper, control_upper, control_lower, lower = _parenthesis_curve(
                midpoint,
                perp_x,
                perp_y,
                inner_dir_x,
                inner_dir_y,
                0.52,
            )
            curve_points = _draw_cubic_curve(
                drawer,
                upper,
                control_upper,
                control_lower,
                lower,
                segments=16,
            )
        curve_draw_points: List[Point2D] = []
        for point in curve_points:
            draw_point = drawer.GetDrawCoords(point)
            draw_points.append(draw_point)
            curve_draw_points.append(draw_point)
            bounds = _update_bounds(bounds, draw_point.x, draw_point.y)
        draw_curves.append(curve_draw_points)

    if not draw_points:
        return bounds
    if not show_count:
        return bounds
    right_curve = max(
        draw_curves,
        key=lambda points: max(point.x for point in points),
    )
    label_point = Point2D(
        max(point.x for point in right_curve) + 3.0,
        max(point.y for point in right_curve) + font_size * 0.32,
    )
    bounds = _draw_repeat_count(
        drawer,
        count,
        label_point,
        font_size,
        bounds,
    )
    return bounds


def _virtual_arcs(groups: str) -> List[Tuple[int, str, int, int]]:
    arcs: List[Tuple[int, str, int, int]] = []
    seen_indices: set[int] = set()
    groups = re.sub(r"<s>.*?</s>", "", groups, flags=re.DOTALL)
    pattern = r"<v>(?P<idx>\d+):(?P<name>[^:]*):\[(?P<start>\d+):(?P<end>\d+)\]</v>"
    for match in re.finditer(pattern, groups, re.DOTALL):
        arc_idx = int(match.group("idx"))
        if arc_idx in seen_indices or len(arcs) >= 64:
            continue
        seen_indices.add(arc_idx)
        arcs.append(
            (
                arc_idx,
                match.group("name"),
                int(match.group("start")),
                int(match.group("end")),
            )
        )
    return arcs


def _virtual_arc_substituents(groups: str) -> Dict[int, List[str]]:
    substituents: Dict[int, List[str]] = {}
    groups = re.sub(r"<s>.*?</s>", "", groups, flags=re.DOTALL)
    pattern = r"<r><v>(?P<idx>\d+):(?P<label>.+?)</r>"
    for match in re.finditer(pattern, groups, re.DOTALL):
        idx = int(match.group("idx"))
        label = match.group("label")
        label_match = re.match(
            r"(?P<symbol>[^\[\?\'\"]+)"
            r"(?P<script>\[[^\]]+\])?"
            r"(?P<prime>[\'\"]?)"
            r"(?P<multiple>\?(?:[a-z]|\d+|\d+-\d+)?)?$",
            label,
        )
        if label_match:
            symbol = label_match.group("symbol") or ""
            script = label_match.group("script")
            prime = label_match.group("prime") or ""
            multiple = label_match.group("multiple") or ""
            if script:
                symbol += f"<sub>{script[1:-1]}</sub>"
            label = symbol + prime
            if multiple:
                label = f"({label})<sub>{multiple[1:]}</sub>"
        labels = substituents.setdefault(idx, [])
        if len(labels) < 16:
            labels.append(label)
    return substituents


def _path_centroid(
    mol: Mol,
    conf,
    start_idx: int,
    end_idx: int,
) -> Tuple[float, float]:
    try:
        path = Chem.rdmolops.GetShortestPath(mol, start_idx, end_idx)
    except Exception:
        path = (start_idx, end_idx)
    if not path:
        path = (start_idx, end_idx)
    x = sum(conf.GetAtomPosition(idx).x for idx in path) / len(path)
    y = sum(conf.GetAtomPosition(idx).y for idx in path) / len(path)
    return x, y


def _virtual_arc_control_point(
    mol: Mol,
    conf,
    start_idx: int,
    end_idx: int,
    curvature: float = 0.55,
) -> Tuple[float, float]:
    start_pos = conf.GetAtomPosition(start_idx)
    end_pos = conf.GetAtomPosition(end_idx)
    mid_x = (start_pos.x + end_pos.x) / 2.0
    mid_y = (start_pos.y + end_pos.y) / 2.0
    dx = end_pos.x - start_pos.x
    dy = end_pos.y - start_pos.y
    length = math.hypot(dx, dy)
    if length <= 1e-6:
        return mid_x, mid_y

    perp_x = -dy / length
    perp_y = dx / length
    centroid_x, centroid_y = _path_centroid(mol, conf, start_idx, end_idx)
    away_x = mid_x - centroid_x
    away_y = mid_y - centroid_y
    if away_x * perp_x + away_y * perp_y < 0:
        perp_x = -perp_x
        perp_y = -perp_y
    return mid_x + perp_x * length * curvature, mid_y + perp_y * length * curvature


def _quadratic_bezier_point(
    start: Point2D,
    control: Point2D,
    end: Point2D,
    t: float,
) -> Point2D:
    one_minus_t = 1.0 - t
    x = one_minus_t * one_minus_t * start.x + 2 * one_minus_t * t * control.x + t * t * end.x
    y = one_minus_t * one_minus_t * start.y + 2 * one_minus_t * t * control.y + t * t * end.y
    return Point2D(x, y)


def _draw_quadratic_arc(
    drawer,
    start: Point2D,
    control: Point2D,
    end: Point2D,
    segments: int = 18,
) -> None:
    prev = start
    for step in range(1, segments + 1):
        t = step / segments
        current = _quadratic_bezier_point(start, control, end, t)
        drawer.DrawLine(prev, current)
        prev = current


def _draw_virtual_arcs_legacy(
    drawer,
    mol: Mol,
    virtual_arcs: List[Tuple[int, str, int, int]],
    arc_substituents: Dict[int, List[str]],
    drawing_config: DrawingConfig,
    bounds: Optional[List[float]],
) -> Optional[List[float]]:
    """Run the original named-virtualArc path without changing its SVG."""

    conf = mol.GetConformer()
    drawer.SetLineWidth(1)
    drawer.SetColour((0.1, 0.1, 0.1))
    for arc_idx, arc_name, start_idx, end_idx in virtual_arcs:
        if not (
            0 <= start_idx < mol.GetNumAtoms()
            and 0 <= end_idx < mol.GetNumAtoms()
        ):
            continue
        try:
            start_pos = conf.GetAtomPosition(start_idx)
            end_pos = conf.GetAtomPosition(end_idx)
            start_pt = Point2D(start_pos.x, start_pos.y)
            end_pt = Point2D(end_pos.x, end_pos.y)
            ctrl_x, ctrl_y = _virtual_arc_control_point(
                mol,
                conf,
                start_idx,
                end_idx,
            )
            ctrl_pt = Point2D(ctrl_x, ctrl_y)
            _draw_quadratic_arc(drawer, start_pt, ctrl_pt, end_pt)

            mid_x = (start_pos.x + end_pos.x) / 2.0
            mid_y = (start_pos.y + end_pos.y) / 2.0
            label_x = ctrl_x * 0.12 + mid_x * 0.88
            label_y = ctrl_y * 0.12 + mid_y * 0.88
            label_text = arc_name
            label_pt = Point2D(label_x, label_y)
            if label_text:
                drawer.DrawString(label_text, label_pt, 1)

            layout_points = [start_pt, ctrl_pt, end_pt]
            if label_text:
                layout_points.append(label_pt)
            for point in layout_points:
                draw_point = drawer.GetDrawCoords(point)
                bounds = _update_bounds(bounds, draw_point.x, draw_point.y)
            if label_text:
                draw_label = drawer.GetDrawCoords(label_pt)
                min_x, min_y, max_x, max_y = _estimate_text_bounds(
                    label_text,
                    draw_label.x,
                    draw_label.y,
                    font_px=float(drawing_config.visual.fixedFontSize),
                    margin=4.0,
                )
                bounds = _update_bounds(bounds, min_x, min_y)
                bounds = _update_bounds(bounds, max_x, max_y)

            outward_x = ctrl_x - mid_x
            outward_y = ctrl_y - mid_y
            outward_len = math.hypot(outward_x, outward_y)
            if outward_len <= 1e-6:
                outward_x, outward_y = 0.0, -1.0
                outward_len = 1.0
            outward_x /= outward_len
            outward_y /= outward_len

            labels = arc_substituents.get(arc_idx, [])
            for sub_idx, sub_label in enumerate(labels):
                lateral = (sub_idx - (len(labels) - 1) / 2.0) * 0.35
                lateral_x = -outward_y * lateral
                lateral_y = outward_x * lateral
                conn_start = Point2D(
                    label_x + lateral_x * 0.25,
                    label_y + lateral_y * 0.25,
                )
                conn_end = Point2D(
                    ctrl_x + outward_x * 0.65 + lateral_x,
                    ctrl_y + outward_y * 0.65 + lateral_y,
                )
                sub_label_pt = Point2D(
                    ctrl_x + outward_x * 0.9 + lateral_x,
                    ctrl_y + outward_y * 0.9 + lateral_y,
                )
                drawer.DrawLine(conn_start, conn_end)
                drawer.DrawString(sub_label, sub_label_pt, 1)
                for point in (conn_start, conn_end, sub_label_pt):
                    draw_point = drawer.GetDrawCoords(point)
                    bounds = _update_bounds(bounds, draw_point.x, draw_point.y)
                draw_label = drawer.GetDrawCoords(sub_label_pt)
                min_x, min_y, max_x, max_y = _estimate_text_bounds(
                    sub_label,
                    draw_label.x,
                    draw_label.y,
                    font_px=float(drawing_config.visual.fixedFontSize),
                    margin=4.0,
                )
                bounds = _update_bounds(bounds, min_x, min_y)
                bounds = _update_bounds(bounds, max_x, max_y)
        except Exception as exc:
            logger.debug("Failed to draw virtualArc %s: %s", arc_idx, exc)
    return bounds


def _cubic_bezier_point(
    start: Point2D,
    control_a: Point2D,
    control_b: Point2D,
    end: Point2D,
    t: float,
) -> Point2D:
    one_minus_t = 1.0 - t
    x = (
        one_minus_t**3 * start.x
        + 3 * one_minus_t * one_minus_t * t * control_a.x
        + 3 * one_minus_t * t * t * control_b.x
        + t**3 * end.x
    )
    y = (
        one_minus_t**3 * start.y
        + 3 * one_minus_t * one_minus_t * t * control_a.y
        + 3 * one_minus_t * t * t * control_b.y
        + t**3 * end.y
    )
    return Point2D(x, y)


def _draw_cubic_curve(
    drawer,
    start: Point2D,
    control_a: Point2D,
    control_b: Point2D,
    end: Point2D,
    segments: int = 14,
) -> List[Point2D]:
    points = [start]
    previous = start
    for step in range(1, segments + 1):
        current = _cubic_bezier_point(
            start,
            control_a,
            control_b,
            end,
            step / segments,
        )
        drawer.DrawLine(previous, current)
        points.append(current)
        previous = current
    return points


def _trim_virtual_arc_endpoint(
    mol: Mol,
    conf,
    atom_idx: int,
    control: Point2D,
) -> Point2D:
    """Keep an overlaid virtual arc clear of a visible endpoint label."""

    position = conf.GetAtomPosition(atom_idx)
    point = Point2D(position.x, position.y)
    atom = mol.GetAtomWithIdx(atom_idx)
    display_label = atom.GetProp("_displayLabel") if atom.HasProp("_displayLabel") else ""
    has_visible_label = atom.GetSymbol() != "C" or bool(display_label)
    if not has_visible_label:
        return point

    clearance = 0.32
    if atom.GetSymbol() in {"N", "O"} and atom.GetTotalNumHs() > 0:
        clearance = 0.40
    dx = control.x - point.x
    dy = control.y - point.y
    distance = math.hypot(dx, dy)
    if distance <= 1e-6:
        return point
    shift = min(clearance, distance * 0.35)
    return Point2D(point.x + dx / distance * shift, point.y + dy / distance * shift)


def _virtual_arc_layout(
    mol: Mol,
    conf,
    start_idx: int,
    end_idx: int,
    arc_name: str,
    substituents: List[str],
    font_px: float,
    draw_scale: float,
    label_spacing_scale: Optional[float] = None,
) -> Tuple[
    Point2D,
    Point2D,
    Point2D,
    Point2D,
    List[Tuple[Point2D, Point2D, Point2D]],
]:
    """Lay out one virtual arc, its inner name, and outward substituents."""

    ctrl_x, ctrl_y = _virtual_arc_control_point(mol, conf, start_idx, end_idx)
    control = Point2D(ctrl_x, ctrl_y)
    start = _trim_virtual_arc_endpoint(mol, conf, start_idx, control)
    end = _trim_virtual_arc_endpoint(mol, conf, end_idx, control)
    mid = Point2D((start.x + end.x) / 2.0, (start.y + end.y) / 2.0)
    apex = _quadratic_bezier_point(start, control, end, 0.5)

    outward_x = control.x - mid.x
    outward_y = control.y - mid.y
    outward_len = math.hypot(outward_x, outward_y)
    if outward_len <= 1e-6:
        outward_x, outward_y, outward_len = 0.0, -1.0, 1.0
    outward_x /= outward_len
    outward_y /= outward_len
    lateral_x, lateral_y = -outward_y, outward_x
    chord_x, chord_y = end.x - start.x, end.y - start.y
    if lateral_x * chord_x + lateral_y * chord_y < 0:
        lateral_x, lateral_y = -lateral_x, -lateral_y

    interior_radius = math.hypot(apex.x - mid.x, apex.y - mid.y)
    font_mol = font_px / max(draw_scale, 1e-6)
    desired_gap = min(0.65 * font_mol, 0.45 * interior_radius)
    name_distance = min(0.65 * interior_radius, interior_radius - desired_gap)
    name_distance = max(min(name_distance, 0.80 * interior_radius), 0.12)
    if interior_radius > 1e-6:
        name_distance = min(name_distance, 0.80 * interior_radius)
    if label_spacing_scale is not None:
        name_distance = max(0.12, name_distance / label_spacing_scale)

    name_visible = re.sub(r"<[^>]+>", "", arc_name)
    name_width = 0.60 * font_mol * max(len(name_visible), 1)
    name_height = font_mol

    def name_collision_penalty(candidate: Point2D) -> float:
        penalty = 0.0
        for atom in mol.GetAtoms():
            atom_idx = atom.GetIdx()
            atom_pos = conf.GetAtomPosition(atom_idx)
            display = atom.GetProp("_displayLabel") if atom.HasProp("_displayLabel") else ""
            if display:
                visible = re.sub(r"<[^>]+>", "", display)
            elif atom.GetSymbol() != "C":
                hydrogens = atom.GetTotalNumHs()
                visible = atom.GetSymbol()
                if hydrogens:
                    visible += "H" + (str(hydrogens) if hydrogens > 1 else "")
            else:
                visible = ""

            if visible:
                atom_width = 0.60 * font_mol * max(len(visible), 1)
                atom_height = font_mol
                if atom.GetSymbol() in {"N", "O"} and atom.GetTotalNumHs() > 0:
                    atom_width += 0.25
                    atom_height += 0.25
                overlap_x = (
                    (name_width + atom_width) / 2.0
                    + 0.14
                    - abs(candidate.x - atom_pos.x)
                )
                overlap_y = (
                    (name_height + atom_height) / 2.0
                    + 0.12
                    - abs(candidate.y - atom_pos.y)
                )
                if overlap_x > 0 and overlap_y > 0:
                    penalty += 100.0 + overlap_x * overlap_y
            elif math.hypot(candidate.x - atom_pos.x, candidate.y - atom_pos.y) < 0.34:
                penalty += 5.0
        return penalty

    lateral_cap = max(0.20, min(0.24 * math.hypot(chord_x, chord_y), 0.85))
    candidate_offsets = (0.0, -0.70 * font_mol, 0.70 * font_mol, -1.20 * font_mol, 1.20 * font_mol)
    candidates: List[Tuple[float, Point2D]] = []
    for radial_scale in (1.0, 1.08, 0.90, 1.15, 0.82):
        radial = min(name_distance * radial_scale, 0.80 * interior_radius)
        for raw_lateral in candidate_offsets:
            lateral = min(max(raw_lateral, -lateral_cap), lateral_cap)
            candidate = Point2D(
                mid.x + outward_x * radial + lateral_x * lateral,
                mid.y + outward_y * radial + lateral_y * lateral,
            )
            displacement = abs(radial - name_distance) + abs(lateral) * 0.2
            candidates.append((name_collision_penalty(candidate) * 1000.0 + displacement, candidate))
    name_point = min(candidates, key=lambda item: item[0])[1]

    visible_lengths = [
        max(len(re.sub(r"<[^>]+>", "", label)), 1) for label in substituents
    ]
    widest_label = max(visible_lengths, default=1)
    label_width_mol = max(
        0.5,
        0.60 * font_px * widest_label / max(draw_scale, 1e-6),
    )
    lateral_step = max(0.75, label_width_mol + 0.18)
    if label_spacing_scale is not None:
        lateral_step *= label_spacing_scale
    chord_length = max(math.hypot(chord_x, chord_y), 1e-6)
    connectors: List[Tuple[Point2D, Point2D, Point2D]] = []
    for sub_idx in range(len(substituents)):
        lateral = (sub_idx - (len(substituents) - 1) / 2.0) * lateral_step
        t = min(max(0.5 + lateral / chord_length, 0.22), 0.78)
        anchor = _quadratic_bezier_point(start, control, end, t)
        connector_end = Point2D(
            anchor.x
            + outward_x
            * 0.70
            * (label_spacing_scale if label_spacing_scale is not None else 1.0),
            anchor.y
            + outward_y
            * 0.70
            * (label_spacing_scale if label_spacing_scale is not None else 1.0),
        )
        label_point = Point2D(
            anchor.x
            + outward_x
            * (label_spacing_scale if label_spacing_scale is not None else 1.0),
            anchor.y
            + outward_y
            * (label_spacing_scale if label_spacing_scale is not None else 1.0),
        )
        connectors.append((anchor, connector_end, label_point))

    return start, control, end, name_point, connectors


def _serialize_svg(root: ET.Element) -> str:
    ET.register_namespace("", _SVG_NAMESPACE)
    svg_text = ET.tostring(root, encoding="unicode")
    # ElementTree may still emit ns-prefixed tags if another registration wins.
    return (
        svg_text.replace("<ns0:svg", "<svg")
        .replace("</ns0:svg>", "</svg>")
        .replace("ns0:", "")
        .replace(f'xmlns:ns0="{_SVG_NAMESPACE}"', f'xmlns="{_SVG_NAMESPACE}"')
    )


def _expand_svg_viewbox(svg_text: str, bounds: List[float], pad: float) -> str:
    try:
        root = ET.fromstring(svg_text)
    except Exception:
        return svg_text

    view_box = root.get("viewBox")
    if view_box:
        parts = view_box.split()
        if len(parts) == 4:
            try:
                vb_x, vb_y, vb_w, vb_h = map(float, parts)
            except Exception:
                return svg_text
        else:
            return svg_text
    else:
        width_attr = root.get("width", "")
        height_attr = root.get("height", "")
        width_match = re.search(r"[\d.]+", width_attr or "")
        height_match = re.search(r"[\d.]+", height_attr or "")
        if not width_match or not height_match:
            return svg_text
        vb_x, vb_y = 0.0, 0.0
        vb_w = float(width_match.group(0))
        vb_h = float(height_match.group(0))

    min_x, min_y, max_x, max_y = bounds
    new_min_x = min(vb_x, min_x - pad)
    new_min_y = min(vb_y, min_y - pad)
    new_max_x = max(vb_x + vb_w, max_x + pad)
    new_max_y = max(vb_y + vb_h, max_y + pad)
    new_w = new_max_x - new_min_x
    new_h = new_max_y - new_min_y
    root.set("viewBox", f"{new_min_x} {new_min_y} {new_w} {new_h}")
    root.set("width", f"{new_w}px")
    root.set("height", f"{new_h}px")
    root.set("overflow", "visible")

    for child in root.iter():
        tag = child.tag.rsplit("}", 1)[-1]
        if tag != "rect":
            continue
        style = child.get("style", "")
        if "fill:#FFFFFF" in style or "fill: #FFFFFF" in style:
            child.set("x", str(new_min_x))
            child.set("y", str(new_min_y))
            child.set("width", str(new_w))
            child.set("height", str(new_h))
            break

    for clip in root.iter():
        tag = clip.tag.rsplit("}", 1)[-1]
        if tag != "clipPath":
            continue
        for rect in list(clip):
            rect_tag = rect.tag.rsplit("}", 1)[-1]
            if rect_tag == "rect":
                rect.set("x", str(new_min_x))
                rect.set("y", str(new_min_y))
                rect.set("width", str(new_w))
                rect.set("height", str(new_h))
    return _serialize_svg(root)


class _DrawingPatterns:
    """More lenient regexes used only for drawing."""

    grp_content = re.compile(
        rf"(?P<{TextType.SYMBOL.value}>[^\[\?\'\"]*?)"
        + rf"(?P<{TextType.SCRIPT.value}>(\[\S+\])?)"
        + rf"(?P<{TextType.PRIME.value}>[\'\"]?)"
        + rf"(?P<{TextType.MULTIPLE.value}>(\?([a-z]|\d+|\d+-\d+)$)?)"
    )
    grp_pattern = re.compile(
        rf"({Tokens.atom_start}|{Tokens.circ_start}|{Tokens.dummy_start}|{Tokens.ring_start}|{Tokens.ring_start}{Tokens.circ_start}|{Tokens.ring_start}{Tokens.virtual_start})"
        + r"(\d+:\S*?)"
        + rf"({Tokens.atom_end}|{Tokens.circ_end}|{Tokens.dummy_end}|{Tokens.ring_end})"
    )
    trail_pattern = Patterns.trail_pattern


class _DrawingTranslator:
    """Tokenize / parse / render molecule from caption."""

    @classmethod
    def _preserve_stereochemistry(cls, original_smiles: str, mol: Mol) -> None:
        try:
            stereo_mol = Chem.MolFromSmiles(original_smiles, sanitize=False)
            if stereo_mol is None:
                return
            if stereo_mol.GetNumAtoms() != mol.GetNumAtoms():
                return
            for i in range(mol.GetNumAtoms()):
                original_atom = stereo_mol.GetAtomWithIdx(i)
                target_atom = mol.GetAtomWithIdx(i)
                if (
                    target_atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED
                    and original_atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
                ):
                    target_atom.SetChiralTag(original_atom.GetChiralTag())
                for bond in original_atom.GetBonds():
                    begin_idx = bond.GetBeginAtomIdx()
                    end_idx = bond.GetEndAtomIdx()
                    target_bond = mol.GetBondBetweenAtoms(begin_idx, end_idx)
                    if (
                        target_bond
                        and target_bond.GetStereo() == Chem.BondStereo.STEREONONE
                        and bond.GetStereo() != Chem.BondStereo.STEREONONE
                    ):
                        target_bond.SetStereo(bond.GetStereo())
        except Exception:
            pass

    @classmethod
    def parse_caption(
        cls, caption: str, return_mol: bool = False, error_msg: bool = False
    ) -> Optional[Tuple[Union[Mol, str], str, str]]:
        if Tokens.separator not in caption:
            if error_msg:
                logger.warning("No `%s` found in caption: %s", Tokens.separator, caption)
            return
        smi, trailing = caption.split(Tokens.separator, 1)

        if error_msg:
            RDLogger.EnableLog("rdApp.*")
        mol = Chem.MolFromSmiles(smi)
        RDLogger.DisableLog("rdApp.*")
        if mol is None:
            if error_msg:
                logger.warning("Invalid SMILES: %s", smi)
            return

        cls._preserve_stereochemistry(smi, mol)
        groups, ext = cls.parse_trailing(trailing)
        if return_mol:
            return mol, groups, ext
        return smi, groups, ext

    @classmethod
    def parse_trailing(cls, trailing: str):
        matched = re.match(_DrawingPatterns.trail_pattern, trailing)
        if matched is None:
            return "", ""
        content = matched.groupdict()
        return content.get("groups") or "", content.get("extension") or ""

    @classmethod
    def parse_groups(cls, seq: str) -> List[GroupDesc]:
        if seq == "":
            return []
        seq = re.sub(
            rf"{Tokens.substruct_start}.*?{Tokens.substruct_end}|"
            rf"{Tokens.sgroup_start}.*?{Tokens.sgroup_end}|"
            rf"{Tokens.virtual_start}.*?{Tokens.virtual_end}",
            "",
            seq,
            flags=re.DOTALL,
        )
        descriptions: List[GroupDesc] = []
        for grp_start, grp_content, _ in re.findall(_DrawingPatterns.grp_pattern, seq):
            parsed = cls.parse_group(grp_content)
            if parsed is None:
                continue
            idx, grp_text = parsed
            if grp_start in (Tokens.atom_start, Tokens.dummy_start):
                grp_desc = GroupDesc(id=AtomIndex(idx))
                if len(grp_text) == 0 or grp_start == Tokens.dummy_start:
                    grp_desc.is_dummy = True
            elif grp_start == Tokens.circ_start:
                grp_desc = GroupDesc(id=AtomIndex(idx), is_circle=True)
            elif grp_start in (
                f"{Tokens.ring_start}{Tokens.circ_start}",
                f"{Tokens.ring_start}{Tokens.virtual_start}",
            ):
                grp_desc = GroupDesc(id=RingIndex(idx, virtual=True))
            elif grp_start == Tokens.ring_start:
                grp_desc = GroupDesc(id=RingIndex(idx))
            else:
                continue
            symbol = grp_text.get(TextType.SYMBOL)
            if symbol:
                symbol = _format_symbol_with_subscripts(symbol)
            grp_desc.symbol = symbol
            grp_desc.script = grp_text.get(TextType.SCRIPT)
            grp_desc.prime = grp_text.get(TextType.PRIME)
            grp_desc.multiple = grp_text.get(TextType.MULTIPLE)
            descriptions.append(grp_desc)
        return descriptions

    @classmethod
    def parse_group(cls, group: str) -> Optional[Tuple[int, Dict[TextType, str]]]:
        items = group.split(":")
        if len(items) != 2:
            return
        idx, content = items
        if not idx.isdigit():
            return
        idx = int(idx)
        if content == Tokens.dummy or content == "":
            return idx, {}
        grp_text = cls.get_group_texts(content)
        if len(grp_text) == 0:
            return
        return idx, grp_text

    @classmethod
    def get_group_texts(cls, content: str) -> Dict[TextType, str]:
        texts: Dict[TextType, str] = {}
        bracket_matches = list(re.finditer(r"\[([^\]]+)\]", content))
        if bracket_matches:
            first_bracket_start = bracket_matches[0].start()
            last_bracket_end = bracket_matches[-1].end()
            symbol_part = content[:first_bracket_start] + content[last_bracket_end:]

            script_text = ""
            for i, match in enumerate(bracket_matches):
                script_text += match.group(1)
                if i < len(bracket_matches) - 1:
                    next_start = bracket_matches[i + 1].start()
                    script_text += content[match.end():next_start]
            texts[TextType.SCRIPT] = script_text
        else:
            symbol_part = content

        multiple_match = re.search(r"\?([a-z]|\d+|\d+-\d+)$", symbol_part)
        if multiple_match:
            texts[TextType.MULTIPLE] = multiple_match.group(1)
            symbol_part = re.sub(r"\?([a-z]|\d+|\d+-\d+)$", "", symbol_part)

        prime_match = re.search(r"[\'\"]$", symbol_part)
        if prime_match:
            texts[TextType.PRIME] = prime_match.group()
            symbol_part = symbol_part[:-1]

        if symbol_part:
            texts[TextType.SYMBOL] = symbol_part
        return texts

    @classmethod
    def reconstruct_mol(
        cls,
        mol: Union[Mol, str],
        groups: Optional[str] = None,
        write_to: Optional[str] = None,
        drawer_type: Literal["SVG", "PNG"] = "SVG",
        config: Optional[Union[Dict[str, Any], DrawingConfig]] = None,
    ) -> Optional[Union[str, bytes]]:
        if isinstance(config, DrawingConfig):
            drawing_config = config
        else:
            drawing_config = _validate_config(config)

        if drawer_type not in ("SVG", "PNG"):
            logger.error("Invalid drawer type `%s`", drawer_type)
            return
        extension = ""
        if isinstance(mol, Mol):
            assert groups is not None, "Empty group information"
        elif isinstance(mol, str):
            assert groups is None, "Excessive group information"
            parsed = cls.parse_caption(caption=mol, return_mol=True)
            if parsed is None:
                return
            mol, groups, extension = parsed
        else:
            raise TypeError(f"Invalid argument type: `{mol.__class__.__name__}`")

        grp_descriptions = cls.parse_groups(groups)
        candidate_single_ch2_repeats = _single_ch2_repeat_records(
            mol,
            grp_descriptions,
        )
        candidate_sgroup_records = _sgroup_records(groups)
        candidate_whole_sru_count = _sru_count(extension)
        repeat_setting = drawing_config.features.repeat_units
        repeat_units = (
            bool(
                candidate_single_ch2_repeats
                or candidate_sgroup_records
                or candidate_whole_sru_count is not None
            )
            if repeat_setting is None
            else repeat_setting
        )
        single_ch2_repeats = (
            candidate_single_ch2_repeats if repeat_units else []
        )
        single_ch2_indices = {record[0] for record in single_ch2_repeats}
        whole_sru_count = candidate_whole_sru_count if repeat_units else None
        virtual_arcs = _virtual_arcs(groups)
        academic_arc_setting = drawing_config.features.academic_virtual_arcs
        academic_virtual_arcs = (
            len(virtual_arcs) > 1 or any(not name for _idx, name, _start, _end in virtual_arcs)
            if academic_arc_setting is None
            else academic_arc_setting
        )
        endpoint_balls_enabled = drawing_config.features.endpoint_balls is not False
        whole_sru_atoms = (
            _whole_sru_atoms(mol, groups)
            if whole_sru_count is not None
            else None
        )
        whole_sru_ports = (
            _whole_sru_ports(mol, groups)
            if whole_sru_atoms is not None
            else None
        )
        whole_sru_endpoints = set(
            _whole_sru_endpoints(mol, groups) or []
        ) if whole_sru_atoms is not None else set()
        circ_indices = [
            int(desc.id)
            for desc in grp_descriptions
            if desc.is_circle and desc.symbol is not None
        ]

        ring_captions: List[str] = []
        ring_annotations: List[Tuple[int, str]] = []
        virtual_ring_annotations: List[Tuple[int, str]] = []
        endpoint_balls: Dict[int, str] = {}
        for desc in grp_descriptions:
            i = int(desc.id)
            label = None
            if isinstance(desc.id, AtomIndex):
                if not 0 <= i < mol.GetNumAtoms():
                    continue
                atom = mol.GetAtomWithIdx(i)
                ball_token = None
                if atom.GetSymbol() == "*":
                    if i in single_ch2_indices:
                        label = ""
                    else:
                        ball_token = (
                            _endpoint_ball_token(desc)
                            if endpoint_balls_enabled
                            else None
                        )
                        if ball_token is not None:
                            label = ""
                            endpoint_balls[i] = ball_token
                        elif desc.is_dummy:
                            label = "" if drawing_config.styling.use_modern_symbols else "dum"
                        elif desc.symbol is None:
                            continue
                        elif desc.is_circle:
                            label = (
                                desc.symbol
                                if drawing_config.styling.use_modern_symbols
                                else f"c{desc.symbol}"
                            )
                        else:
                            label = _format_group_label(desc)
                elif desc.symbol == Tokens.special_id:
                    # A non-ball special id is an opaque graphical group label,
                    # independent of the underlying atom element.
                    label = _format_group_label(desc)
                elif atom.GetSymbol() in ("C", "O"):
                    if desc.symbol is None and desc.multiple is not None:
                        label = f"({atom.GetSymbol()})<sub>{desc.multiple}</sub>"
                if label is not None:
                    atom.SetProp("_displayLabel", label)
                    if ball_token is None:
                        endpoint_balls.pop(i, None)
            elif isinstance(desc.id, RingIndex):
                if (
                    (not desc.id.virtual)
                    and (not 0 <= i < mol.GetRingInfo().NumRings())
                    or (desc.id.virtual)
                    and (i not in circ_indices)
                    or desc.symbol is None
                ):
                    continue
                if desc.symbol is not None:
                    if desc.multiple:
                        ring_label_text = "("
                        ring_label_text += desc.symbol
                        if desc.script:
                            ring_label_text += "<sub>" + desc.script + "</sub>"
                        if desc.prime:
                            ring_label_text += desc.prime
                        ring_label_text += ")"
                        ring_label_text += "<sub>" + desc.multiple + "</sub>"
                    else:
                        ring_label_text = desc.symbol
                        if desc.script:
                            ring_label_text += "<sub>" + desc.script + "</sub>"
                        if desc.prime:
                            ring_label_text += desc.prime
                    if desc.id.virtual:
                        virtual_ring_annotations.append((i, ring_label_text))
                    else:
                        ring_annotations.append((i, ring_label_text))

        # A legacy whole-SRU has no <d> records to suppress the two terminal
        # dummy labels.  Once its endpoints have passed the deliberately
        # strict, unambiguous inference above, depict them as ordinary
        # continuation bonds crossing the parentheses, just like explicit
        # E-SMILES 2.0 endpoints.
        if whole_sru_count is not None:
            for atom_idx in whole_sru_endpoints:
                mol.GetAtomWithIdx(atom_idx).SetProp("_displayLabel", "")

        mol.RemoveAllConformers()
        params = Chem.rdCoordGen.CoordGenParams()
        params.minimizerPrecision = params.sketcherBestPrecision
        Chem.rdCoordGen.AddCoords(mol, params)
        Chem.rdDepictor.StraightenDepiction(mol, 0)

        dummy_info: List[Tuple[int, int]] = []
        if drawing_config.features.dummy_atoms:
            for desc in grp_descriptions:
                if isinstance(desc.id, AtomIndex) and desc.is_dummy:
                    atom_idx = int(desc.id)
                    # A valid whole-SRU uses its terminal dummy atoms only as
                    # machine-readable endpoints.  The ordinary terminal
                    # bonds already cross the repeat parentheses, so adding
                    # wavy connection glyphs would be visually redundant.
                    if atom_idx in whole_sru_endpoints:
                        continue
                    if atom_idx < mol.GetNumAtoms():
                        atom = mol.GetAtomWithIdx(atom_idx)
                        connected_atoms = [
                            bond.GetOtherAtomIdx(atom_idx) for bond in atom.GetBonds()
                        ]
                        if connected_atoms:
                            dummy_info.append((atom_idx, connected_atoms[0]))

        palette_setting = drawing_config.styling.palette
        if isinstance(palette_setting, str):
            use_bw_palette = palette_setting == "bw"
        else:
            bw_prob = min(max(palette_setting, 0.0), 1.0)
            use_bw_palette = random.random() < bw_prob

        highlight_atoms: List[int] = []
        highlight_colors: Dict[int, Tuple[float, float, float]] = {}
        highlight_atom_radii: Dict[int, float] = {}
        highlight_bonds: List[int] = []
        highlight_bond_colors: Dict[int, Tuple[float, float, float]] = {}

        if drawing_config.styling.enhanced_contrast and not use_bw_palette:
            use_atom_hl = _resolve_probability(drawing_config.styling.highlight_atoms)
            use_bond_hl = _resolve_probability(drawing_config.styling.highlight_bonds)
            custom_colors_active = _resolve_probability(drawing_config.styling.custom_colors)
            if custom_colors_active:
                dummy_bond_color = _random_color()
                rgroup_bond_color = _random_color()
                rgroup_atom_color = _random_color()
            else:
                dummy_bond_color = (0.5, 0.7, 1.0)
                rgroup_bond_color = (0.5, 0.8, 0.5)
                rgroup_atom_color = (0.2, 0.7, 0.2)
            for desc in grp_descriptions:
                if isinstance(desc.id, AtomIndex):
                    atom_idx = int(desc.id)
                    if atom_idx < mol.GetNumAtoms():
                        atom = mol.GetAtomWithIdx(atom_idx)
                        if (
                            use_atom_hl
                            and desc.symbol is not None
                            and not desc.is_circle
                            and not desc.is_dummy
                            and atom_idx not in endpoint_balls
                        ):
                            highlight_atoms.append(atom_idx)
                            highlight_colors[atom_idx] = rgroup_atom_color
                        if use_bond_hl:
                            for bond in atom.GetBonds():
                                other_idx = bond.GetOtherAtomIdx(atom_idx)
                                bond_idx = mol.GetBondBetweenAtoms(atom_idx, other_idx).GetIdx()
                                if bond_idx in highlight_bonds:
                                    continue
                                highlight_bonds.append(bond_idx)
                                if desc.is_dummy:
                                    highlight_bond_colors[bond_idx] = dummy_bond_color
                                elif desc.symbol is not None:
                                    highlight_bond_colors[bond_idx] = rgroup_bond_color

        # RDKit sizes its dynamic canvas before the custom ball overlay is
        # painted.  A white, same-radius highlight reserves the required bounds
        # in both SVG and Cairo and is completely covered by the final ball.
        for atom_idx in endpoint_balls:
            if atom_idx not in highlight_atoms:
                highlight_atoms.append(atom_idx)
            highlight_colors[atom_idx] = (1.0, 1.0, 1.0)
            highlight_atom_radii[atom_idx] = _ENDPOINT_BALL_RADIUS

        if write_to is not None:
            drawer_type = "PNG"
        if drawer_type == "SVG":
            drawer = rdMolDraw2D.MolDraw2DSVG(width=-1, height=-1)
        else:
            drawer = rdMolDraw2D.MolDraw2DCairo(width=-1, height=-1)
        dopts = rdMolDraw2D.MolDrawOptions()
        if use_bw_palette:
            dopts.useBWAtomPalette()
        else:
            dopts.useCDKAtomPalette()
        visual_options = drawing_config.visual.model_dump()
        # This is consumed only by the custom virtual-arc overlay below.
        visual_options.pop("virtualArcLabelSpacingScale", None)
        # Assigning None to RDKit's string-valued fontFile option raises, and
        # omitting it is what preserves the historical built-in-font default.
        if visual_options.get("fontFile") is None:
            visual_options.pop("fontFile", None)
        visual_options["addAtomIndices"] = _resolve_probability(
            visual_options.get("addAtomIndices", False)
        )
        for k, v in visual_options.items():
            setattr(dopts, k, v)
        # Keep RDKit's close-contact diagnostic enabled.  The production
        # renderer treats any resulting red collision glyph as a rejected
        # image-caption pair; it must never be hidden by disabling this warning.
        dopts.flagCloseContactsDist = 3
        drawer.SetDrawOptions(dopts)
        rdMolDraw2D.PrepareAndDrawMolecule(
            drawer,
            mol,
            highlightAtoms=highlight_atoms,
            highlightAtomColors=highlight_colors,
            highlightAtomRadii=highlight_atom_radii,
            highlightBonds=highlight_bonds,
            highlightBondColors=highlight_bond_colors,
            legend=" | ".join(ring_captions),
        )

        ring_bounds: Optional[List[float]] = None

        if endpoint_balls:
            conf = mol.GetConformer()
            for atom_idx, token in endpoint_balls.items():
                try:
                    atom_pos = conf.GetAtomPosition(atom_idx)
                    center = Point2D(atom_pos.x, atom_pos.y)
                    _draw_endpoint_ball(
                        drawer,
                        center,
                        token,
                        drawing_config.visual.bondLineWidth,
                    )
                    if drawer_type == "SVG":
                        for corner in (
                            Point2D(
                                center.x - _ENDPOINT_BALL_RADIUS,
                                center.y - _ENDPOINT_BALL_RADIUS,
                            ),
                            Point2D(
                                center.x + _ENDPOINT_BALL_RADIUS,
                                center.y + _ENDPOINT_BALL_RADIUS,
                            ),
                        ):
                            draw_corner = drawer.GetDrawCoords(corner)
                            ring_bounds = _update_bounds(
                                ring_bounds,
                                draw_corner.x,
                                draw_corner.y,
                            )
                except Exception as e:
                    logger.debug(
                        "Failed to draw endpoint ball `%s` at atom %s: %s",
                        token,
                        atom_idx,
                        e,
                    )

        if drawing_config.features.dummy_atoms and dummy_info and drawer_type == "SVG":
            conf = mol.GetConformer()
            drawer.SetLineWidth(drawing_config.visual.bondLineWidth)
            for dummy_idx, connected_idx in dummy_info:
                try:
                    dummy_pos = conf.GetAtomPosition(dummy_idx)
                    connected_pos = conf.GetAtomPosition(connected_idx)
                    dx = connected_pos.x - dummy_pos.x
                    dy = connected_pos.y - dummy_pos.y
                    length = math.sqrt(dx * dx + dy * dy)
                    if length > 0:
                        dx /= length
                        dy /= length
                        perp_x = dy
                        perp_y = -dx
                        extension_length = 0.35
                        start_x = dummy_pos.x - perp_x * extension_length
                        start_y = dummy_pos.y - perp_y * extension_length
                        end_x = dummy_pos.x + perp_x * extension_length
                        end_y = dummy_pos.y + perp_y * extension_length
                        wavy_color = (0.0, 0.0, 0.0)
                        drawer.DrawWavyLine(
                            Point2D(start_x, start_y),
                            Point2D(end_x, end_y),
                            wavy_color,
                            wavy_color,
                            nSegments=8,
                            vertOffset=0.08,
                        )
                except Exception as e:
                    logger.debug("Failed to draw wavy line for dummy atom %s: %s", dummy_idx, e)

        if drawing_config.features.circled_r_groups and drawer_type == "SVG":
            conf = mol.GetConformer()
            for desc in grp_descriptions:
                if isinstance(desc.id, AtomIndex) and desc.is_circle:
                    atom_idx = int(desc.id)
                    if atom_idx < mol.GetNumAtoms():
                        try:
                            atom_pos = conf.GetAtomPosition(atom_idx)
                            r = 0.5
                            drawer.SetColour((0, 0, 0))
                            drawer.SetFillPolys(False)
                            drawer.SetLineWidth(drawing_config.visual.bondLineWidth)
                            drawer.DrawEllipse(
                                Point2D(atom_pos.x - r, atom_pos.y - r),
                                Point2D(atom_pos.x + r, atom_pos.y + r),
                            )
                        except Exception as e:
                            logger.debug(
                                "Failed to draw circle for circled R-group at atom %s: %s",
                                atom_idx,
                                e,
                            )

        if (
            drawing_config.features.ring_annotations
            and virtual_ring_annotations
            and drawer_type == "SVG"
        ):
            conf = mol.GetConformer()
            by_atom: Dict[int, List[str]] = {}
            for atom_idx, label_text in virtual_ring_annotations:
                by_atom.setdefault(atom_idx, []).append(label_text)

            for atom_idx, labels in by_atom.items():
                if not 0 <= atom_idx < mol.GetNumAtoms():
                    continue
                try:
                    atom_pos = conf.GetAtomPosition(atom_idx)
                    # Spread labels around the abstract ring placeholder.
                    start_angle = -math.pi / 2
                    if len(labels) > 1:
                        start_angle -= (len(labels) - 1) * math.pi / 8
                    for label_idx, label_text in enumerate(labels):
                        angle = start_angle + label_idx * math.pi / 4
                        dir_x = math.cos(angle)
                        dir_y = math.sin(angle)
                        line_start = Point2D(
                            atom_pos.x + dir_x * 0.45,
                            atom_pos.y + dir_y * 0.45,
                        )
                        line_end = Point2D(
                            atom_pos.x + dir_x * 0.9,
                            atom_pos.y + dir_y * 0.9,
                        )
                        label_pt = Point2D(
                            atom_pos.x + dir_x * 1.15,
                            atom_pos.y + dir_y * 1.15,
                        )
                        drawer.SetColour((0, 0, 0))
                        drawer.SetLineWidth(1)
                        drawer.DrawLine(line_start, line_end)
                        drawer.DrawString(label_text, label_pt, 1)

                        for pt in (line_start, line_end, label_pt):
                            draw_pt = drawer.GetDrawCoords(pt)
                            ring_bounds = _update_bounds(ring_bounds, draw_pt.x, draw_pt.y)
                        draw_label = drawer.GetDrawCoords(label_pt)
                        min_x, min_y, max_x, max_y = _estimate_text_bounds(
                            label_text,
                            draw_label.x,
                            draw_label.y,
                            font_px=float(drawing_config.visual.fixedFontSize),
                            margin=4.0,
                        )
                        ring_bounds = _update_bounds(ring_bounds, min_x, min_y)
                        ring_bounds = _update_bounds(ring_bounds, max_x, max_y)
                except Exception as e:
                    logger.debug("Failed to draw virtual ring annotations at atom %s: %s", atom_idx, e)

        if (
            drawing_config.features.ring_annotations
            and ring_annotations
            and drawer_type == "SVG"
        ):
            conf = mol.GetConformer()
            ring_info = mol.GetRingInfo()
            atom_rings = ring_info.AtomRings()
            bond_rings = ring_info.BondRings()
            drawer.SetLineWidth(1)
            ring_annotations_by_idx: Dict[int, List[str]] = {}
            for ring_idx, label_text in ring_annotations:
                ring_annotations_by_idx.setdefault(ring_idx, []).append(label_text)
            for ring_idx, label_texts in ring_annotations_by_idx.items():
                if not (0 <= ring_idx < len(atom_rings)):
                    continue
                ring_atoms = atom_rings[ring_idx]
                if not ring_atoms:
                    continue
                center_x = sum(conf.GetAtomPosition(idx).x for idx in ring_atoms) / len(ring_atoms)
                center_y = sum(conf.GetAtomPosition(idx).y for idx in ring_atoms) / len(ring_atoms)

                direction_x = direction_y = None
                ring_radius = None
                midpoint_distance = None
                bond_mid_data: List[Tuple[float, float, float]] = []
                bond_ids = bond_rings[ring_idx] if 0 <= ring_idx < len(bond_rings) else []
                for bond_idx in bond_ids:
                    bond = mol.GetBondWithIdx(bond_idx)
                    pos1 = conf.GetAtomPosition(bond.GetBeginAtomIdx())
                    pos2 = conf.GetAtomPosition(bond.GetEndAtomIdx())
                    mid_x = (pos1.x + pos2.x) / 2.0
                    mid_y = (pos1.y + pos2.y) / 2.0
                    dx = mid_x - center_x
                    dy = mid_y - center_y
                    dist = math.hypot(dx, dy)
                    if dist > 1e-6:
                        bond_mid_data.append((dist, dx / dist, dy / dist))
                if bond_mid_data:
                    midpoint_distance, direction_x, direction_y = max(
                        bond_mid_data, key=lambda item: item[0]
                    )
                    ring_radius = midpoint_distance
                else:
                    for atom_idx in ring_atoms:
                        pos = conf.GetAtomPosition(atom_idx)
                        dx = pos.x - center_x
                        dy = pos.y - center_y
                        dist = math.hypot(dx, dy)
                        if dist > 1e-6:
                            direction_x = dx / dist
                            direction_y = dy / dist
                            ring_radius = dist
                            midpoint_distance = dist
                            break

                if ring_radius is None or direction_x is None:
                    continue
                if midpoint_distance is None:
                    midpoint_distance = ring_radius
                extension = drawing_config.styling.ring_connector_extension
                total_length = midpoint_distance + extension

                base_angle = math.atan2(direction_y, direction_x)
                spread = math.radians(26)
                label_offset = 0.25
                for label_idx, label_text in enumerate(label_texts):
                    angle_offset = (label_idx - (len(label_texts) - 1) / 2.0) * spread
                    angle = base_angle + angle_offset
                    label_dir_x = math.cos(angle)
                    label_dir_y = math.sin(angle)

                    drawer.SetColour((0, 0, 0))
                    if drawing_config.styling.ring_connector_style == "solid":
                        drawer.DrawLine(
                            Point2D(center_x, center_y),
                            Point2D(
                                center_x + label_dir_x * total_length,
                                center_y + label_dir_y * total_length,
                            ),
                        )
                    else:
                        dash_length = 0.15
                        gap_length = 0.1
                        offset = 0.0
                        while offset < total_length:
                            end_offset = min(offset + dash_length, total_length)
                            if end_offset > offset:
                                drawer.DrawLine(
                                    Point2D(
                                        center_x + label_dir_x * offset,
                                        center_y + label_dir_y * offset,
                                    ),
                                    Point2D(
                                        center_x + label_dir_x * end_offset,
                                        center_y + label_dir_y * end_offset,
                                    ),
                                )
                            offset = end_offset + gap_length

                    label_pt = Point2D(
                        center_x + label_dir_x * (total_length + label_offset),
                        center_y + label_dir_y * (total_length + label_offset),
                    )
                    drawer.SetColour((0, 0, 0))
                    drawer.DrawString(label_text, label_pt, 1)

                    draw_center = drawer.GetDrawCoords(Point2D(center_x, center_y))
                    draw_label = drawer.GetDrawCoords(label_pt)
                    ring_bounds = _update_bounds(ring_bounds, draw_center.x, draw_center.y)
                    ring_bounds = _update_bounds(ring_bounds, draw_label.x, draw_label.y)
                    mol_end = Point2D(
                        center_x + label_dir_x * total_length,
                        center_y + label_dir_y * total_length,
                    )
                    draw_end = drawer.GetDrawCoords(mol_end)
                    ring_bounds = _update_bounds(ring_bounds, draw_end.x, draw_end.y)

                    font_px = float(drawing_config.visual.fixedFontSize)
                    min_x, min_y, max_x, max_y = _estimate_text_bounds(
                        label_text,
                        draw_label.x,
                        draw_label.y,
                        font_px=font_px,
                        margin=4.0,
                    )
                    ring_bounds = _update_bounds(ring_bounds, min_x, min_y)
                    ring_bounds = _update_bounds(ring_bounds, max_x, max_y)

        rendered_sgroups: Optional[set[str]] = None
        if repeat_units:
            rendered_sgroups = set()
            for body, ports, count in candidate_sgroup_records:
                repeat_atoms = _sgroup_repeat_atoms(mol, ports)
                if repeat_atoms is None:
                    continue
                ring_bounds = _draw_sgroup_boundary_parentheses(
                    drawer,
                    mol,
                    ports,
                    count,
                    drawing_config.visual.bondLineWidth,
                    float(drawing_config.visual.fixedFontSize),
                    ring_bounds,
                    count != "n" or drawing_config.styling.show_default_sru_count,
                    drawing_config.styling.sgroup_bracket_style,
                )
                rendered_sgroups.add(body)

        for _atom_idx, ports, count in single_ch2_repeats:
            ring_bounds = _draw_sgroup_boundary_parentheses(
                drawer,
                mol,
                ports,
                count,
                drawing_config.visual.bondLineWidth,
                float(drawing_config.visual.fixedFontSize),
                ring_bounds,
                True,
                "round",
            )

        if whole_sru_count is not None and whole_sru_ports is not None:
            ring_bounds = _draw_sgroup_boundary_parentheses(
                drawer,
                mol,
                whole_sru_ports,
                whole_sru_count,
                drawing_config.visual.bondLineWidth,
                float(drawing_config.visual.fixedFontSize),
                ring_bounds,
                whole_sru_count != "n"
                or drawing_config.styling.show_default_sru_count,
                "round",
            )

        if virtual_arcs and not academic_virtual_arcs:
            ring_bounds = _draw_virtual_arcs_legacy(
                drawer,
                mol,
                virtual_arcs,
                _virtual_arc_substituents(groups),
                drawing_config,
                ring_bounds,
            )
        if virtual_arcs and academic_virtual_arcs:
            conf = mol.GetConformer()
            arc_substituents = _virtual_arc_substituents(groups)
            drawer.SetLineWidth(drawing_config.visual.bondLineWidth)
            drawer.SetColour((0.1, 0.1, 0.1))
            draw_origin = drawer.GetDrawCoords(Point2D(0.0, 0.0))
            draw_unit = drawer.GetDrawCoords(Point2D(1.0, 0.0))
            draw_scale = max(
                math.hypot(draw_unit.x - draw_origin.x, draw_unit.y - draw_origin.y),
                1e-6,
            )
            for arc_idx, arc_name, start_idx, end_idx in virtual_arcs:
                if not (0 <= start_idx < mol.GetNumAtoms() and 0 <= end_idx < mol.GetNumAtoms()):
                    continue
                try:
                    substituents = arc_substituents.get(arc_idx, [])
                    start_pt, ctrl_pt, end_pt, label_pt, connector_layout = (
                        _virtual_arc_layout(
                            mol,
                            conf,
                            start_idx,
                            end_idx,
                            arc_name,
                            substituents,
                            float(drawing_config.visual.fixedFontSize),
                            draw_scale,
                            drawing_config.visual.virtualArcLabelSpacingScale,
                        )
                    )
                    _draw_quadratic_arc(drawer, start_pt, ctrl_pt, end_pt)

                    label_text = arc_name
                    if label_text:
                        drawer.DrawString(label_text, label_pt, 1)

                    layout_points = [start_pt, ctrl_pt, end_pt]
                    if label_text:
                        layout_points.append(label_pt)
                    for pt in layout_points:
                        draw_pt = drawer.GetDrawCoords(pt)
                        ring_bounds = _update_bounds(ring_bounds, draw_pt.x, draw_pt.y)
                    if label_text:
                        draw_label = drawer.GetDrawCoords(label_pt)
                        min_x, min_y, max_x, max_y = _estimate_text_bounds(
                            label_text,
                            draw_label.x,
                            draw_label.y,
                            font_px=float(drawing_config.visual.fixedFontSize),
                            margin=4.0,
                        )
                        ring_bounds = _update_bounds(ring_bounds, min_x, min_y)
                        ring_bounds = _update_bounds(ring_bounds, max_x, max_y)

                    for sub_label, (conn_start, conn_end, sub_label_pt) in zip(
                        substituents,
                        connector_layout,
                    ):
                        drawer.DrawLine(conn_start, conn_end)
                        drawer.DrawString(sub_label, sub_label_pt, 1)
                        for pt in (conn_start, conn_end, sub_label_pt):
                            draw_pt = drawer.GetDrawCoords(pt)
                            ring_bounds = _update_bounds(ring_bounds, draw_pt.x, draw_pt.y)
                        draw_label = drawer.GetDrawCoords(sub_label_pt)
                        min_x, min_y, max_x, max_y = _estimate_text_bounds(
                            sub_label,
                            draw_label.x,
                            draw_label.y,
                            font_px=float(drawing_config.visual.fixedFontSize),
                            margin=4.0,
                        )
                        ring_bounds = _update_bounds(ring_bounds, min_x, min_y)
                        ring_bounds = _update_bounds(ring_bounds, max_x, max_y)
                except Exception as e:
                    logger.debug("Failed to draw virtualArc %s: %s", arc_idx, e)

        precompat_notes = _precompat_annotations(groups, rendered_sgroups)
        if precompat_notes:
            conf = mol.GetConformer()
            xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
            ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
            note_x = min(xs) if xs else 0.0
            note_y = (max(ys) if ys else 0.0) + 0.8
            drawer.SetColour((0.1, 0.1, 0.1))
            for idx, note in enumerate(precompat_notes):
                note_pt = Point2D(note_x, note_y + idx * 0.45)
                drawer.DrawString(note, note_pt, 1)
                draw_note = drawer.GetDrawCoords(note_pt)
                ring_bounds = _update_bounds(ring_bounds, draw_note.x, draw_note.y)
                min_x, min_y, max_x, max_y = _estimate_text_bounds(
                    note,
                    draw_note.x,
                    draw_note.y,
                    font_px=float(drawing_config.visual.fixedFontSize),
                    margin=4.0,
                )
                ring_bounds = _update_bounds(ring_bounds, min_x, min_y)
                ring_bounds = _update_bounds(ring_bounds, max_x, max_y)

        drawer.FinishDrawing()
        if write_to is not None:
            Path(write_to).parent.mkdir(parents=True, exist_ok=True)
            drawer.WriteDrawingText(write_to)
            return
        svg_text = drawer.GetDrawingText()
        if ring_bounds is not None:
            svg_text = _expand_svg_viewbox(svg_text, ring_bounds, pad=6.0)
        return svg_text


def _needs_svg_overlay(smi: str, config: DrawingConfig) -> bool:
    """Return whether PNG must pass through the completed SVG overlay."""

    repeat_syntax = bool(
        "<g>" in smi
        or "|Sg:" in smi
        or re.search(r"<a>\d+:CH2\?(?:[a-z]|\d+|\d+-\d+)</a>", smi)
    )
    repeat_setting = config.features.repeat_units
    repeat_overlay = repeat_setting is True or (
        repeat_setting is None and repeat_syntax
    )

    ball_pattern = "|".join(re.escape(token) for token in _ENDPOINT_BALL_STYLES)
    ball_syntax = re.search(rf"<id>\[(?:{ball_pattern})\]", smi) is not None
    ball_setting = config.features.endpoint_balls
    ball_overlay = ball_setting is True or (
        ball_setting is None and ball_syntax
    )

    arc_names = re.findall(
        r"<v>\d+:([^:]*):\[\d+:\d+\]</v>",
        smi,
        flags=re.DOTALL,
    )
    academic_arc_syntax = len(arc_names) > 1 or any(not name for name in arc_names)
    arc_setting = config.features.academic_virtual_arcs
    arc_overlay = arc_setting is True or (
        arc_setting is None and academic_arc_syntax
    )
    return repeat_overlay or ball_overlay or arc_overlay


def draw(
    smi: str,
    config: Optional[Union[Dict[str, Any], DrawingConfig]] = None,
    output_format: Literal["svg", "png"] = "svg",
) -> Union[str, bytes]:
    """Draw a SMILES / E-SMILES string and return SVG text or PNG bytes."""
    drawing_config = (
        config if isinstance(config, DrawingConfig) else _validate_config(config)
    )
    requested_format = output_format.lower()
    svg_to_png = None
    if requested_format == "png" and _needs_svg_overlay(smi, drawing_config):
        try:
            # CairoSVG is a required project dependency.  Rendering PNG from
            # the completed SVG keeps custom Markush connectors, label paths,
            # margins, and stereobonds identical across the two public formats.
            from cairosvg import svg2png as svg_to_png
        except (ImportError, OSError):
            # Keep the historical RDKit-Cairo fallback for environments that
            # import the source tree without installing its declared extras.
            svg_to_png = None
        if svg_to_png is None and Tokens.separator in smi and any(
            token in smi
            for token in ("<s>", "<g>", "<v>", "|Sg:", "<r>", "<d>", "<c>")
        ):
            raise RuntimeError(
                "CairoSVG is required for PNG output containing SVG-only "
                "E-SMILES overlays; the direct Cairo fallback cannot render "
                "or expand them safely"
            )

    drawer_type: Literal["SVG", "PNG"] = (
        "SVG" if requested_format == "svg" or svg_to_png is not None else "PNG"
    )
    if Tokens.separator in smi:
        drawing = _DrawingTranslator.reconstruct_mol(
            mol=smi,
            drawer_type=drawer_type,
            config=drawing_config,
        )
    else:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smi}")
        drawing = _DrawingTranslator.reconstruct_mol(
            mol=mol,
            groups="",
            drawer_type=drawer_type,
            config=drawing_config,
        )
    if drawing is None:
        raise RuntimeError(f"Failed to draw molecule {smi}")
    if svg_to_png is not None:
        if not isinstance(drawing, str):
            raise RuntimeError("SVG renderer returned non-text output")
        return svg_to_png(bytestring=drawing.encode("utf-8"))
    return drawing
