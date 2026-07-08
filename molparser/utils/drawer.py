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


class VisualConfig(BaseModel):
    padding: float = Field(default=0.1, ge=0.0, le=1.0)
    additionalAtomLabelPadding: float = Field(default=0.05, ge=0.0, le=0.5)
    fixedFontSize: int = Field(default=10, ge=1, le=100)
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


class FeaturesConfig(BaseModel):
    dummy_atoms: bool = Field(default=True)
    circled_r_groups: bool = Field(default=True)
    ring_annotations: bool = Field(default=True)


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
    label_text: str, x: float, y: float, font_px: float, margin: float
) -> Tuple[float, float, float, float]:
    visible = re.sub(r"<[^>]+>", "", label_text)
    text_len = max(len(visible), 1)
    text_w = font_px * 0.6 * text_len + 2 * margin
    text_h = font_px + 2 * margin
    return x - text_w, y - text_h, x + text_w, y + text_h


def _format_symbol_with_subscripts(symbol: str) -> str:
    if "<" in symbol and ">" in symbol:
        return symbol
    return re.sub(r"(\d+)", r"<sub>\1</sub>", symbol)


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


def _truncate_annotation(text: str, limit: int = 64) -> str:
    compact = re.sub(r"\s+", " ", text).strip()
    return compact if len(compact) <= limit else compact[: limit - 3] + "..."


def _precompat_annotations(groups: str) -> List[str]:
    annotations: List[str] = []
    for tag, label in (("s", "substructure"), ("g", "s-group")):
        for match in re.finditer(rf"<{tag}>(.*?)</{tag}>", groups, re.DOTALL):
            annotations.append(f"{label}: {_truncate_annotation(match.group(1))}")
    return annotations


def _virtual_arcs(groups: str) -> List[Tuple[int, str, int, int]]:
    arcs: List[Tuple[int, str, int, int]] = []
    pattern = r"<v>(?P<idx>\d+):(?P<name>[^:]+):\[(?P<start>\d+):(?P<end>\d+)\]</v>"
    for match in re.finditer(pattern, groups, re.DOTALL):
        arcs.append(
            (
                int(match.group("idx")),
                match.group("name"),
                int(match.group("start")),
                int(match.group("end")),
            )
        )
    return arcs


def _virtual_arc_substituents(groups: str) -> Dict[int, List[str]]:
    substituents: Dict[int, List[str]] = {}
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
        substituents.setdefault(idx, []).append(label)
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
        if isinstance(mol, Mol):
            assert groups is not None, "Empty group information"
        elif isinstance(mol, str):
            assert groups is None, "Excessive group information"
            parsed = cls.parse_caption(caption=mol, return_mol=True)
            if parsed is None:
                return
            mol, groups, _ = parsed
        else:
            raise TypeError(f"Invalid argument type: `{mol.__class__.__name__}`")

        grp_descriptions = cls.parse_groups(groups)
        circ_indices = [
            int(desc.id)
            for desc in grp_descriptions
            if desc.is_circle and desc.symbol is not None
        ]

        ring_captions: List[str] = []
        ring_annotations: List[Tuple[int, str]] = []
        virtual_ring_annotations: List[Tuple[int, str]] = []
        for desc in grp_descriptions:
            i = int(desc.id)
            label = None
            if isinstance(desc.id, AtomIndex):
                if not 0 <= i < mol.GetNumAtoms():
                    continue
                atom = mol.GetAtomWithIdx(i)
                if atom.GetSymbol() == "*":
                    if desc.is_dummy:
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
                elif atom.GetSymbol() in ("C", "O"):
                    if desc.symbol == Tokens.special_id:
                        label = _format_group_label(desc)
                    elif desc.symbol is None and desc.multiple is not None:
                        label = f"({atom.GetSymbol()})<sub>{desc.multiple}</sub>"
                if label is not None:
                    atom.SetProp("_displayLabel", label)
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
        visual_options["addAtomIndices"] = _resolve_probability(
            visual_options.get("addAtomIndices", False)
        )
        for k, v in visual_options.items():
            setattr(dopts, k, v)
        drawer.SetDrawOptions(dopts)
        rdMolDraw2D.PrepareAndDrawMolecule(
            drawer,
            mol,
            highlightAtoms=highlight_atoms,
            highlightAtomColors=highlight_colors,
            highlightBonds=highlight_bonds,
            highlightBondColors=highlight_bond_colors,
            legend=" | ".join(ring_captions),
        )

        ring_bounds: Optional[List[float]] = None

        if drawing_config.features.dummy_atoms and dummy_info and drawer_type == "SVG":
            conf = mol.GetConformer()
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
                            drawer.SetLineWidth(1)
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

        virtual_arcs = _virtual_arcs(groups)
        if virtual_arcs:
            conf = mol.GetConformer()
            arc_substituents = _virtual_arc_substituents(groups)
            drawer.SetLineWidth(1)
            drawer.SetColour((0.1, 0.1, 0.1))
            for arc_idx, arc_name, start_idx, end_idx in virtual_arcs:
                if not (0 <= start_idx < mol.GetNumAtoms() and 0 <= end_idx < mol.GetNumAtoms()):
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
                    drawer.DrawString(label_text, label_pt, 1)

                    for pt in (start_pt, ctrl_pt, end_pt, label_pt):
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

                    outward_x = ctrl_x - mid_x
                    outward_y = ctrl_y - mid_y
                    outward_len = math.hypot(outward_x, outward_y)
                    if outward_len <= 1e-6:
                        outward_x, outward_y = 0.0, -1.0
                        outward_len = 1.0
                    outward_x /= outward_len
                    outward_y /= outward_len

                    for sub_idx, sub_label in enumerate(arc_substituents.get(arc_idx, [])):
                        lateral = (sub_idx - (len(arc_substituents.get(arc_idx, [])) - 1) / 2.0) * 0.35
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

        precompat_notes = _precompat_annotations(groups)
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


def draw(
    smi: str,
    config: Optional[Union[Dict[str, Any], DrawingConfig]] = None,
    output_format: Literal["svg", "png"] = "svg",
) -> Union[str, bytes]:
    """Draw a SMILES / E-SMILES string and return SVG text or PNG bytes."""
    drawer_type: Literal["SVG", "PNG"] = "SVG" if output_format.lower() == "svg" else "PNG"
    if Tokens.separator in smi:
        drawing = _DrawingTranslator.reconstruct_mol(
            mol=smi,
            drawer_type=drawer_type,
            config=config,
        )
    else:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smi}")
        drawing = _DrawingTranslator.reconstruct_mol(
            mol=mol,
            groups="",
            drawer_type=drawer_type,
            config=config,
        )
    if drawing is None:
        raise RuntimeError(f"Failed to draw molecule {smi}")
    return drawing
