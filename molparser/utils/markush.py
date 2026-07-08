"""Markush substituent expansion for E-SMILES captions."""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from itertools import combinations
from typing import Iterable

from rdkit import Chem, RDLogger

try:
    from . import chem_utils
    from .translator import AtomIndex, GroupDesc, RingIndex, Tokens, Translator
except ImportError:  # Support running from package directory as working directory.
    import chem_utils
    from translator import AtomIndex, GroupDesc, RingIndex, Tokens, Translator


DefinitionValue = int | str | Sequence[str]
_PRECOMPAT_RECORD_PATTERN = re.compile(
    r"<s>.*?</s>|<g>.*?</g>|<r><v>\d+:.+?</r>|<v>.*?</v>",
    re.DOTALL,
)
_SUBSTRUCT_RECORD_PATTERN = re.compile(r"^<s>(?P<body>.*)</s>$", re.DOTALL)


def _normalize_label(label: str) -> str:
    return label.replace("[", "").replace("]", "").strip()


def _definition_lookup(definitions: Mapping[str, DefinitionValue]) -> dict[str, DefinitionValue]:
    lookup: dict[str, DefinitionValue] = {}
    for key, value in definitions.items():
        clean_key = str(key).strip()
        lookup[clean_key] = value
        lookup[_normalize_label(clean_key)] = value
    return lookup


def _group_labels(desc: GroupDesc) -> list[str]:
    if not desc.symbol:
        return []
    labels = [desc.symbol]
    if desc.script:
        labels = [f"{desc.symbol}[{desc.script}]", f"{desc.symbol}{desc.script}"]
    return labels


def _as_values(value: DefinitionValue) -> list[str]:
    if isinstance(value, (int, str)):
        return [str(value)]
    return [str(item) for item in value]


def _resolve_group_smiles(
    desc: GroupDesc,
    definitions: Mapping[str, DefinitionValue],
) -> list[str]:
    definition_lookup = _definition_lookup(definitions)
    for label in _group_labels(desc):
        value = definition_lookup.get(label)
        if value is None:
            value = definition_lookup.get(_normalize_label(label))
        if value is not None:
            return [_resolve_fragment_smiles(v) for v in _as_values(value)]

    lookup_symbol = desc.symbol or ""
    if desc.script:
        lookup_symbol += desc.script
    src = chem_utils.get_abbrev_smi().get(lookup_symbol)
    if src is not None:
        return [src]
    return []


def _resolve_fragment_smiles(value: str) -> str:
    value = str(value).strip()
    return chem_utils.get_abbrev_smi().get(value, value)


def _repeat_counts(desc: GroupDesc, site_count: int) -> list[int]:
    if not desc.multiple:
        return [1]
    if desc.multiple.isdigit():
        return [int(desc.multiple)]
    if desc.multiple == "n":
        return list(range(1, site_count + 1))
    if "-" in desc.multiple:
        start, end = desc.multiple.split("-", 1)
        if start.isdigit() and end.isdigit():
            return list(range(int(start), int(end) + 1))
    raise ValueError(f"Unsupported multiplicity: ?{desc.multiple}")


def _resolve_repetition_count(
    desc: GroupDesc,
    definitions: Mapping[str, DefinitionValue],
) -> int | None:
    if not desc.multiple:
        return None
    if desc.multiple.isdigit():
        return int(desc.multiple)

    value = _definition_lookup(definitions).get(desc.multiple)
    if isinstance(value, Sequence) and not isinstance(value, str):
        raise ValueError(f"Multiplicity ?{desc.multiple} must resolve to one integer")
    if value is None:
        return None
    try:
        count = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Multiplicity ?{desc.multiple} must resolve to one integer") from exc
    if count < 1:
        raise ValueError(f"Multiplicity ?{desc.multiple} must be positive")
    return count


def _is_carbon_chain_repeat(desc: GroupDesc) -> bool:
    return (desc.symbol == "CH2") or (desc.symbol == "CH" and desc.script == "2")


def _is_special_id_label(desc: GroupDesc) -> bool:
    return desc.symbol == Tokens.special_id


def _find_atom_by_source_index(mol: Chem.rdchem.Mol, source_idx: int) -> int | None:
    map_num = source_idx + 1
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() == map_num:
            return atom.GetIdx()
    return None


def _source_attachment(src_mol: Chem.rdchem.Mol) -> tuple[int, int | None]:
    dummy_atoms = [atom for atom in src_mol.GetAtoms() if atom.GetSymbol() == "*"]
    if not dummy_atoms:
        return 0, None
    if len(dummy_atoms) != 1:
        raise ValueError("Substituent SMILES must contain at most one `*` attachment atom")
    dummy = dummy_atoms[0]
    neighbors = dummy.GetNeighbors()
    if len(neighbors) != 1:
        raise ValueError("Substituent `*` attachment atom must have exactly one neighbor")
    return neighbors[0].GetIdx(), dummy.GetIdx()


def _copy_fragment(
    tgt_mol: Chem.rdchem.RWMol,
    src_mol: Chem.rdchem.Mol,
    omitted_idx: int | None,
) -> dict[int, int]:
    idx_map: dict[int, int] = {}
    for atom in src_mol.GetAtoms():
        if atom.GetIdx() == omitted_idx:
            continue
        new_atom = Chem.Atom(atom)
        new_atom.SetAtomMapNum(0)
        idx_map[atom.GetIdx()] = tgt_mol.AddAtom(new_atom)

    for bond in src_mol.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin == omitted_idx or end == omitted_idx:
            continue
        tgt_mol.AddBond(idx_map[begin], idx_map[end], bond.GetBondType())
    return idx_map


def _attach_fragment(
    mol: Chem.rdchem.RWMol,
    attach_idx: int,
    fragment_smiles: str,
) -> None:
    src_mol = Chem.MolFromSmiles(fragment_smiles)
    if src_mol is None:
        raise ValueError(f"Invalid substituent SMILES: {fragment_smiles}")

    src_attach_idx, src_dummy_idx = _source_attachment(src_mol)
    attach_atom = mol.GetAtomWithIdx(attach_idx)
    replace_dummy = attach_atom.GetSymbol() == "*"
    target_idx = attach_idx

    if replace_dummy:
        neighbors = attach_atom.GetNeighbors()
        if len(neighbors) != 1:
            raise ValueError("Target `*` attachment atom must have exactly one neighbor")
        neighbor = neighbors[0]
        bond = mol.GetBondBetweenAtoms(attach_idx, neighbor.GetIdx())
        if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
            raise ValueError("Target `*` attachment atom must link through a single bond")
        target_idx = neighbor.GetIdx()
        mol.RemoveBond(attach_idx, target_idx)

    idx_map = _copy_fragment(mol, src_mol, src_dummy_idx)
    mol.AddBond(target_idx, idx_map[src_attach_idx], Chem.BondType.SINGLE)

    if replace_dummy:
        mol.RemoveAtom(attach_idx)


def _canonical_smiles(mol: Chem.rdchem.Mol) -> str:
    output = Chem.Mol(mol)
    for atom in output.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.SanitizeMol(output)
    return Chem.MolToSmiles(output, canonical=True, isomericSmiles=True)


def _collect_valid_smiles(states: list[Chem.rdchem.RWMol]) -> list[str]:
    RDLogger.DisableLog("rdApp.*")
    try:
        smiles = set()
        for state in states:
            try:
                smiles.add(_canonical_smiles(state))
            except Exception:
                continue
        return sorted(smiles)
    finally:
        RDLogger.EnableLog("rdApp.*")


def _ring_atom_indices(
    mol: Chem.rdchem.Mol,
    source_ring: Iterable[int],
) -> list[int]:
    indices: list[int] = []
    for source_atom_idx in source_ring:
        atom_idx = _find_atom_by_source_index(mol, source_atom_idx)
        if atom_idx is not None:
            indices.append(atom_idx)
    return indices


def _expand_atom_group(
    states: list[Chem.rdchem.RWMol],
    desc: GroupDesc,
    fragments: list[str],
) -> list[Chem.rdchem.RWMol]:
    counts = _repeat_counts(desc, 1)
    if counts != [1]:
        raise ValueError("Atom-indexed Markush groups only support single substitution")

    next_states: list[Chem.rdchem.RWMol] = []
    for state in states:
        attach_idx = _find_atom_by_source_index(state, int(desc.id))
        if attach_idx is None:
            raise ValueError(f"Atom index {int(desc.id)} is not present in the molecule")
        for fragment in fragments:
            new_state = Chem.RWMol(state)
            _attach_fragment(new_state, attach_idx, fragment)
            next_states.append(new_state)
    return next_states


def _source_ring_for_atom(
    source_rings: tuple[tuple[int, ...], ...],
    source_atom_idx: int,
) -> tuple[int, ...] | None:
    for source_ring in source_rings:
        if source_atom_idx in source_ring:
            return source_ring
    return None


def _expand_atom_group_copies(
    states: list[Chem.rdchem.RWMol],
    desc: GroupDesc,
    fragments: list[str],
    definitions: Mapping[str, DefinitionValue],
    source_rings: tuple[tuple[int, ...], ...],
    max_outputs: int,
) -> list[Chem.rdchem.RWMol] | None:
    count = _resolve_repetition_count(desc, definitions)
    if count is None:
        return None
    if count == 1:
        return _expand_atom_group(states, desc, fragments)

    next_states: list[Chem.rdchem.RWMol] = []
    for state in states:
        attach_idx = _find_atom_by_source_index(state, int(desc.id))
        if attach_idx is None:
            raise ValueError(f"Atom index {int(desc.id)} is not present in the molecule")

        attach_atom = state.GetAtomWithIdx(attach_idx)
        if attach_atom.GetSymbol() == "*":
            neighbors = attach_atom.GetNeighbors()
            if len(neighbors) != 1:
                raise ValueError("Target `*` attachment atom must have exactly one neighbor")
            ring_anchor_idx = neighbors[0].GetIdx()
        else:
            ring_anchor_idx = attach_idx

        anchor_map_num = state.GetAtomWithIdx(ring_anchor_idx).GetAtomMapNum()
        source_ring = _source_ring_for_atom(source_rings, anchor_map_num - 1)
        if source_ring is None:
            raise ValueError(
                f"Cannot copy atom-indexed group `{str(desc)}` without an anchor ring"
            )

        sites = [
            site
            for site in _ring_atom_indices(state, source_ring)
            if site != ring_anchor_idx
        ]
        if count - 1 > len(sites):
            raise ValueError(
                f"Cannot place {count} copies of `{str(desc)}` on the anchor ring"
            )

        for extra_sites in combinations(sites, count - 1):
            for fragment in fragments:
                new_state = Chem.RWMol(state)
                for site in extra_sites:
                    _attach_fragment(new_state, site, fragment)
                _attach_fragment(new_state, attach_idx, fragment)
                next_states.append(new_state)
                if len(next_states) > max_outputs:
                    raise ValueError(f"Markush expansion exceeded max_outputs={max_outputs}")
    return next_states


def _expand_atom_repetition(
    states: list[Chem.rdchem.RWMol],
    desc: GroupDesc,
    definitions: Mapping[str, DefinitionValue],
) -> list[Chem.rdchem.RWMol]:
    count = _resolve_repetition_count(desc, definitions)
    if count is None:
        raise ValueError(f"No repetition count found for group `{str(desc)}`")

    repeated_desc = GroupDesc(
        id=desc.id,
        symbol=desc.symbol,
        script=desc.script,
        prime=desc.prime,
        multiple=str(count),
        is_circle=desc.is_circle,
        is_dummy=desc.is_dummy,
    )
    next_states: list[Chem.rdchem.RWMol] = []
    for state in states:
        atom_idx = _find_atom_by_source_index(state, int(desc.id))
        if atom_idx is None:
            raise ValueError(f"Atom index {int(desc.id)} is not present in the molecule")
        new_state = Chem.RWMol(state)
        is_markush = chem_utils.carbon_chain_repetition_process(
            new_state,
            atom_idx,
            repeated_desc,
            is_markush=False,
        )
        if is_markush:
            raise ValueError(f"Failed to expand carbon-chain repeat `{str(desc)}`")
        next_states.append(new_state)
    return next_states


def _expand_ring_group(
    states: list[Chem.rdchem.RWMol],
    desc: GroupDesc,
    fragments: list[str],
    source_rings: tuple[tuple[int, ...], ...],
    max_outputs: int,
) -> list[Chem.rdchem.RWMol]:
    source_ring = source_rings[int(desc.id)]
    next_states: list[Chem.rdchem.RWMol] = []

    for state in states:
        sites = _ring_atom_indices(state, source_ring)
        if not sites:
            raise ValueError(f"Ring index {int(desc.id)} is not present in the molecule")
        for count in _repeat_counts(desc, len(sites)):
            if count < 1 or count > len(sites):
                continue
            for site_group in combinations(sites, count):
                for fragment in fragments:
                    new_state = Chem.RWMol(state)
                    for site in site_group:
                        _attach_fragment(new_state, site, fragment)
                    next_states.append(new_state)
                    if len(next_states) > max_outputs:
                        raise ValueError(
                            f"Markush expansion exceeded max_outputs={max_outputs}"
                        )
    return next_states


def _format_substitution_outputs(
    smiles: list[str],
    annotation_variants: list[str],
    ext: str,
    force_esmiles: bool,
) -> list[str]:
    if not force_esmiles and not annotation_variants:
        return smiles

    annotations = annotation_variants or [""]
    return [
        Translator.build_esmi(smile, annotation, ext)
        for smile in smiles
        for annotation in annotations
    ]


def _format_preserved_group_outputs(
    states: list[Chem.rdchem.RWMol],
    preserved_groups: str,
    source_rings: tuple[tuple[int, ...], ...],
    annotation_variants: list[str],
    ext: str,
) -> list[str]:
    outputs: set[str] = set()
    annotations = annotation_variants or [""]
    for state in states:
        smiles = _canonical_smiles(state)
        remap_mol = Chem.Mol(state)
        Chem.SanitizeMol(remap_mol)
        remapped_groups = chem_utils.remap_groups(
            remap_mol,
            preserved_groups,
            source_rings,
        )
        for annotation in annotations:
            outputs.add(Translator.build_esmi(smiles, remapped_groups + annotation, ext))
    return sorted(outputs)


def _substitute_precompat_records(
    groups: str,
    definitions: Mapping[str, DefinitionValue],
    *,
    max_outputs: int,
    error_msg: bool,
) -> list[str]:
    variants = [""]
    found_record = False
    for match in _PRECOMPAT_RECORD_PATTERN.finditer(groups):
        found_record = True
        record = match.group(0)
        substruct_match = _SUBSTRUCT_RECORD_PATTERN.match(record)
        if substruct_match:
            body = substruct_match.group("body")
            body_variants = _substitute_markush_outputs(
                body,
                definitions,
                max_outputs=max_outputs,
                error_msg=error_msg,
                force_esmiles=True,
            )
            record_variants = [f"<s>{body_variant}</s>" for body_variant in body_variants]
        else:
            record_variants = [record]

        variants = [
            prefix + record_variant
            for prefix in variants
            for record_variant in record_variants
        ]
        if len(variants) > max_outputs:
            raise ValueError(f"Markush expansion exceeded max_outputs={max_outputs}")
    return variants if found_record else []


def _substitute_markush_outputs(
    caption: str,
    definitions: Mapping[str, DefinitionValue],
    *,
    max_outputs: int,
    error_msg: bool,
    force_esmiles: bool = False,
) -> list[str]:
    raw_caption = str(caption).strip()
    parsed = Translator.parse_caption(raw_caption, return_mol=True, error_msg=error_msg)
    if parsed is None:
        raise ValueError(f"Invalid E-SMILES caption: {caption}")

    mol, groups, ext = parsed
    annotation_variants = _substitute_precompat_records(
        groups,
        definitions,
        max_outputs=max_outputs,
        error_msg=error_msg,
    )
    groups = Translator.repair_atom_group_indices(mol, groups, error_msg=error_msg)
    source_rings = mol.GetRingInfo().AtomRings()
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)

    states = [Chem.RWMol(mol)]
    special_id_groups: list[str] = []
    for desc in Translator.parse_groups(groups):
        if desc.is_dummy or desc.is_circle:
            continue
        if _is_special_id_label(desc):
            special_id_groups.append(f"{Tokens.atom_start}{int(desc.id)}:{str(desc)}{Tokens.atom_end}")
            continue
        if isinstance(desc.id, RingIndex) and desc.id.virtual:
            continue
        if isinstance(desc.id, AtomIndex) and _is_carbon_chain_repeat(desc):
            states = _expand_atom_repetition(states, desc, definitions)
            continue

        fragments = _resolve_group_smiles(desc, definitions)
        if not fragments:
            raise ValueError(f"No Markush definition found for group `{str(desc)}`")

        if isinstance(desc.id, AtomIndex):
            copied_states = None
            if desc.multiple:
                copied_states = _expand_atom_group_copies(
                    states,
                    desc,
                    fragments,
                    definitions,
                    source_rings,
                    max_outputs,
                )
            states = copied_states if copied_states is not None else _expand_atom_group(
                states,
                desc,
                fragments,
            )
        elif isinstance(desc.id, RingIndex):
            if int(desc.id) >= len(source_rings):
                raise ValueError(f"Ring index {int(desc.id)} is out of range")
            states = _expand_ring_group(
                states,
                desc,
                fragments,
                source_rings,
                max_outputs=max_outputs,
            )

    smiles = _collect_valid_smiles(states)
    if not smiles:
        raise ValueError("No valid SMILES generated from Markush substitution")
    if special_id_groups:
        return _format_preserved_group_outputs(
            states,
            "".join(special_id_groups),
            source_rings,
            annotation_variants,
            ext,
        )
    return _format_substitution_outputs(smiles, annotation_variants, ext, force_esmiles)


def substitute_markush(
    caption: str,
    definitions: Mapping[str, DefinitionValue],
    *,
    max_outputs: int = 1024,
    error_msg: bool = False,
) -> str | list[str]:
    """Substitute Markush definitions into an E-SMILES caption.

    ``definitions`` maps labels such as ``R1`` or ``R[1]`` to either an
    abbreviation (``Me``) or a SMILES fragment. If the fragment contains ``*``,
    that atom is treated as the attachment point and removed during merging.
    Ring-indexed groups are expanded over possible ring attachment positions.
    """
    outputs = _substitute_markush_outputs(
        str(caption).strip(),
        definitions,
        max_outputs=max_outputs,
        error_msg=error_msg,
    )
    if len(outputs) == 1:
        return outputs[0]
    return outputs


__all__ = ["substitute_markush"]
