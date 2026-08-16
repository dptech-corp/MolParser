"""Markush substituent expansion for E-SMILES captions."""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from itertools import combinations
from typing import Iterable, Literal

from rdkit import Chem, RDLogger

try:
    from . import chem_utils
    from .translator import AtomIndex, GroupDesc, RingIndex, Tokens, Translator
except ImportError:  # Support running from package directory as working directory.
    import chem_utils
    from translator import AtomIndex, GroupDesc, RingIndex, Tokens, Translator


DefinitionValue = int | str | Sequence[str]
RepeatPolicy = Literal["preserve", "best_effort", "strict"]
TerminalPolicy = Literal["preserve", "hydrogen"]


class _UnexpandableRepeat(ValueError):
    """A valid repeat annotation whose physical topology is under-specified."""

_PRECOMPAT_RECORD_PATTERN = re.compile(
    r"<s>.*?</s>|<g>.*?</g>|<r><v>\d+:.+?</r>|<v>.*?</v>",
    re.DOTALL,
)
_SUBSTRUCT_RECORD_PATTERN = re.compile(r"^<s>(?P<body>.*)</s>$", re.DOTALL)
_SGROUP_COUNT_PATTERN = re.compile(r"\|Sg:(?P<count>[^|]+)\|")
_SGROUP_SYMBOL_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
_PRESERVED_ATOM_RECORD_PATTERN = re.compile(
    r"<(?P<tag>a|d|c)>(?P<index>\d+):(?P<value>.*?)</(?P=tag)>",
    re.DOTALL,
)
_LOCAL_SGROUP_RECORD_PATTERN = re.compile(r"<g>(?P<body>.*?)</g>", re.DOTALL)
_LOCAL_SGROUP_PORT_PATTERN = re.compile(r"\[(?P<inner>\d+):(?P<outer>\d+)\]")


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


def _resolve_sgroup_count(
    count: str,
    definitions: Mapping[str, DefinitionValue],
) -> str:
    """Resolve a symbolic ``|Sg:...|`` count without expanding its graph."""
    clean_count = count.strip()
    if not _SGROUP_SYMBOL_PATTERN.fullmatch(clean_count):
        return count

    definition_lookup = _definition_lookup(definitions)
    if clean_count not in definition_lookup:
        return count
    value = definition_lookup[clean_count]

    if isinstance(value, Sequence) and not isinstance(value, str):
        raise ValueError(
            f"S-group count `{clean_count}` must resolve to one positive integer"
        )
    if isinstance(value, bool):
        raise ValueError(
            f"S-group count `{clean_count}` must resolve to one positive integer"
        )
    if isinstance(value, int):
        resolved = value
    elif isinstance(value, str) and value.strip().isdigit():
        resolved = int(value.strip())
    else:
        raise ValueError(
            f"S-group count `{clean_count}` must resolve to one positive integer"
        )
    if resolved < 1:
        raise ValueError(
            f"S-group count `{clean_count}` must resolve to one positive integer"
        )
    return str(resolved)


def _substitute_sgroup_counts(
    text: str,
    definitions: Mapping[str, DefinitionValue],
) -> str:
    def replace(match: re.Match) -> str:
        resolved = _resolve_sgroup_count(match.group("count"), definitions)
        return f"|Sg:{resolved}|"

    return _SGROUP_COUNT_PATTERN.sub(replace, text)


def _concrete_sgroup_count(
    count: str,
    definitions: Mapping[str, DefinitionValue],
) -> int | None:
    """Return one concrete positive count, or ``None`` for a residual symbol/range."""
    clean_count = count.strip()
    if clean_count.isdigit():
        value = int(clean_count)
        if value < 1:
            raise ValueError("S-group repeat count must be one positive integer")
        return value
    numeric_range = re.fullmatch(r"(\d+)-(\d+)", clean_count)
    if numeric_range is not None:
        start, end = (int(value) for value in numeric_range.groups())
        if start < 1 or end < start:
            raise ValueError("S-group repeat range must be positive and increasing")
        return None
    if not _SGROUP_SYMBOL_PATTERN.fullmatch(clean_count):
        return None
    if clean_count not in _definition_lookup(definitions):
        return None
    return int(_resolve_sgroup_count(clean_count, definitions))


def _parse_local_sgroup_record(
    record: str,
    definitions: Mapping[str, DefinitionValue],
) -> tuple[tuple[tuple[int, int], tuple[int, int]], int] | None:
    matched = _LOCAL_SGROUP_RECORD_PATTERN.fullmatch(record)
    if matched is None:
        return None
    body = matched.group("body")
    count_matches = list(_SGROUP_COUNT_PATTERN.finditer(body))
    ports = [
        (int(port.group("inner")), int(port.group("outer")))
        for port in _LOCAL_SGROUP_PORT_PATTERN.finditer(body)
    ]
    if len(count_matches) != 1 or len(ports) != 2:
        return None
    residue = _LOCAL_SGROUP_PORT_PATTERN.sub("", body)
    residue = _SGROUP_COUNT_PATTERN.sub("", residue)
    if residue.strip(" \t\r\n:"):
        return None
    count = _concrete_sgroup_count(count_matches[0].group("count"), definitions)
    if count is None:
        return None
    return (ports[0], ports[1]), count


def _top_level_dummy_records(groups: str) -> list[tuple[int, int, int]]:
    protected = [
        (match.start(), match.end())
        for match in _PRECOMPAT_RECORD_PATTERN.finditer(groups)
    ]
    records: list[tuple[int, int, int]] = []
    for match in _PRESERVED_ATOM_RECORD_PATTERN.finditer(groups):
        if any(start <= match.start() and match.end() <= end for start, end in protected):
            continue
        if match.group("value") == Tokens.dummy:
            records.append((match.start(), match.end(), int(match.group("index"))))
    return records


def _remove_spans(text: str, spans: Iterable[tuple[int, int]]) -> str:
    output = text
    for start, end in sorted(spans, reverse=True):
        output = output[:start] + output[end:]
    return output


def _preserved_top_level_atom_records(groups: str) -> str:
    """Keep dummy and special-id records outside nested pre-compatible records."""
    top_level_groups = _PRECOMPAT_RECORD_PATTERN.sub("", groups)
    records: list[str] = []
    for match in _PRESERVED_ATOM_RECORD_PATTERN.finditer(top_level_groups):
        value = match.group("value")
        if (
            match.group("tag") == "c"
            or value == Tokens.dummy
            or value.startswith(f"{Tokens.special_id}[")
        ):
            records.append(match.group(0))
    return "".join(records)


def _group_record(desc: GroupDesc) -> str:
    """Serialize one parsed top-level group for best-effort preservation."""
    if isinstance(desc.id, AtomIndex):
        tag = "d" if desc.is_dummy else "a"
    elif isinstance(desc.id, RingIndex) and not desc.id.virtual:
        tag = "r"
    else:
        return ""
    return f"<{tag}>{int(desc.id)}:{str(desc)}</{tag}>"


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
        count = int(desc.multiple)
        if count < 1:
            raise ValueError(f"Multiplicity ?{desc.multiple} must be positive")
        return count

    value = _definition_lookup(definitions).get(desc.multiple)
    if isinstance(value, Sequence) and not isinstance(value, str):
        raise ValueError(f"Multiplicity ?{desc.multiple} must resolve to one integer")
    if value is None:
        return None
    if isinstance(value, bool):
        raise ValueError(f"Multiplicity ?{desc.multiple} must resolve to one integer")
    try:
        count = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Multiplicity ?{desc.multiple} must resolve to one integer") from exc
    if count < 1:
        raise ValueError(f"Multiplicity ?{desc.multiple} must be positive")
    return count


def _resolve_repetition_counts(
    desc: GroupDesc,
    definitions: Mapping[str, DefinitionValue],
) -> list[int] | None:
    if desc.multiple and re.fullmatch(r"\d+-\d+", desc.multiple):
        start_text, end_text = desc.multiple.split("-", 1)
        start, end = int(start_text), int(end_text)
        if start < 1 or end < start:
            raise ValueError(f"Unsupported multiplicity: ?{desc.multiple}")
        return list(range(start, end + 1))
    count = _resolve_repetition_count(desc, definitions)
    return None if count is None else [count]


def _is_carbon_chain_repeat(desc: GroupDesc) -> bool:
    return (desc.symbol == "CH2") or (desc.symbol == "CH" and desc.script == "2")


def _is_safe_carbon_chain_repeat_target(atom: Chem.rdchem.Atom) -> bool:
    """Reject targets whose chemistry would be destroyed by chain insertion."""
    if atom.GetSymbol() not in {"*", "C"}:
        return False
    if atom.GetDegree() not in {1, 2} or atom.GetIsAromatic() or atom.IsInRing():
        return False
    if (
        atom.GetIsotope() != 0
        or atom.GetFormalCharge() != 0
        or atom.GetNumRadicalElectrons() != 0
        or atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
    ):
        return False
    return all(
        bond.GetBondType() == Chem.BondType.SINGLE and not bond.GetIsAromatic()
        for bond in atom.GetBonds()
    )


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
    omitted_replacement_idx: int | None = None,
) -> tuple[
    dict[int, int],
    list[tuple[int, int, tuple[int, int], Chem.rdchem.BondStereo]],
]:
    idx_map: dict[int, int] = {}
    for atom in src_mol.GetAtoms():
        if atom.GetIdx() == omitted_idx:
            continue
        new_atom = Chem.Atom(atom)
        new_atom.SetAtomMapNum(0)
        idx_map[atom.GetIdx()] = tgt_mol.AddAtom(new_atom)

    copied_bonds: list[tuple[Chem.rdchem.Bond, Chem.rdchem.Bond]] = []
    for bond in src_mol.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin == omitted_idx or end == omitted_idx:
            continue
        tgt_mol.AddBond(idx_map[begin], idx_map[end], bond.GetBondType())
        copied = tgt_mol.GetBondBetweenAtoms(idx_map[begin], idx_map[end])
        copied.SetBondDir(bond.GetBondDir())
        copied.SetIsAromatic(bond.GetIsAromatic())
        copied.SetIsConjugated(bond.GetIsConjugated())
        copied_bonds.append((bond, copied))

    pending_stereo: list[
        tuple[int, int, tuple[int, int], Chem.rdchem.BondStereo]
    ] = []
    for source_bond, copied_bond in copied_bonds:
        stereo_atoms = tuple(source_bond.GetStereoAtoms())
        if stereo_atoms:
            mapped_stereo: list[int] = []
            for atom_idx in stereo_atoms:
                if atom_idx == omitted_idx:
                    if omitted_replacement_idx is None:
                        raise _UnexpandableRepeat(
                            "Cannot preserve double-bond stereo across the attachment"
                        )
                    mapped_stereo.append(omitted_replacement_idx)
                elif atom_idx in idx_map:
                    mapped_stereo.append(idx_map[atom_idx])
                else:
                    raise _UnexpandableRepeat(
                        "Cannot remap double-bond stereo across the attachment"
                    )
            if len(mapped_stereo) != 2:
                raise _UnexpandableRepeat(
                    "Double-bond stereo requires two mapped reference atoms"
                )
            pending_stereo.append(
                (
                    copied_bond.GetBeginAtomIdx(),
                    copied_bond.GetEndAtomIdx(),
                    (mapped_stereo[0], mapped_stereo[1]),
                    source_bond.GetStereo(),
                )
            )
        elif source_bond.GetStereo() != Chem.BondStereo.STEREONONE:
            copied_bond.SetStereo(source_bond.GetStereo())

    if omitted_idx is not None:
        for source_atom in src_mol.GetAtoms():
            neighbors = [neighbor.GetIdx() for neighbor in source_atom.GetNeighbors()]
            if omitted_idx not in neighbors:
                continue
            copied_atom = tgt_mol.GetAtomWithIdx(idx_map[source_atom.GetIdx()])
            if copied_atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
                continue
            omitted_position = neighbors.index(omitted_idx)
            if (len(neighbors) - 1 - omitted_position) % 2:
                copied_atom.InvertChirality()
    return idx_map, pending_stereo


def _restore_stereo_snapshots(
    mol: Chem.rdchem.RWMol,
    snapshots: list[tuple[int, int, tuple[int, int], Chem.rdchem.BondStereo]],
    removed_idx: int | None = None,
) -> None:
    """Restore stereo metadata after attachment bonds/atoms have been replaced."""
    def adjusted(atom_idx: int) -> int:
        if removed_idx is not None and atom_idx > removed_idx:
            return atom_idx - 1
        return atom_idx

    for begin_idx, end_idx, stereo_atoms, stereo in snapshots:
        bond = mol.GetBondBetweenAtoms(adjusted(begin_idx), adjusted(end_idx))
        if bond is None:
            raise _UnexpandableRepeat(
                "Cannot locate double bond after attachment replacement"
            )
        mapped = (adjusted(stereo_atoms[0]), adjusted(stereo_atoms[1]))
        if mapped[0] == mapped[1]:
            raise _UnexpandableRepeat(
                "Cannot preserve double-bond stereo across the attachment"
            )
        bond.SetStereoAtoms(mapped[0], mapped[1])
        bond.SetStereo(stereo)


def _target_stereo_snapshots(
    mol: Chem.rdchem.RWMol,
    attachment_idx: int,
) -> list[tuple[int, int, tuple[int, int], Chem.rdchem.BondStereo]]:
    """Capture alkene stereo before deleting its directional attachment bond."""
    snapshots: list[tuple[int, int, tuple[int, int], Chem.rdchem.BondStereo]] = []
    for bond in mol.GetBonds():
        stereo_atoms = tuple(bond.GetStereoAtoms())
        if attachment_idx not in stereo_atoms:
            continue
        if len(stereo_atoms) != 2:
            raise _UnexpandableRepeat(
                "Double-bond stereo requires two reference atoms"
            )
        snapshots.append(
            (
                bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                (stereo_atoms[0], stereo_atoms[1]),
                bond.GetStereo(),
            )
        )
    return snapshots


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
    invert_target_chirality = False
    target_stereo = []

    if replace_dummy:
        neighbors = attach_atom.GetNeighbors()
        if len(neighbors) != 1:
            raise ValueError("Target `*` attachment atom must have exactly one neighbor")
        neighbor = neighbors[0]
        bond = mol.GetBondBetweenAtoms(attach_idx, neighbor.GetIdx())
        if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
            raise ValueError("Target `*` attachment atom must link through a single bond")
        target_stereo = _target_stereo_snapshots(mol, attach_idx)
        target_idx = neighbor.GetIdx()
        target_neighbors = [item.GetIdx() for item in neighbor.GetNeighbors()]
        if neighbor.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            dummy_position = target_neighbors.index(attach_idx)
            invert_target_chirality = (
                (len(target_neighbors) - 1 - dummy_position) % 2 == 1
            )
        mol.RemoveBond(attach_idx, target_idx)

    idx_map, source_stereo = _copy_fragment(
        mol,
        src_mol,
        src_dummy_idx,
        omitted_replacement_idx=target_idx if src_dummy_idx is not None else None,
    )
    mol.AddBond(target_idx, idx_map[src_attach_idx], Chem.BondType.SINGLE)

    if replace_dummy:
        mapped_target_stereo = [
            (
                begin_idx,
                end_idx,
                tuple(
                    idx_map[src_attach_idx] if atom_idx == attach_idx else atom_idx
                    for atom_idx in stereo_atoms
                ),
                stereo,
            )
            for begin_idx, end_idx, stereo_atoms, stereo in target_stereo
        ]
        if invert_target_chirality:
            mol.GetAtomWithIdx(target_idx).InvertChirality()
        mol.RemoveAtom(attach_idx)
        _restore_stereo_snapshots(
            mol,
            mapped_target_stereo + source_stereo,
            removed_idx=attach_idx,
        )
    else:
        _restore_stereo_snapshots(mol, source_stereo)


def _canonical_smiles(mol: Chem.rdchem.Mol) -> str:
    output = Chem.Mol(mol)
    for atom in output.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.SanitizeMol(output)
    Chem.SetDoubleBondNeighborDirections(output)
    return Chem.MolToSmiles(output, canonical=True, isomericSmiles=True)


def _copy_atom_subset(
    target: Chem.rdchem.RWMol,
    source: Chem.rdchem.Mol,
    atom_indices: set[int],
) -> dict[int, int] | None:
    index_map: dict[int, int] = {}
    for atom_idx in sorted(atom_indices):
        atom = Chem.Atom(source.GetAtomWithIdx(atom_idx))
        atom.SetAtomMapNum(0)
        index_map[atom_idx] = target.AddAtom(atom)
    copied_bonds: list[tuple[Chem.rdchem.Bond, Chem.rdchem.Bond]] = []
    for bond in source.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin not in atom_indices or end not in atom_indices:
            continue
        target.AddBond(index_map[begin], index_map[end], bond.GetBondType())
        copied = target.GetBondBetweenAtoms(index_map[begin], index_map[end])
        copied.SetBondDir(bond.GetBondDir())
        copied.SetIsAromatic(bond.GetIsAromatic())
        copied.SetIsConjugated(bond.GetIsConjugated())
        copied_bonds.append((bond, copied))
    for source_bond, copied_bond in copied_bonds:
        stereo_atoms = tuple(source_bond.GetStereoAtoms())
        if stereo_atoms:
            if len(stereo_atoms) != 2 or any(
                atom_idx not in index_map for atom_idx in stereo_atoms
            ):
                return None
            copied_bond.SetStereoAtoms(
                index_map[stereo_atoms[0]],
                index_map[stereo_atoms[1]],
            )
        copied_bond.SetStereo(source_bond.GetStereo())
    return index_map


def _component_without_edges(
    mol: Chem.rdchem.Mol,
    start: int,
    excluded_edges: set[frozenset[int]],
) -> set[int]:
    visited: set[int] = set()
    pending = [start]
    while pending:
        atom_idx = pending.pop()
        if atom_idx in visited:
            continue
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if frozenset((atom_idx, neighbor_idx)) in excluded_edges:
                continue
            if neighbor_idx not in visited:
                pending.append(neighbor_idx)
    return visited


def _terminal_dummy_indices(
    mol: Chem.rdchem.Mol,
    source_indices: Sequence[int],
) -> list[int] | None:
    if len(source_indices) != 2:
        return None
    resolved: list[int] = []
    for source_idx in source_indices:
        atom_idx = _find_atom_by_source_index(mol, source_idx)
        if atom_idx is None:
            return None
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() != "*" or atom.GetDegree() != 1:
            return None
        bond = atom.GetBonds()[0]
        if bond.GetBondType() != Chem.BondType.SINGLE:
            return None
        resolved.append(atom_idx)
    if resolved[0] == resolved[1]:
        return None
    return resolved


def _convert_terminal_dummies_to_hydrogen(
    mol: Chem.rdchem.RWMol,
    atom_indices: Iterable[int],
) -> bool:
    """Turn terminal ``*`` atoms into explicit H before RDKit removes them."""
    for atom_idx in atom_indices:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() != "*" or atom.GetDegree() != 1:
            return False
        atom.SetAtomicNum(1)
        atom.SetIsotope(0)
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(True)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        atom.SetAtomMapNum(0)
    return True


def _expand_whole_sru_state(
    state: Chem.rdchem.RWMol,
    endpoint_source_indices: Sequence[int],
    count: int,
    terminal_policy: TerminalPolicy,
) -> Chem.rdchem.RWMol | None:
    source = Chem.Mol(state)
    endpoints = _terminal_dummy_indices(source, endpoint_source_indices)
    if endpoints is None:
        return None
    left_dummy, right_dummy = endpoints
    left_neighbor = source.GetAtomWithIdx(left_dummy).GetNeighbors()[0].GetIdx()
    right_neighbor = source.GetAtomWithIdx(right_dummy).GetNeighbors()[0].GetIdx()
    left_bond = source.GetBondBetweenAtoms(left_dummy, left_neighbor)
    right_bond = source.GetBondBetweenAtoms(right_dummy, right_neighbor)
    if left_bond.GetBondType() != right_bond.GetBondType():
        return None
    if left_bond.GetBondType() != Chem.BondType.SINGLE:
        return None
    endpoint_set = {left_dummy, right_dummy}
    endpoint_dependent_stereo = any(
        bond.GetStereo() != Chem.BondStereo.STEREONONE
        and endpoint_set.intersection(bond.GetStereoAtoms())
        for bond in source.GetBonds()
    )
    if endpoint_dependent_stereo and (
        count > 1 or terminal_policy == "hydrogen"
    ):
        # Removing/converting a dummy used as a double-bond stereo reference
        # would silently erase E/Z information.  Keep the E-SMILES annotation
        # until an explicit stereo-reference stitching rule is available.
        return None
    core_atoms = set(range(source.GetNumAtoms())) - {left_dummy, right_dummy}
    if not core_atoms:
        return None
    terminal_edges = {
        frozenset((left_dummy, left_neighbor)),
        frozenset((right_dummy, right_neighbor)),
    }
    if (
        _component_without_edges(source, next(iter(core_atoms)), terminal_edges)
        != core_atoms
    ):
        return None

    combined = Chem.Mol(source)
    for _ in range(count - 1):
        combined = Chem.CombineMols(combined, source)
    editable = Chem.RWMol(combined)
    atom_count = source.GetNumAtoms()
    for copy_idx in range(count - 1):
        editable.AddBond(
            copy_idx * atom_count + right_neighbor,
            (copy_idx + 1) * atom_count + left_neighbor,
            left_bond.GetBondType(),
        )

    to_remove = {
        copy_idx * atom_count + right_dummy for copy_idx in range(count - 1)
    }
    to_remove.update(
        copy_idx * atom_count + left_dummy for copy_idx in range(1, count)
    )
    if terminal_policy == "hydrogen":
        if not _convert_terminal_dummies_to_hydrogen(
            editable,
            (left_dummy, (count - 1) * atom_count + right_dummy),
        ):
            return None
    for atom_idx in sorted(to_remove, reverse=True):
        editable.RemoveAtom(atom_idx)
    try:
        Chem.SanitizeMol(editable)
        if terminal_policy == "hydrogen":
            editable = Chem.RWMol(Chem.RemoveHs(Chem.Mol(editable)))
            Chem.SanitizeMol(editable)
    except Exception:
        return None
    return editable


def _expand_local_sgroup_state(
    state: Chem.rdchem.RWMol,
    ports: tuple[tuple[int, int], tuple[int, int]],
    count: int,
    terminal_source_indices: Sequence[int],
    terminal_policy: TerminalPolicy,
) -> Chem.rdchem.RWMol | None:
    source = Chem.Mol(state)
    resolved_ports: list[tuple[int, int]] = []
    for inner_source, outer_source in ports:
        inner = _find_atom_by_source_index(source, inner_source)
        outer = _find_atom_by_source_index(source, outer_source)
        if inner is None or outer is None or inner == outer:
            return None
        resolved_ports.append((inner, outer))
    (left_inner, left_outer), (right_inner, right_outer) = resolved_ports
    if len({left_inner, left_outer, right_inner, right_outer}) != 4:
        return None
    left_bond = source.GetBondBetweenAtoms(left_inner, left_outer)
    right_bond = source.GetBondBetweenAtoms(right_inner, right_outer)
    if left_bond is None or right_bond is None:
        return None
    if left_bond.GetBondType() != right_bond.GetBondType():
        return None
    if left_bond.GetBondType() != Chem.BondType.SINGLE:
        return None
    cut_edges = {
        frozenset((left_inner, left_outer)),
        frozenset((right_inner, right_outer)),
    }
    repeat_atoms = _component_without_edges(source, left_inner, cut_edges)
    if right_inner not in repeat_atoms:
        return None
    if left_outer in repeat_atoms or right_outer in repeat_atoms:
        return None

    crossing_edges = {
        frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        for bond in source.GetBonds()
        if (bond.GetBeginAtomIdx() in repeat_atoms)
        != (bond.GetEndAtomIdx() in repeat_atoms)
    }
    if crossing_edges != cut_edges:
        return None

    terminal_indices = _terminal_dummy_indices(source, terminal_source_indices)
    if terminal_policy == "hydrogen" and terminal_indices is None:
        return None
    if terminal_indices is not None and any(
        atom_idx in repeat_atoms for atom_idx in terminal_indices
    ):
        return None

    editable = Chem.RWMol(source)
    editable.RemoveBond(right_inner, right_outer)
    current_right = right_inner
    for _ in range(count - 1):
        copied = _copy_atom_subset(editable, source, repeat_atoms)
        if copied is None:
            return None
        editable.AddBond(current_right, copied[left_inner], left_bond.GetBondType())
        current_right = copied[right_inner]
    editable.AddBond(current_right, right_outer, right_bond.GetBondType())

    if terminal_policy == "hydrogen" and terminal_indices is not None:
        if not _convert_terminal_dummies_to_hydrogen(editable, terminal_indices):
            return None
    try:
        Chem.SanitizeMol(editable)
        if terminal_policy == "hydrogen":
            editable = Chem.RWMol(Chem.RemoveHs(Chem.Mol(editable)))
            Chem.SanitizeMol(editable)
    except Exception:
        return None
    return editable


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


def _top_level_local_sgroup_matches(groups: str) -> list[re.Match[str]]:
    substruct_spans = [
        (match.start(), match.end())
        for match in re.finditer(r"<s>.*?</s>", groups, re.DOTALL)
    ]
    return [
        match
        for match in _LOCAL_SGROUP_RECORD_PATTERN.finditer(groups)
        if not any(
            start <= match.start() and match.end() <= end
            for start, end in substruct_spans
        )
    ]


def _repeat_fallback(
    policy: RepeatPolicy,
    states: list[Chem.rdchem.RWMol],
    groups: str,
    ext: str,
) -> tuple[list[Chem.rdchem.RWMol], str, str]:
    if policy == "strict":
        raise ValueError("Repeat annotation cannot be expanded unambiguously")
    return states, groups, ext


def _has_index_sensitive_residual(groups: str) -> bool:
    return bool(
        _PRECOMPAT_RECORD_PATTERN.search(groups)
        or _preserved_top_level_atom_records(groups)
    )


def _apply_physical_repeats(
    states: list[Chem.rdchem.RWMol],
    groups: str,
    ext: str,
    definitions: Mapping[str, DefinitionValue],
    *,
    repeat_policy: RepeatPolicy,
    terminal_policy: TerminalPolicy,
    max_outputs: int,
    has_unresolved_groups: bool = False,
) -> tuple[list[Chem.rdchem.RWMol], str, str]:
    if repeat_policy == "preserve":
        return states, groups, ext

    local_matches = _top_level_local_sgroup_matches(groups)
    # A residual indexed group may point inside the repeat unit.  Without a
    # scope rule saying whether it applies once or to every copy, retain the
    # repeat annotation rather than emitting a chemically misleading graph.
    if has_unresolved_groups and (local_matches or ext):
        return _repeat_fallback(repeat_policy, states, groups, ext)
    if ext and local_matches:
        return _repeat_fallback(repeat_policy, states, groups, ext)
    if len(local_matches) > 1:
        return _repeat_fallback(repeat_policy, states, groups, ext)

    if local_matches:
        local_match = local_matches[0]
        parsed = _parse_local_sgroup_record(local_match.group(0), definitions)
        if parsed is None:
            return _repeat_fallback(repeat_policy, states, groups, ext)
        ports, count = parsed
        if count > max_outputs:
            raise ValueError(f"Markush expansion exceeded max_outputs={max_outputs}")
        dummy_records = _top_level_dummy_records(groups)
        dummy_sources = [record[2] for record in dummy_records]
        valid_terminal_records = len(dummy_records) == 2 and all(
            _terminal_dummy_indices(Chem.Mol(state), dummy_sources) is not None
            for state in states
        )
        expanded: list[Chem.rdchem.RWMol] = []
        for state in states:
            output = _expand_local_sgroup_state(
                state,
                ports,
                count,
                dummy_sources if valid_terminal_records else (),
                terminal_policy,
            )
            if output is None:
                return _repeat_fallback(repeat_policy, states, groups, ext)
            expanded.append(output)
        consumed_spans: list[tuple[int, int]] = [
            (local_match.start(), local_match.end())
        ]
        if valid_terminal_records:
            consumed_spans.extend((start, end) for start, end, _ in dummy_records)
        residual_groups = _remove_spans(groups, consumed_spans)
        if _has_index_sensitive_residual(residual_groups):
            return _repeat_fallback(repeat_policy, states, groups, ext)
        return expanded, residual_groups, ext

    if ext:
        matched = _SGROUP_COUNT_PATTERN.fullmatch(ext)
        if matched is None:
            return _repeat_fallback(repeat_policy, states, groups, ext)
        count = _concrete_sgroup_count(matched.group("count"), definitions)
        if count is None:
            return _repeat_fallback(repeat_policy, states, groups, ext)
        if count > max_outputs:
            raise ValueError(f"Markush expansion exceeded max_outputs={max_outputs}")
        dummy_records = _top_level_dummy_records(groups)
        if len(dummy_records) != 2:
            return _repeat_fallback(repeat_policy, states, groups, ext)
        dummy_sources = [record[2] for record in dummy_records]
        expanded = []
        for state in states:
            output = _expand_whole_sru_state(
                state,
                dummy_sources,
                count,
                terminal_policy,
            )
            if output is None:
                return _repeat_fallback(repeat_policy, states, groups, ext)
            expanded.append(output)
        residual_groups = _remove_spans(
            groups,
            [(start, end) for start, end, _ in dummy_records],
        )
        if _has_index_sensitive_residual(residual_groups):
            return _repeat_fallback(repeat_policy, states, groups, ext)
        return expanded, residual_groups, ""

    return states, groups, ext


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


def _available_ring_atom_indices(
    mol: Chem.rdchem.Mol,
    source_ring: Iterable[int],
) -> list[int]:
    """Return ring atoms with one implicit hydrogen available for substitution."""
    sanitized = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(sanitized)
    except Exception:
        return []
    return [
        atom_idx
        for atom_idx in _ring_atom_indices(sanitized, source_ring)
        if sanitized.GetAtomWithIdx(atom_idx).GetNumImplicitHs() > 0
    ]


def _expand_atom_group(
    states: list[Chem.rdchem.RWMol],
    desc: GroupDesc,
    fragments: list[str],
    max_outputs: int,
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
            if len(next_states) > max_outputs:
                raise ValueError(
                    f"Markush expansion exceeded max_outputs={max_outputs}"
                )
    return next_states


def _source_rings_for_atom(
    source_rings: tuple[tuple[int, ...], ...],
    source_atom_idx: int,
) -> list[tuple[int, ...]]:
    return [
        source_ring
        for source_ring in source_rings
        if source_atom_idx in source_ring
    ]


def _expand_atom_group_copies(
    states: list[Chem.rdchem.RWMol],
    desc: GroupDesc,
    fragments: list[str],
    definitions: Mapping[str, DefinitionValue],
    source_rings: tuple[tuple[int, ...], ...],
    max_outputs: int,
) -> list[Chem.rdchem.RWMol] | None:
    counts = _resolve_repetition_counts(desc, definitions)
    if counts is None:
        return None
    if max(counts) > max_outputs or len(counts) > max_outputs:
        raise ValueError(f"Markush expansion exceeded max_outputs={max_outputs}")

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
        matching_rings = _source_rings_for_atom(source_rings, anchor_map_num - 1)
        if len(matching_rings) != 1:
            raise _UnexpandableRepeat(
                f"Cannot copy atom-indexed group `{str(desc)}` without one unique anchor ring"
            )
        source_ring = matching_rings[0]
        if (
            attach_atom.GetSymbol() != "*"
            and ring_anchor_idx
            not in _available_ring_atom_indices(state, source_ring)
        ):
            raise _UnexpandableRepeat(
                f"Atom-indexed group `{str(desc)}` has no available ring hydrogen"
            )

        sites = [
            site
            for site in _available_ring_atom_indices(state, source_ring)
            if site != ring_anchor_idx
        ]
        for count in counts:
            if count - 1 > len(sites):
                raise _UnexpandableRepeat(
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
                        raise ValueError(
                            f"Markush expansion exceeded max_outputs={max_outputs}"
                        )
    return next_states


def _expand_atom_repetition(
    states: list[Chem.rdchem.RWMol],
    desc: GroupDesc,
    definitions: Mapping[str, DefinitionValue],
    max_outputs: int,
) -> list[Chem.rdchem.RWMol]:
    counts = _resolve_repetition_counts(desc, definitions) or []
    if not counts:
        raise ValueError(f"No repetition count found for group `{str(desc)}`")
    if max(counts) > max_outputs or len(states) * len(counts) > max_outputs:
        raise ValueError(f"Markush expansion exceeded max_outputs={max_outputs}")
    next_states: list[Chem.rdchem.RWMol] = []
    for state in states:
        atom_idx = _find_atom_by_source_index(state, int(desc.id))
        if atom_idx is None:
            raise ValueError(f"Atom index {int(desc.id)} is not present in the molecule")
        if not _is_safe_carbon_chain_repeat_target(state.GetAtomWithIdx(atom_idx)):
            raise _UnexpandableRepeat(
                f"Carbon-chain repeat target `{str(desc)}` is not a safe aliphatic site"
            )
        for count in counts:
            repeated_desc = GroupDesc(
                id=desc.id,
                symbol=desc.symbol,
                script=desc.script,
                prime=desc.prime,
                multiple=str(count),
                is_circle=desc.is_circle,
                is_dummy=desc.is_dummy,
            )
            new_state = Chem.RWMol(state)
            is_markush = chem_utils.carbon_chain_repetition_process(
                new_state,
                atom_idx,
                repeated_desc,
                is_markush=False,
            )
            if is_markush:
                raise _UnexpandableRepeat(
                    f"Failed to expand carbon-chain repeat `{str(desc)}`"
                )
            next_states.append(new_state)
            if len(next_states) > max_outputs:
                raise ValueError(f"Markush expansion exceeded max_outputs={max_outputs}")
    return next_states


def _expand_ring_group(
    states: list[Chem.rdchem.RWMol],
    desc: GroupDesc,
    fragments: list[str],
    source_rings: tuple[tuple[int, ...], ...],
    definitions: Mapping[str, DefinitionValue],
    max_outputs: int,
) -> list[Chem.rdchem.RWMol] | None:
    source_ring = source_rings[int(desc.id)]
    next_states: list[Chem.rdchem.RWMol] = []
    counts = (
        [1]
        if not desc.multiple
        else _resolve_repetition_counts(desc, definitions)
    )
    if counts is None:
        return None

    for state in states:
        sites = _available_ring_atom_indices(state, source_ring)
        if not sites:
            continue
        for count in counts:
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
    return next_states or None


def _format_substitution_outputs(
    smiles: list[str],
    annotation_variants: list[str],
    ext: str,
    force_esmiles: bool,
    max_outputs: int,
) -> list[str]:
    if not force_esmiles and not annotation_variants and not ext:
        return smiles

    annotations = annotation_variants or [""]
    outputs: list[str] = []
    for smile in smiles:
        for annotation in annotations:
            outputs.append(Translator.build_esmi(smile, annotation, ext))
            if len(outputs) > max_outputs:
                raise ValueError(
                    f"Markush expansion exceeded max_outputs={max_outputs}"
                )
    return outputs


def _format_preserved_group_outputs(
    states: list[Chem.rdchem.RWMol],
    preserved_groups: str,
    source_rings: tuple[tuple[int, ...], ...],
    annotation_variants: list[str],
    ext: str,
    max_outputs: int,
) -> list[str]:
    outputs: set[str] = set()
    annotations = annotation_variants or [""]
    for state in states:
        smiles = _canonical_smiles(state)
        for annotation in annotations:
            remap_mol = Chem.Mol(state)
            Chem.SanitizeMol(remap_mol)
            remapped_groups = chem_utils.remap_groups(
                remap_mol,
                preserved_groups + annotation,
                source_rings,
            )
            outputs.add(Translator.build_esmi(smiles, remapped_groups, ext))
            if len(outputs) > max_outputs:
                raise ValueError(
                    f"Markush expansion exceeded max_outputs={max_outputs}"
                )
    return sorted(outputs)


def _substitute_precompat_records(
    groups: str,
    definitions: Mapping[str, DefinitionValue],
    *,
    max_outputs: int,
    error_msg: bool,
    repeat_policy: RepeatPolicy,
    terminal_policy: TerminalPolicy,
) -> list[str]:
    variants = [""]
    found_record = False
    for match in _PRECOMPAT_RECORD_PATTERN.finditer(groups):
        found_record = True
        record = match.group(0)
        substruct_match = _SUBSTRUCT_RECORD_PATTERN.match(record)
        if substruct_match:
            body = substruct_match.group("body")
            if repeat_policy == "strict" and (
                "|Sg:" in body or "<g>" in body
            ):
                raise ValueError(
                    "Nested <s> repeat has no unambiguous outer connection semantics"
                )
            body_variants = _substitute_markush_outputs(
                body,
                definitions,
                max_outputs=max_outputs,
                error_msg=error_msg,
                force_esmiles=True,
                repeat_policy=repeat_policy,
                terminal_policy=terminal_policy,
                allow_physical_repeats=False,
            )
            record_variants = [f"<s>{body_variant}</s>" for body_variant in body_variants]
        else:
            record_variants = [_substitute_sgroup_counts(record, definitions)]

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
    repeat_policy: RepeatPolicy = "best_effort",
    terminal_policy: TerminalPolicy = "preserve",
    allow_physical_repeats: bool = True,
) -> list[str]:
    raw_caption = str(caption).strip()
    parsed = Translator.parse_caption(raw_caption, return_mol=True, error_msg=error_msg)
    if parsed is None:
        raise ValueError(f"Invalid E-SMILES caption: {caption}")

    mol, groups, ext = parsed
    ext = _substitute_sgroup_counts(ext, definitions)
    groups = Translator.repair_atom_group_indices(mol, groups, error_msg=error_msg)
    source_rings = mol.GetRingInfo().AtomRings()
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)

    states = [Chem.RWMol(mol)]
    unresolved_group_records: list[str] = []
    for desc in Translator.parse_groups(groups):
        if desc.is_dummy or desc.is_circle:
            continue
        if _is_special_id_label(desc):
            continue
        if isinstance(desc.id, RingIndex) and desc.id.virtual:
            continue
        if isinstance(desc.id, AtomIndex) and _is_carbon_chain_repeat(desc):
            try:
                states = _expand_atom_repetition(
                    states,
                    desc,
                    definitions,
                    max_outputs,
                )
            except _UnexpandableRepeat:
                if repeat_policy != "best_effort":
                    raise
                unresolved_group_records.append(_group_record(desc))
            except ValueError as exc:
                if (
                    repeat_policy == "best_effort"
                    and "No repetition count found" in str(exc)
                ):
                    unresolved_group_records.append(_group_record(desc))
                else:
                    raise
            continue

        fragments = _resolve_group_smiles(desc, definitions)
        if not fragments:
            if repeat_policy == "best_effort":
                unresolved_group_records.append(_group_record(desc))
                continue
            raise ValueError(f"No Markush definition found for group `{str(desc)}`")

        if isinstance(desc.id, AtomIndex):
            copied_states = None
            if desc.multiple:
                try:
                    copied_states = _expand_atom_group_copies(
                        states,
                        desc,
                        fragments,
                        definitions,
                        source_rings,
                        max_outputs,
                    )
                except _UnexpandableRepeat:
                    if repeat_policy == "best_effort":
                        unresolved_group_records.append(_group_record(desc))
                        continue
                    raise
                if copied_states is None:
                    if repeat_policy == "best_effort":
                        unresolved_group_records.append(_group_record(desc))
                        continue
                    raise ValueError(
                        f"No repetition count found for group `{str(desc)}`"
                    )
            states = copied_states if copied_states is not None else _expand_atom_group(
                states,
                desc,
                fragments,
                max_outputs,
            )
        elif isinstance(desc.id, RingIndex):
            if int(desc.id) >= len(source_rings):
                raise ValueError(f"Ring index {int(desc.id)} is out of range")
            expanded_ring_states = _expand_ring_group(
                states,
                desc,
                fragments,
                source_rings,
                definitions,
                max_outputs=max_outputs,
            )
            if expanded_ring_states is None:
                if repeat_policy == "best_effort":
                    unresolved_group_records.append(_group_record(desc))
                    continue
                raise ValueError(
                    f"No repetition count found for group `{str(desc)}`"
                )
            states = expanded_ring_states

    if allow_physical_repeats:
        states, groups, ext = _apply_physical_repeats(
            states,
            groups,
            ext,
            definitions,
            repeat_policy=repeat_policy,
            terminal_policy=terminal_policy,
            max_outputs=max_outputs,
            has_unresolved_groups=bool(unresolved_group_records),
        )
    annotation_variants = _substitute_precompat_records(
        groups,
        definitions,
        max_outputs=max_outputs,
        error_msg=error_msg,
        repeat_policy=repeat_policy,
        terminal_policy=terminal_policy,
    )
    preserved_atom_groups = (
        _preserved_top_level_atom_records(groups)
        + "".join(record for record in unresolved_group_records if record)
    )
    smiles = _collect_valid_smiles(states)
    if not smiles:
        raise ValueError("No valid SMILES generated from Markush substitution")
    if preserved_atom_groups or annotation_variants:
        try:
            return _format_preserved_group_outputs(
                states,
                preserved_atom_groups,
                source_rings,
                annotation_variants,
                ext,
                max_outputs,
            )
        except chem_utils.UnmappableAnnotationError:
            if repeat_policy == "best_effort":
                # The successful substitution removed an atom still used by a
                # residual annotation.  Revert the complete branch instead of
                # returning an E-SMILES record that silently points elsewhere.
                return [_substitute_sgroup_counts(raw_caption, definitions)]
            raise
    return _format_substitution_outputs(
        smiles,
        annotation_variants,
        ext,
        force_esmiles,
        max_outputs,
    )


def substitute_markush(
    caption: str,
    definitions: Mapping[str, DefinitionValue],
    *,
    max_outputs: int = 1024,
    error_msg: bool = False,
    repeat_policy: RepeatPolicy = "best_effort",
    terminal_policy: TerminalPolicy = "preserve",
) -> str | list[str]:
    """Substitute Markush definitions into an E-SMILES caption.

    ``definitions`` maps labels such as ``R1`` or ``R[1]`` to either an
    abbreviation (``Me``) or a SMILES fragment. If the fragment contains ``*``,
    that atom is treated as the attachment point and removed during merging.
    Ring-indexed groups are expanded over possible ring attachment positions.

    ``repeat_policy="best_effort"`` physically expands unambiguous concrete
    Whole-SRU and local ``<g>`` repeats. A fully resolved branch is returned as
    RDKit-valid plain SMILES; a branch with residual annotations remains
    E-SMILES. ``"preserve"`` keeps repeat annotations, while ``"strict"``
    rejects an active repeat that cannot be expanded without guessing.
    ``terminal_policy="hydrogen"`` caps consumed terminal dummy atoms with
    implicit hydrogens; the default retains terminal ``*`` atoms.
    """
    if (
        isinstance(max_outputs, bool)
        or not isinstance(max_outputs, int)
        or max_outputs < 1
    ):
        raise ValueError("max_outputs must be one positive integer")
    if repeat_policy not in {"preserve", "best_effort", "strict"}:
        raise ValueError(
            "repeat_policy must be 'preserve', 'best_effort', or 'strict'"
        )
    if terminal_policy not in {"preserve", "hydrogen"}:
        raise ValueError("terminal_policy must be 'preserve' or 'hydrogen'")
    outputs = _substitute_markush_outputs(
        str(caption).strip(),
        definitions,
        max_outputs=max_outputs,
        error_msg=error_msg,
        repeat_policy=repeat_policy,
        terminal_policy=terminal_policy,
    )
    if len(outputs) == 1:
        return outputs[0]
    return outputs


__all__ = ["substitute_markush"]
