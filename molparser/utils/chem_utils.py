"""Helpers used by ``translator.py`` for E-SMILES caption refactoring."""

from __future__ import annotations

import ast
import csv
import re
from pathlib import Path
from typing import Dict, List, Tuple

from rdkit import Chem


_ABBREV_CSV = Path(__file__).parent / "abbrevs_example.csv"
_IONIC_MONOVALENT_METALS = frozenset({"Li", "Na", "K"})


def parse_smiles(smiles: str, *, sanitize: bool = True) -> Chem.rdchem.Mol | None:
    """Keep explicit atoms so E-SMILES indices survive every parse boundary."""
    params = Chem.SmilesParserParams()
    params.removeHs = False
    params.parseName = False
    params.sanitize = sanitize
    return Chem.MolFromSmiles(str(smiles), params)


def normalize_ionic_smiles(smiles: str) -> str:
    """Charge neutral O-alkali bonds as an ionic salt when unambiguous."""
    mol = parse_smiles(smiles)
    if mol is None:
        return str(smiles)
    editable = Chem.RWMol(mol)
    changed = False
    for bond in editable.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        begin, end = bond.GetBeginAtom(), bond.GetEndAtom()
        if begin.GetSymbol() not in _IONIC_MONOVALENT_METALS and end.GetSymbol() not in _IONIC_MONOVALENT_METALS:
            continue
        oxygen = begin if begin.GetSymbol() == "O" else end if end.GetSymbol() == "O" else None
        metal = begin if begin.GetSymbol() in _IONIC_MONOVALENT_METALS else end if end.GetSymbol() in _IONIC_MONOVALENT_METALS else None
        if oxygen is None or metal is None:
            continue
        if (
            metal.GetDegree() != 1
            or oxygen.GetDegree() > 2
            or oxygen.GetNumExplicitHs() > 0
            or any(
                atom.GetSymbol() in _IONIC_MONOVALENT_METALS and atom.GetIdx() != metal.GetIdx()
                for atom in oxygen.GetNeighbors()
            )
        ):
            continue
        if oxygen.GetFormalCharge() == 0 and metal.GetFormalCharge() == 0:
            oxygen.SetFormalCharge(-1)
            metal.SetFormalCharge(1)
            changed = True
    if not changed:
        return str(smiles)
    try:
        Chem.SanitizeMol(editable)
        # Atom 0 is the implicit attachment site. Canonical output can move it
        # to oxygen or lithium; serialization can also pin its H as in [SH].
        rewritten = Chem.MolToSmiles(editable, canonical=False, isomericSmiles=True)
        round_tripped = parse_smiles(rewritten)
        if round_tripped is None or _attachment_slot(mol) != _attachment_slot(round_tripped):
            return str(smiles)
        return rewritten
    except Exception:
        return str(smiles)


def _attachment_slot(mol: Chem.rdchem.Mol) -> Tuple[str, int]:
    """Track the implicit attachment atom and any pinned explicit hydrogens."""
    if mol.GetNumAtoms() == 0:
        return ("", 0)
    atom = mol.GetAtomWithIdx(0)
    return (atom.GetSymbol(), atom.GetNumExplicitHs() if atom.GetNoImplicit() else 0)


class UnmappableAnnotationError(ValueError):
    """A retained annotation references an atom removed by substitution."""


def _load_abbrev_smi() -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with _ABBREV_CSV.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            symbol = row.get("symbol")
            smiles = row.get("smiles")
            if symbol is None or smiles is None:
                continue
            mapping[symbol] = normalize_ionic_smiles(smiles)
    return mapping


_abbrev_smi: Dict[str, str] = _load_abbrev_smi()


def get_abbrev_smi() -> Dict[str, str]:
    return _abbrev_smi


def get_mol(smi: str) -> Chem.rdchem.Mol:
    mol = parse_smiles(smi)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smi}")
    return mol


_PRESERVED_RECORD_PATTERN = re.compile(
    r"<s>.*?</s>|<g>.*?</g>|<r><v>\d+:.+?</r>|<v>.*?</v>|"
    r"<c>.*?</c>|<x>.*?</x>",
    re.DOTALL,
)
_PORT_PAIR_PATTERN = re.compile(r"\[(?P<left>\d+):(?P<right>\d+)\]")
_ATOM_RECORD_INDEX_PATTERN = re.compile(r"(?<=<c>)(?P<index>\d+)(?=:)")


def _preserved_records(
    groups: str,
    atom_index_map: Dict[int, int] | None = None,
) -> str:
    """Return pre-compatible records, remapping top-level graph indices.

    Nested ``<s>`` records own a separate molecule namespace and are therefore
    kept byte-for-byte.  Top-level local-SRU ports and virtual-arc endpoints,
    however, refer to the outer graph and must follow canonical atom ordering.
    """

    def remap_pair(match: re.Match, ordered: bool = False) -> str:
        if atom_index_map is None:
            return match.group(0)
        left, right = int(match.group("left")), int(match.group("right"))
        if left not in atom_index_map or right not in atom_index_map:
            raise UnmappableAnnotationError(
                "Unable to remap a preserved E-SMILES atom index"
            )
        left, right = atom_index_map[left], atom_index_map[right]
        if ordered:
            left, right = sorted((left, right))
        return f"[{left}:{right}]"

    records: List[str] = []
    for match in _PRESERVED_RECORD_PATTERN.finditer(groups):
        record = match.group(0)
        if record.startswith(("<v>", "<x>")):
            record = _PORT_PAIR_PATTERN.sub(lambda pair: remap_pair(pair, ordered=True), record)
        elif record.startswith("<g>"):
            record = _PORT_PAIR_PATTERN.sub(remap_pair, record)
        elif record.startswith("<c>") and atom_index_map is not None:
            indexed = _ATOM_RECORD_INDEX_PATTERN.search(record)
            old_index = int(indexed.group("index")) if indexed is not None else -1
            if old_index not in atom_index_map:
                raise UnmappableAnnotationError(
                    "Unable to remap a preserved E-SMILES atom index"
                )
            record = _ATOM_RECORD_INDEX_PATTERN.sub(
                str(atom_index_map[old_index]),
                record,
                count=1,
            )
        records.append(record)
    return "".join(records)


def _strip_preserved_records(groups: str) -> str:
    return _PRESERVED_RECORD_PATTERN.sub("", groups)


def _ring_index_key(index_text: str) -> tuple[int, ...] | None:
    if re.fullmatch(r"\d+(?:,\d+)*", index_text) is None:
        return None
    return tuple(dict.fromkeys(int(part) for part in index_text.split(",")))


def split_groups(
    groups: str,
) -> Tuple[Dict[int, Dict[str, str]], Dict[tuple[int, ...], List[str]]]:
    pattern = r"<(?P<type>r|a|d)>(?P<content>.*?)</(?P=type)>"
    matches = re.finditer(pattern, _strip_preserved_records(groups))
    a_groups: Dict[int, Dict[str, str]] = {}
    r_groups: Dict[tuple[int, ...], List[str]] = {}

    for m in matches:
        type_ = m.group("type")
        text = m.group("content")
        if ":" in text:
            ind, content = text.split(":", 1)
        else:
            ind, content = "", text
        indexes = _ring_index_key(ind)
        if indexes is None:
            continue
        if type_ == "r":
            r_groups.setdefault(indexes, []).append(content)
        elif len(indexes) == 1:
            a_groups[indexes[0]] = {"content": content, "type": type_}
    return a_groups, r_groups


def get_groups_str(
    a_groups: Dict[int, Dict],
    r_groups: Dict[tuple[int, ...], Dict] | None = None,
) -> str:
    groups_str = ""
    for k, v in sorted(a_groups.items(), key=lambda x: x[0]):
        tag = v.get("type", "a")
        groups_str += f"<{tag}>{k}:{v['content']}</{tag}>"
    if r_groups is not None:
        for k, v in sorted(r_groups.items(), key=lambda x: x[0]):
            contents = v["content"] if isinstance(v["content"], list) else [v["content"]]
            index_text = ",".join(str(index) for index in k)
            for content in contents:
                groups_str += f"<r>{index_text}:{content}</r>"
    return groups_str


def _canonical_output_order(mol: Chem.rdchem.Mol) -> Tuple[Dict[int, int], Chem.rdchem.Mol | None]:
    output_mol = Chem.Mol(mol)
    for atom in output_mol.GetAtoms():
        atom.SetAtomMapNum(0)

    smi = Chem.MolToSmiles(output_mol, canonical=True, isomericSmiles=True)
    final_mol = parse_smiles(smi)
    try:
        order = ast.literal_eval(output_mol.GetProp("_smilesAtomOutputOrder"))
    except Exception:
        order = list(range(output_mol.GetNumAtoms()))
    return {internal_idx: output_idx for output_idx, internal_idx in enumerate(order)}, final_mol


def is_single_bond(bond: Chem.rdchem.Bond) -> bool:
    return bond.GetBondType() == Chem.BondType.SINGLE


def alter_atom(
    atom: Chem.rdchem.Atom,
    smiles: None | str,
    element: None | str = None,
) -> None:
    if smiles:
        src_mol_check = parse_smiles(smiles)
        rep_atom = src_mol_check.GetAtomWithIdx(0)
        atom.SetAtomicNum(rep_atom.GetAtomicNum())
        atom.SetFormalCharge(rep_atom.GetFormalCharge())
        atom.SetIsotope(rep_atom.GetIsotope())
    elif element:
        atom.SetAtomicNum(Chem.GetPeriodicTable().GetAtomicNumber(element))
        atom.SetFormalCharge(0)
        atom.SetIsotope(0)
    else:
        raise ValueError("Either smiles or element must be provided")

    atom.SetNumExplicitHs(0)
    atom.SetNoImplicit(False)
    if atom.GetDegree() == 1:
        atom.SetIsAromatic(False)


def merge_group(tgt_mol: Chem.rdchem.RWMol, src: str, attach_idx: int) -> None:
    dummy_atom = tgt_mol.GetAtomWithIdx(attach_idx)
    assert dummy_atom.GetSymbol() == "*", "Non-dummy attachment point"

    src_mol = get_mol(src)

    neighbors = dummy_atom.GetNeighbors()
    assert len(neighbors) == 1, "Attachment point must have exactly one neighbor"
    nb_idx = neighbors[0].GetIdx()

    conn_bond = tgt_mol.GetBondBetweenAtoms(attach_idx, nb_idx)
    assert conn_bond.GetBondType() == Chem.BondType.SINGLE, (
        "Attachment point must link to single bond"
    )

    idx_map: Dict[int, int] = {}
    for atom in src_mol.GetAtoms():
        i = tgt_mol.AddAtom(atom)
        idx_map[atom.GetIdx()] = i
        tgt_atom = tgt_mol.GetAtomWithIdx(i)
        tgt_atom.SetChiralTag(atom.GetChiralTag())
    for bond in src_mol.GetBonds():
        j = idx_map[bond.GetBeginAtomIdx()]
        k = idx_map[bond.GetEndAtomIdx()]
        tgt_mol.AddBond(j, k, bond.GetBondType())
    tgt_mol.RemoveBond(attach_idx, nb_idx)
    tgt_mol.AddBond(nb_idx, idx_map[0], Chem.BondType.SINGLE)


def remap_groups(mol: Chem.rdchem.Mol, groups: str, ring_info: tuple) -> str:
    old_ind2agroup, old_ind2rgroup = split_groups(groups)
    new_ind2agroup: Dict[int, Dict] = {}
    internal_ind_map: Dict[int, int] = {}
    old_to_output: Dict[int, int] = {}
    internal_to_output, output_mol = _canonical_output_order(mol)

    for atom in mol.GetAtoms():
        s = atom.GetSmarts()
        atom.SetAtomMapNum(0)
        if ":" not in s:
            continue
        try:
            old_ind = int(s.split(":")[1].replace("]", "")) - 1
        except Exception:
            continue
        internal_ind = atom.GetIdx()
        internal_ind_map[old_ind] = internal_ind
        new_ind = internal_to_output.get(internal_ind, internal_ind)
        old_to_output[old_ind] = new_ind
        group = old_ind2agroup.get(old_ind)
        if group is not None:
            new_ind2agroup[new_ind] = group

    new_ind2rgroup: Dict[tuple[int, ...], Dict] = {}
    if old_ind2rgroup:
        new_ring_info = output_mol.GetRingInfo().AtomRings() if output_mol is not None else ()
        new_ring_sets = [set(ring) for ring in new_ring_info]
        ring_map: Dict[int, int | None] = {}

        for orig_idx, orig_ring in enumerate(ring_info):
            mapped = {
                internal_to_output.get(internal_ind_map.get(atom, atom), internal_ind_map.get(atom, atom))
                for atom in orig_ring
            }
            best_match = None
            best_overlap = 0
            for new_idx, new_ring in enumerate(new_ring_sets):
                overlap = len(mapped.intersection(new_ring))
                if overlap > best_overlap:
                    best_overlap = overlap
                    best_match = new_idx
            ring_map[orig_idx] = best_match if best_overlap > 0 else None

        for old_inds, group_list in old_ind2rgroup.items():
            new_inds: list[int] = []
            for old_ind in old_inds:
                new_ring_idx = ring_map.get(old_ind)
                if new_ring_idx is None:
                    new_inds = []
                    break
                if new_ring_idx not in new_inds:
                    new_inds.append(new_ring_idx)
            if not new_inds:
                continue
            entry = new_ind2rgroup.setdefault(
                tuple(new_inds), {"content": [], "type": "r"}
            )
            entry["content"].extend(group_list)

    return (
        get_groups_str(new_ind2agroup, new_ind2rgroup)
        + _preserved_records(groups, old_to_output)
    )


def carbon_chain_repetition_process(mol, atom_id, desc, is_markush, error_msg=False):
    atom = mol.GetAtomWithIdx(atom_id)
    count = 1
    if desc.multiple and desc.multiple.isdigit():
        count = int(desc.multiple)

    try:
        atom.SetAtomicNum(6)
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(False)
        atom.SetIsAromatic(False)
        atom.SetIsotope(0)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

        for bond in atom.GetBonds():
            bond.SetBondType(Chem.BondType.SINGLE)
            bond.SetIsAromatic(False)

        if count > 1:
            nbrs = atom.GetNeighbors()
            degree = len(nbrs)

            if degree == 2:
                nb_right = nbrs[1]
                right_idx = nb_right.GetIdx()
                mol.RemoveBond(atom_id, right_idx)

                current_end_idx = atom_id
                for _ in range(count - 1):
                    new_atom = Chem.Atom(6)
                    new_atom.SetIsAromatic(False)
                    new_idx = mol.AddAtom(new_atom)
                    mol.AddBond(current_end_idx, new_idx, Chem.BondType.SINGLE)
                    current_end_idx = new_idx
                mol.AddBond(current_end_idx, right_idx, Chem.BondType.SINGLE)

            elif degree == 1:
                current_end_idx = atom_id
                for _ in range(count - 1):
                    new_atom = Chem.Atom(6)
                    new_atom.SetIsAromatic(False)
                    new_idx = mol.AddAtom(new_atom)
                    mol.AddBond(current_end_idx, new_idx, Chem.BondType.SINGLE)
                    current_end_idx = new_idx

            else:
                is_markush = True

    except Exception:
        is_markush = True

    return is_markush


__all__ = [
    "alter_atom",
    "carbon_chain_repetition_process",
    "get_abbrev_smi",
    "normalize_ionic_smiles",
    "get_groups_str",
    "get_mol",
    "is_single_bond",
    "merge_group",
    "remap_groups",
    "split_groups",
]
