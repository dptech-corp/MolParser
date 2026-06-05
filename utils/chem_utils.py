"""Helpers used by ``translator.py`` for E-SMILES caption refactoring."""

from __future__ import annotations

import ast
import csv
import re
from pathlib import Path
from typing import Dict, Tuple

from rdkit import Chem


_ABBREV_CSV = Path(__file__).parent / "abbrevs_example.csv"


def _load_abbrev_smi() -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with _ABBREV_CSV.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            symbol = row.get("symbol")
            smiles = row.get("smiles")
            if symbol is None or smiles is None:
                continue
            mapping[symbol] = smiles
    return mapping


_abbrev_smi: Dict[str, str] = _load_abbrev_smi()


def get_abbrev_smi() -> Dict[str, str]:
    return _abbrev_smi


def get_mol(smi: str) -> Chem.rdchem.Mol:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smi}")
    return mol


def split_groups(groups: str) -> Tuple[Dict[int, str], Dict[int, str]]:
    pattern = r"<(?P<type>r|a)>(?P<content>.*?)</(?P=type)>"
    matches = re.finditer(pattern, groups)
    a_groups: Dict[int, str] = {}
    r_groups: Dict[int, str] = {}

    for m in matches:
        type_ = m.group("type")
        text = m.group("content")
        if ":" in text:
            ind, content = text.split(":", 1)
        else:
            ind, content = "", text
        if type_ == "r":
            r_groups[int(ind)] = content
        else:
            a_groups[int(ind)] = content
    return a_groups, r_groups


def get_groups_str(a_groups: Dict[int, Dict], r_groups: Dict[int, Dict] | None = None) -> str:
    groups_str = ""
    for k, v in sorted(a_groups.items(), key=lambda x: x[0]):
        groups_str += f"<a>{k}:{v['content']}</a>"
    if r_groups is not None:
        for k, v in sorted(r_groups.items(), key=lambda x: x[0]):
            groups_str += f"<r>{k}:{v['content']}</r>"
    return groups_str


def _canonical_output_order(mol: Chem.rdchem.Mol) -> Tuple[Dict[int, int], Chem.rdchem.Mol | None]:
    output_mol = Chem.Mol(mol)
    for atom in output_mol.GetAtoms():
        atom.SetAtomMapNum(0)

    smi = Chem.MolToSmiles(output_mol, canonical=True, isomericSmiles=True)
    final_mol = Chem.MolFromSmiles(smi)
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
        src_mol_check = Chem.MolFromSmiles(smiles)
        rep_atom = src_mol_check.GetAtomWithIdx(0)
        atom.SetAtomicNum(rep_atom.GetAtomicNum())
        atom.SetFormalCharge(rep_atom.GetFormalCharge())
    elif element:
        atom.SetAtomicNum(Chem.GetPeriodicTable().GetAtomicNumber(element))
        atom.SetFormalCharge(0)
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
        group = old_ind2agroup.get(old_ind)
        if group is not None:
            new_ind2agroup[new_ind] = {"content": group, "type": "a"}

    new_ind2rgroup: Dict[int, Dict] = {}
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

        for old_ind, group in old_ind2rgroup.items():
            new_ring_idx = ring_map.get(old_ind)
            if new_ring_idx is not None:
                new_ind2rgroup[new_ring_idx] = {"content": group, "type": "r"}

    return get_groups_str(new_ind2agroup, new_ind2rgroup)


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
    "get_groups_str",
    "get_mol",
    "is_single_bond",
    "merge_group",
    "remap_groups",
    "split_groups",
]
