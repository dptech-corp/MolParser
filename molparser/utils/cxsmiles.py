"""Best-effort conversion from E-SMILES captions to CXSMILES."""

from __future__ import annotations

import ast
import re

from rdkit import Chem

try:
    from .translator import AtomIndex, Translator
except ImportError:  # Support running from package directory as working directory.
    from translator import AtomIndex, Translator


_RAW_GROUP_PATTERN = re.compile(r"<(?P<tag>a|c|r)>(?P<idx>\d+):(?P<label>.+?)</(?P=tag)>")


def _raw_groups(groups: str) -> list[tuple[str, int, str]]:
    return [
        (match.group("tag"), int(match.group("idx")), match.group("label"))
        for match in _RAW_GROUP_PATTERN.finditer(groups)
    ]


def _caption_groups(caption: str) -> str:
    if "<sep>" not in caption:
        return ""
    trailing = caption.split("<sep>", 1)[1]
    groups, _ = Translator.parse_trailing(trailing)
    return groups


def _apply_atom_labels(mol: Chem.rdchem.Mol, groups: str) -> None:
    """Preserve atom-indexed E-SMILES labels when CXSMILES can carry them."""
    for desc in Translator.parse_groups(groups):
        if not isinstance(desc.id, AtomIndex):
            continue
        atom_idx = int(desc.id)
        if atom_idx >= mol.GetNumAtoms() or desc.is_dummy:
            continue
        label = str(desc)
        if label:
            mol.GetAtomWithIdx(atom_idx).SetProp("atomLabel", label)


def _apply_raw_atom_labels(mol: Chem.rdchem.Mol, groups: str) -> None:
    for tag, atom_idx, label in _raw_groups(groups):
        if atom_idx >= mol.GetNumAtoms():
            continue
        if tag == "c":
            mol.GetAtomWithIdx(atom_idx).SetProp("atomLabel", f"c{label}")
        elif tag == "a" and label.startswith("<") and label != "<dum>":
            mol.GetAtomWithIdx(atom_idx).SetProp("atomLabel", label)


def _apply_sru_sgroup(mol: Chem.rdchem.Mol) -> None:
    if mol.GetNumAtoms() == 0:
        return
    sgroup = Chem.CreateMolSubstanceGroup(mol, "SRU")
    sgroup.SetAtoms(list(range(mol.GetNumAtoms())))
    sgroup.SetProp("LABEL", "n")
    sgroup.SetProp("CONNECT", "HT")


def _add_ring_attachments(
    mol: Chem.rdchem.Mol,
    groups: str,
) -> tuple[Chem.rdchem.Mol, list[tuple[int, tuple[int, ...]]]]:
    """Represent ring-indexed groups as variable attachment dummy atoms."""
    ring_info = mol.GetRingInfo().AtomRings()
    attachments: list[tuple[int, tuple[int, ...]]] = []
    ring_counts: dict[int, int] = {}
    rw_mol = Chem.RWMol(mol)

    for tag, ring_idx, label in _raw_groups(groups):
        if tag != "r":
            continue
        if ring_idx >= len(ring_info):
            continue
        ring_atoms = tuple(ring_info[ring_idx])
        if not ring_atoms:
            continue

        atom = Chem.Atom(0)
        if label:
            atom.SetProp("atomLabel", label)
        dummy_idx = rw_mol.AddAtom(atom)
        anchor_idx = ring_atoms[ring_counts.get(ring_idx, 0) % len(ring_atoms)]
        ring_counts[ring_idx] = ring_counts.get(ring_idx, 0) + 1
        rw_mol.AddBond(dummy_idx, anchor_idx, Chem.BondType.SINGLE)
        attachments.append((dummy_idx, ring_atoms))

    if not attachments:
        return mol, attachments
    return rw_mol.GetMol(), attachments


def _append_cx_fields(cxsmiles: str, fields: list[str]) -> str:
    if not fields:
        return cxsmiles
    if " |" not in cxsmiles:
        return f"{cxsmiles} |{','.join(fields)}|"

    prefix, block = cxsmiles.rsplit(" |", 1)
    if not block.endswith("|"):
        return f"{cxsmiles} |{','.join(fields)}|"
    return f"{prefix} |{block[:-1]},{','.join(fields)}|"


def _ring_attachment_fields(
    mol: Chem.rdchem.Mol,
    attachments: list[tuple[int, tuple[int, ...]]],
) -> list[str]:
    if not attachments:
        return []

    Chem.MolToSmiles(mol, canonical=False, isomericSmiles=True)
    try:
        order = ast.literal_eval(mol.GetProp("_smilesAtomOutputOrder"))
    except Exception:
        order = list(range(mol.GetNumAtoms()))
    output_idx = {atom_idx: idx for idx, atom_idx in enumerate(order)}

    fields = []
    for dummy_idx, ring_atoms in attachments:
        if dummy_idx not in output_idx:
            continue
        targets = [output_idx[atom_idx] for atom_idx in ring_atoms if atom_idx in output_idx]
        if targets:
            fields.append(
                f"m:{output_idx[dummy_idx]}:{'.'.join(str(idx) for idx in targets)}"
            )
    return fields


def _convert_refactored_esmi_to_cxsmiles(
    esmi: str,
    source_groups: str = "",
    sru: bool = False,
) -> str:
    """Convert a Translator-refactored E-SMILES string to CXSMILES."""
    parsed = Translator.parse_caption(esmi)
    if parsed is None:
        return esmi.split("<sep>", 1)[0]
    smi, groups, ext = parsed

    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return smi

    ring_groups = source_groups if "<r>" in source_groups else groups
    mol, ring_attachments = _add_ring_attachments(mol, ring_groups)
    _apply_atom_labels(mol, groups)
    _apply_raw_atom_labels(mol, source_groups or groups)
    if sru or Translator.parse_extension(ext) == "Sg:n":
        _apply_sru_sgroup(mol)

    params = Chem.SmilesWriteParams()
    params.canonical = not ring_attachments
    params.doIsomericSmiles = True
    cxsmiles = Chem.MolToCXSmiles(mol, params, Chem.CXSmilesFields.CX_ALL)
    return _append_cx_fields(cxsmiles, _ring_attachment_fields(mol, ring_attachments))
