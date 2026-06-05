"""E-SMILES caption parsing & refactoring used by the inference postprocess."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from enum import Enum, unique
from typing import Dict, List, Optional, Tuple, Union

from rdkit import Chem, RDLogger

try:
    from . import chem_utils
except ImportError:  # Support running from package directory as working directory.
    import chem_utils


# Ambiguous element symbols are omitted.
PERIODIC_TABLE = {
    "H", "He", "Li", "Be", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Th",
    "Pa", "U",
}


LONG_CHAR_ELEMENTS = (
    "Cl", "Br", "Si", "Na", "Al", "Ca", "Sn", "As", "Hg",
    "Fe", "Zn", "Cr", "Se", "Gd", "Au", "Li",
)


logger = logging.getLogger(__name__)


@unique
class TextType(Enum):
    SYMBOL = "symbol"
    SCRIPT = "script"
    MULTIPLE = "multiple"
    PRIME = "prime"


class Index(int):
    def __new__(cls, value: int, **kwds):
        assert isinstance(value, int) and value >= 0
        return super().__new__(cls, value)


class AtomIndex(Index):
    key: int = 1

    def __hash__(self):
        return hash(("atom", int(self)))

    def __eq__(self, other):
        if isinstance(other, AtomIndex):
            return hash(self) == hash(other)
        return False


class RingIndex(Index):
    def __init__(self, value: int, virtual: bool = False) -> None:
        self._virtual = virtual
        self._key = 2_000_000 if self._virtual else 1_000_000

    def __hash__(self):
        return hash(("ring", int(self), self._virtual))

    def __eq__(self, other):
        if isinstance(other, RingIndex):
            return hash(self) == hash(other)
        return False

    @property
    def virtual(self) -> bool:
        return self._virtual

    @property
    def key(self) -> bool | int:
        return self._key


class Tokens:
    """Special token names used in the captioning format."""

    atom_start = "<a>"
    atom_end = "</a>"
    circ_start = "<c>"
    circ_end = "</c>"
    ring_start = "<r>"
    ring_end = "</r>"
    dummy = "<dum>"
    separator = "<sep>"


class Patterns:
    """Compiled regexes for caption parsing."""

    long_char_elements_pattern = re.compile(
        rf'{"|".join(LONG_CHAR_ELEMENTS) + "|."}'
    )
    grp_content = re.compile(
        rf"(?P<{TextType.SYMBOL.value}>[A-Za-z0-9-\(\)]*)"
        + rf"(?P<{TextType.SCRIPT.value}>(\[\S+\])?)"
        + rf"(?P<{TextType.PRIME.value}>[\'\"]?)"
        + rf"(?P<{TextType.MULTIPLE.value}>(\?([a-z]|\d+|\d-\d)$)?)"
    )
    grp_pattern = re.compile(
        rf"({Tokens.atom_start}|{Tokens.circ_start}|{Tokens.ring_start}|{Tokens.ring_start}{Tokens.circ_start})"
        + r"(\d+:\S+?)"
        + rf"({Tokens.atom_end}|{Tokens.circ_end}|{Tokens.ring_end})"
    )
    trail_pattern = re.compile(
        r"(?P<groups>([^|]*)?)(?P<extension>(\|\S+\|)?$)"
    )


@dataclass
class GroupDesc:
    id: Index
    symbol: Optional[str] = None
    script: Optional[str] = None
    prime: Optional[str] = None
    multiple: Optional[str] = None
    is_circle: bool = False
    is_dummy: bool = False

    def __str__(self) -> str:
        if self.is_dummy:
            return "<dum>"
        expr = ""
        if self.is_circle:
            expr = "c"
        if self.symbol is not None:
            expr += self.symbol
        if self.script is not None:
            expr = expr + "[" + self.script + "]"
        if self.prime is not None:
            expr += self.prime
        if self.multiple is not None:
            expr = expr + "?" + self.multiple
        return expr


@dataclass(frozen=True)
class TranslatedMolecule:
    smi: str
    groups: str
    caption: str
    esmi: str
    markush: bool
    sru: bool


class Translator:
    """Caption parsing and refactor (E-SMILES) for inference postprocess."""

    @classmethod
    def canonicalize_smiles(cls, smi: str) -> str:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                return smi
            return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        except Exception:
            return smi

    @classmethod
    def build_esmi(cls, smi: str, groups: str = "", ext: str = "") -> str:
        return f"{smi}{Tokens.separator}{groups}{ext}"

    @classmethod
    def remove_atom_groups(cls, groups: str, atom_indices: set[int]) -> str:
        if not atom_indices:
            return groups
        atom_group_pattern = re.compile(rf"{Tokens.atom_start}(\d+):.+?{Tokens.atom_end}")

        def replace(match: re.Match) -> str:
            return "" if int(match.group(1)) in atom_indices else match.group(0)

        return atom_group_pattern.sub(replace, groups)

    @classmethod
    def parse_caption(
        cls,
        caption: str,
        return_mol: bool = False,
        error_msg: bool = False,
    ) -> Optional[Tuple[Union[Chem.rdchem.Mol, str], str, str]]:
        smi, *trailings = caption.split(Tokens.separator)
        if len(trailings) != 1:
            if error_msg:
                logger.warning(
                    f"{len(trailings)} `{Tokens.separator}` found in caption: {caption}"
                )
            return

        if error_msg:
            RDLogger.EnableLog("rdApp.*")
        mol = Chem.MolFromSmiles(smi)
        RDLogger.DisableLog("rdApp.*")
        if mol is None:
            if error_msg:
                logger.warning(f"Invalid SMILES: {smi}")
            return
        groups, ext = cls.parse_trailing(trailings[0])
        if return_mol:
            return mol, groups, ext
        return smi, groups, ext

    @classmethod
    def parse_trailing(cls, trailing: str):
        matched = re.match(Patterns.trail_pattern, trailing)
        if matched is None:
            return "", ""
        content = matched.groupdict()
        return content.get("groups", ""), content.get("extension", "")

    @classmethod
    def parse_extension(cls, ext: str) -> str:
        return ext.strip("|")

    @classmethod
    def parse_groups(cls, seq: str) -> List[GroupDesc]:
        if seq == "":
            return []
        descriptions: List[GroupDesc] = []
        for grp_start, grp_content, _ in re.findall(Patterns.grp_pattern, seq):
            parsed = cls.parse_group(grp_content)
            if parsed is None:
                continue
            idx, grp_text = parsed
            if grp_start == Tokens.atom_start:
                grp_desc = GroupDesc(id=AtomIndex(idx))
                if len(grp_text) == 0:
                    grp_desc.is_dummy = True
            elif grp_start == Tokens.circ_start:
                grp_desc = GroupDesc(id=AtomIndex(idx), is_circle=True)
            elif grp_start == f"{Tokens.ring_start}{Tokens.circ_start}":
                grp_desc = GroupDesc(id=RingIndex(idx, virtual=True))
            elif grp_start == Tokens.ring_start:
                grp_desc = GroupDesc(id=RingIndex(idx))
            else:
                continue
            grp_desc.symbol = grp_text.get(TextType.SYMBOL)
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
        if content == Tokens.dummy:
            return idx, {}
        grp_text = cls.get_group_texts(content)
        if len(grp_text) == 0:
            return
        return idx, grp_text

    @classmethod
    def get_group_texts(cls, content: str) -> Dict[TextType, str]:
        texts: Dict[TextType, str] = {}
        matched = re.match(Patterns.grp_content, content).groupdict()
        for tt_value, text in matched.items():
            if len(text) == 0:
                continue
            if tt_value == TextType.SYMBOL.value:
                texts[TextType.SYMBOL] = text
            elif tt_value == TextType.SCRIPT.value:
                assert text.startswith("[") and text.endswith("]")
                texts[TextType.SCRIPT] = text[1:-1]
            elif tt_value == TextType.PRIME.value:
                texts[TextType.PRIME] = text
            elif tt_value == TextType.MULTIPLE.value:
                assert text.startswith("?")
                texts[TextType.MULTIPLE] = text[1:]
        return texts

    @classmethod
    def repair_atom_group_indices(
        cls,
        mol: Chem.rdchem.Mol,
        trailing: str,
        error_msg: bool = False,
    ) -> str:
        """Repair atom-group tags that miss their dummy atom."""
        star_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "*"]
        if not star_indices or Tokens.atom_start not in trailing:
            return trailing

        atom_group_pattern = re.compile(rf"{Tokens.atom_start}(\d+):(.+?){Tokens.atom_end}")
        matches = list(atom_group_pattern.finditer(trailing))
        if not matches:
            return trailing

        exact_targets = {
            int(match.group(1))
            for match in matches
            if int(match.group(1)) in star_indices
        }
        used_targets: set[int] = set()

        def replace(match: re.Match) -> str:
            raw_idx, content = match.group(1), match.group(2)
            idx = int(raw_idx)
            if idx in star_indices and idx not in used_targets:
                used_targets.add(idx)
                return match.group(0)

            candidates = [
                star_idx
                for star_idx in star_indices
                if star_idx not in used_targets and star_idx not in exact_targets
            ]
            if not candidates:
                return match.group(0)

            nearest = min(candidates, key=lambda star_idx: (abs(star_idx - idx), star_idx))
            used_targets.add(nearest)
            if error_msg:
                logger.warning(
                    "Repair atom group index: <a>%s:%s</a> -> <a>%s:%s</a>",
                    raw_idx,
                    content,
                    nearest,
                    content,
                )
            return f"{Tokens.atom_start}{nearest}:{content}{Tokens.atom_end}"

        return atom_group_pattern.sub(replace, trailing)

    @classmethod
    def refactor(
        cls,
        caption: str,
        error_msg: bool = False,
    ) -> Optional[TranslatedMolecule]:
        """Refactor E-SMILES and detect Markush/SRU output."""
        if "<sep>" not in caption:
            canonical_smi = cls.canonicalize_smiles(caption)
            return TranslatedMolecule(
                smi=canonical_smi,
                groups="",
                caption=caption,
                esmi=cls.build_esmi(canonical_smi),
                markush=False,
                sru=False,
            )
        if caption.split("<sep>")[-1] == "":
            smi = caption.split("<sep>")[0]
            canonical_smi = cls.canonicalize_smiles(smi)
            return TranslatedMolecule(
                smi=canonical_smi,
                groups="",
                caption=caption,
                esmi=cls.build_esmi(canonical_smi),
                markush=False,
                sru=False,
            )

        smi = caption.split("<sep>")[0]
        trailing = caption.split("<sep>", 1)[1]
        groups, ext = cls.parse_trailing(trailing)

        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                if error_msg:
                    logger.warning(f"Invalid SMILES: {smi}")
                return TranslatedMolecule(
                    smi=smi,
                    groups=groups,
                    caption=caption,
                    esmi=cls.build_esmi(smi, groups, ext),
                    markush=len(groups) > 0,
                    sru=False,
                )
            repaired_groups = cls.repair_atom_group_indices(mol, groups, error_msg=error_msg)
            if repaired_groups != groups:
                groups = repaired_groups
                caption = cls.build_esmi(smi, groups, ext)
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(atom.GetIdx() + 1)
            ring_info = mol.GetRingInfo().AtomRings()
            mapped_smiles = Chem.MolToSmiles(mol, canonical=False, isomericSmiles=True)
            mapped_caption = cls.build_esmi(mapped_smiles, groups, ext)

            parsed = cls.parse_caption(
                mapped_caption, return_mol=True, error_msg=error_msg
            )
            mol, groups, ext = parsed
            mol = Chem.RWMol(mol)

            mapnum2idx: Dict[int, int] = {}
            for a in mol.GetAtoms():
                mn = a.GetAtomMapNum()
                if mn > 0:
                    mapnum2idx[mn] = a.GetIdx()

            to_remove: List[int] = []
            consumed_atom_groups: set[int] = set()
            preserve_dummy_groups = False
            is_markush = False
            is_sru = False

            if cls.parse_extension(ext) == "Sg:n":
                is_sru = True

            abbrev_map = chem_utils.get_abbrev_smi()

            for desc in cls.parse_groups(groups):
                if not isinstance(desc.id, AtomIndex) or desc.is_circle:
                    is_markush = True
                    continue

                atom_idx = mapnum2idx.get(int(desc.id) + 1)
                if atom_idx is None:
                    is_markush = True
                    continue

                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() != "*" and not desc.is_dummy:
                    is_markush = True
                    continue

                if desc.is_dummy:
                    preserve_dummy_groups = True
                    continue

                if not desc.symbol:
                    continue

                # Carbon chain repetition, e.g. (CH2)n.
                if (desc.symbol == "CH2") or (desc.symbol == "CH" and desc.script == "2"):
                    if desc.multiple and not desc.multiple.isdigit():
                        is_markush = True
                        continue
                    is_markush = chem_utils.carbon_chain_repetition_process(
                        mol, atom_idx, desc, is_markush, error_msg=False
                    )
                    consumed_atom_groups.add(int(desc.id))
                    continue

                # Build lookup key, e.g. NO + 2 -> NO2.
                lookup_symbol = desc.symbol
                if desc.script:
                    lookup_symbol += desc.script

                if lookup_symbol == "CN":
                    src_smi = "C(#N)"
                else:
                    src_smi = abbrev_map.get(lookup_symbol)

                if src_smi is not None:
                    src_mol_check = Chem.MolFromSmiles(src_smi)
                    if src_mol_check and src_mol_check.GetNumAtoms() == 1:
                        chem_utils.alter_atom(atom, smiles=src_smi)
                        consumed_atom_groups.add(int(desc.id))
                        continue

                    if atom.GetDegree() != 1:
                        if error_msg:
                            logger.warning(
                                f"Group `{lookup_symbol}` cannot be attached to atom "
                                f"{atom_idx} with degree {atom.GetDegree()}"
                            )
                        is_markush = True
                        continue

                    if any(not chem_utils.is_single_bond(b) for b in atom.GetBonds()):
                        if error_msg:
                            logger.warning(
                                f"Group `{lookup_symbol}` must link to single bond"
                            )
                        is_markush = True
                        continue

                    try:
                        chem_utils.merge_group(tgt_mol=mol, src=src_smi, attach_idx=atom_idx)
                        to_remove.append(atom_idx)
                        consumed_atom_groups.add(int(desc.id))
                    except Exception as e:
                        if error_msg:
                            logger.warning(f"Merge group failed: {e}")
                        is_markush = True
                    continue

                if lookup_symbol in PERIODIC_TABLE and not desc.multiple:
                    try:
                        chem_utils.alter_atom(atom, smiles=None, element=lookup_symbol)
                        consumed_atom_groups.add(int(desc.id))
                    except Exception as e:
                        if error_msg:
                            logger.warning(
                                f"Failed to mutate atom {atom_idx} to {lookup_symbol}: {e}"
                            )
                        is_markush = True
                    continue

                is_markush = True

            for i in sorted(to_remove, reverse=True):
                mol.RemoveAtom(i)

            try:
                Chem.SanitizeMol(mol)
            except Exception as e:
                if error_msg:
                    logger.error(f"Sanitize failed: {e}")
                return TranslatedMolecule(
                    smi=smi,
                    groups=groups,
                    caption=caption,
                    esmi=cls.build_esmi(smi, groups, ext),
                    markush=True,
                    sru=is_sru,
                )

            mol = Chem.MolFromSmiles(
                Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
            )

            if is_markush or preserve_dummy_groups:
                remaining_groups = cls.remove_atom_groups(groups, consumed_atom_groups)
                new_groups = chem_utils.remap_groups(mol, remaining_groups, ring_info)
            else:
                new_groups = ""
                for atom in mol.GetAtoms():
                    atom.SetAtomMapNum(0)

            new_smi = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
            new_esmi = cls.build_esmi(new_smi, new_groups, ext)

            return TranslatedMolecule(
                smi=new_smi,
                groups=new_groups,
                caption=caption,
                esmi=new_esmi,
                markush=is_markush,
                sru=is_sru,
            )
        except Exception as e:
            if error_msg:
                logger.error(f"Error while refactoring SMILES: {smi}")
                logger.error(repr(e))
            return TranslatedMolecule(
                smi=smi,
                groups=groups,
                caption=caption,
                esmi=cls.build_esmi(smi, groups, ext),
                markush=len(groups) > 0,
                sru=False,
            )


__all__ = [
    "AtomIndex",
    "GroupDesc",
    "Index",
    "LONG_CHAR_ELEMENTS",
    "Patterns",
    "PERIODIC_TABLE",
    "RingIndex",
    "TextType",
    "Tokens",
    "TranslatedMolecule",
    "Translator",
]
