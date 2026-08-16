from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

from rdkit import Chem

from molparser import utils as mutils


VALIDATOR_PATH = (
    Path(__file__).parents[1]
    / "skills"
    / "molparser-extended-smiles"
    / "validate_esmiles.py"
)
SPEC = importlib.util.spec_from_file_location("skill_validate_esmiles", VALIDATOR_PATH)
assert SPEC is not None and SPEC.loader is not None
VALIDATOR = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = VALIDATOR
SPEC.loader.exec_module(VALIDATOR)


def errors(caption: str) -> list[str]:
    return [
        message.text
        for message in VALIDATOR.validate(caption, strict=True)
        if message.level == "error"
    ]


class SpecialIdContractTests(unittest.TestCase):
    def test_text_and_endpoint_ball_values_are_valid_atom_labels(self) -> None:
        self.assertEqual(errors("*C<sep><a>0:<id>[DNA]</a>"), [])
        self.assertEqual(errors("*C<sep><a>0:<id>[blue]</a>"), [])

    def test_special_id_is_not_a_ring_or_virtual_arc_payload(self) -> None:
        captions = [
            "c1ccccc1<sep><r>0:<id>[DNA]</r>",
            "*C<sep><c>0:<id>[DNA]</c>",
            "CCC<sep><v>0:A:[0:2]</v><r><v>0:<id>[DNA]</r>",
            "CCC<sep><v>0:<id>[DNA]:[0:2]</v>",
        ]
        for caption in captions:
            with self.subTest(caption=caption):
                self.assertTrue(any("atom-indexed <a>" in item for item in errors(caption)))

    def test_special_id_requires_a_nonempty_whitespace_free_note(self) -> None:
        self.assertTrue(errors("*C<sep><a>0:<id>[]</a>"))
        self.assertTrue(errors("*C<sep><a>0:<id>[two words]</a>"))
        self.assertTrue(errors("*C<sep><a>0:<id>[DNA]?3</a>"))
        self.assertTrue(errors("*C<sep><a>0:X<id>[DNA]</a>"))
        self.assertTrue(errors("c1ccccc1<sep><r>0:R<id>[DNA]</r>"))

    def test_substitution_preserves_special_id(self) -> None:
        value = mutils.substitute_markush(
            "*CC*<sep><a>0:R[1]</a><a>3:<id>[DNA]</a>",
            {"R1": "Me", "DNA": "Cl", "<id>[DNA]": "Br"},
        )
        self.assertIn("<sep>", value)
        self.assertIn("<id>[DNA]", value)
        self.assertNotIn("R[1]", value)
        self.assertIsNotNone(Chem.MolFromSmiles(value.split("<sep>", 1)[0]))

    def test_only_approved_values_on_dummy_atoms_render_as_balls(self) -> None:
        tokens = (
            "ball", "grey", "black", "green", "blue", "yellow", "purple",
            "orange", "pink", "brown",
        )
        for token in tokens:
            with self.subTest(token=token):
                svg = mutils.draw(f"*C<sep><a>0:<id>[{token}]</a>")
                self.assertIn("<ellipse", svg)

        plain_star = mutils.draw("*C<sep>")
        for note in ("DNA", "red", "gray", "Blue"):
            with self.subTest(note=note):
                svg = mutils.draw(f"*C<sep><a>0:<id>[{note}]</a>")
                self.assertNotIn("<ellipse", svg)
                self.assertNotEqual(svg, plain_star)

        disabled_svg = mutils.draw(
            "*C<sep><a>0:<id>[blue]</a>",
            config={"features": {"endpoint_balls": False}},
        )
        self.assertNotIn("<ellipse", disabled_svg)
        self.assertNotEqual(disabled_svg, plain_star)

        for base in ("CC", "N", "P", "S", "Cl"):
            with self.subTest(base=base):
                caption = f"{base}<sep><a>0:<id>[blue]</a>"
                svg = mutils.draw(caption)
                self.assertNotIn("<ellipse", svg)
                self.assertNotEqual(svg, mutils.draw(f"{base}<sep>"))


if __name__ == "__main__":
    unittest.main()
