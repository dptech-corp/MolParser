from __future__ import annotations

import unittest

from rdkit import Chem

from molparser import utils as mutils
from molparser.utils.translator import Translator


class ConcreteRepeatPostprocessTests(unittest.TestCase):
    def assert_parseable(self, value: str) -> None:
        self.assertIsNotNone(Chem.MolFromSmiles(value), value)

    def test_whole_sru_fixed_count_populates_all_formats(self) -> None:
        raw = "*CC*<sep><d>0:<dum></d><d>3:<dum></d>|Sg:3|"
        result = mutils.postprocess_caption(raw)

        self.assertEqual(result["caption"], raw)
        self.assertEqual(result["smi"], "*CCCCCC*")
        self.assertEqual(result["esmi"], "*CCCCCC*<sep>")
        self.assertTrue(result["cxsmiles"].startswith("*CCCCCC*"))
        self.assertNotIn("Sg:", result["cxsmiles"])
        self.assertEqual(result["groups"], "")
        self.assertFalse(result["markush"])
        self.assertTrue(result["sru"])
        self.assert_parseable(result["smi"])
        self.assert_parseable(result["cxsmiles"])

    def test_local_sgroup_fixed_count_populates_all_formats(self) -> None:
        raw = "CCOCC<sep><g>[1:0]:[3:4]:|Sg:3|</g>"
        result = mutils.postprocess_caption(raw)

        self.assertEqual(result["caption"], raw)
        self.assertEqual(result["smi"], "CCOCCOCCOCC")
        self.assertEqual(result["esmi"], "CCOCCOCCOCC<sep>")
        self.assertEqual(result["cxsmiles"], "CCOCCOCCOCC")
        self.assertEqual(result["groups"], "")
        self.assertFalse(result["markush"])
        self.assertFalse(result["sru"])
        self.assert_parseable(result["smi"])
        self.assert_parseable(result["cxsmiles"])

    def test_ch2_fixed_count_remaps_residual_atom_id_once(self) -> None:
        raw = "*CCO<sep><a>0:CH2?3</a><c>3:B</c>"
        result = mutils.postprocess_caption(raw)

        self.assertEqual(result["smi"], "CCCCCO")
        self.assertEqual(result["esmi"], "CCCCCO<sep><c>5:B</c>")
        self.assertEqual(result["groups"], "<c>5:B</c>")
        self.assertEqual(result["cxsmiles"], "CCCCCO |$;;;;;cB$|")
        self.assertNotIn(";;;cB;;cB", result["cxsmiles"])
        self.assert_parseable(result["smi"])
        self.assert_parseable(result["cxsmiles"])

        direct = Translator.esmiles_to_cxsmiles(raw)
        self.assertEqual(direct, "CCCCCO |$;;;;;cB$|")

    def test_ch2_fixed_count_remaps_special_id(self) -> None:
        raw = "*CCO<sep><a>0:CH2?3</a><a>3:<id>[OH]</a>"
        result = mutils.postprocess_caption(raw)

        self.assertEqual(result["esmi"], "CCCCCO<sep><a>5:<id>[OH]</a>")
        self.assertEqual(result["groups"], "<a>5:<id>[OH]</a>")
        self.assertEqual(result["cxsmiles"], "CCCCCO |$;;;;;<id>[OH]$|")
        self.assert_parseable(result["cxsmiles"])

    def test_ambiguous_fixed_abbreviation_stays_esmiles(self) -> None:
        raw = "*CC<sep><a>0:Me?3</a>"
        result = mutils.postprocess_caption(raw)

        self.assertEqual(result["smi"], "*CC")
        self.assertIn("<a>0:Me?3</a>", result["esmi"])
        self.assertIn("<a>0:Me?3</a>", result["groups"])
        self.assertTrue(result["markush"])

    def test_unsafe_ch2_fixed_count_stays_esmiles(self) -> None:
        raw = "C=*C<sep><a>1:CH2?3</a>"
        result = mutils.postprocess_caption(raw)

        self.assertEqual(result["smi"], "C=*C")
        self.assertIn("<a>1:CH2?3</a>", result["esmi"])
        self.assertTrue(result["markush"])

    def test_oversized_ch2_fixed_count_stays_esmiles(self) -> None:
        raw = "*CCO<sep><a>0:CH2?1025</a>"
        result = mutils.postprocess_caption(raw)

        self.assertEqual(result["smi"], "*CCO")
        self.assertIn("CH2?1025", result["esmi"])

    def test_ring_outputs_only_use_available_substitution_sites(self) -> None:
        cxsmiles = Translator.esmiles_to_cxsmiles(
            "CCCc1ccccc1<sep><r>0:Me</r>"
        )
        self.assert_parseable(cxsmiles)
        self.assertIn("m:", cxsmiles)

        values = mutils.substitute_markush(
            "*c1ccccc1*<sep><r>0:Me?2</r>",
            {},
        )
        self.assertIsInstance(values, list)
        for value in values:
            mol = Chem.MolFromSmiles(value)
            self.assertIsNotNone(mol, value)
            self.assertEqual(sum(atom.GetIsAromatic() for atom in mol.GetAtoms()), 6)

    def test_ring_substitution_and_whole_repeat_do_not_leak_invalid_states(self) -> None:
        values = mutils.substitute_markush(
            "*c1ccccc1*<sep><d>0:<dum></d><d>7:<dum></d>"
            "<r>0:Me</r>|Sg:2|",
            {},
        )
        values = [values] if isinstance(values, str) else values
        self.assertTrue(values)
        for value in values:
            self.assert_parseable(value)

    def test_direct_cxsmiles_preserves_unexpanded_sru_label(self) -> None:
        value = Translator.esmiles_to_cxsmiles("*CC*<sep>|Sg:3|")
        self.assertIn(":3:ht", value)
        self.assertNotIn(":n:ht", value)
        self.assert_parseable(value)

    def test_symbolic_multiplicity_rejects_boolean_definition(self) -> None:
        with self.assertRaisesRegex(ValueError, "one integer"):
            mutils.substitute_markush(
                "*CCO<sep><a>0:CH2?n</a>",
                {"n": True},
            )


if __name__ == "__main__":
    unittest.main()
