from __future__ import annotations

import unittest

from rdkit import Chem

from molparser import utils as mutils


class MarkushRepeatRegressionTests(unittest.TestCase):
    def assert_plain_smiles(self, value: str) -> None:
        self.assertNotIn("<sep>", value)
        self.assertIsNotNone(Chem.MolFromSmiles(value), value)

    def assert_same_isomer(self, actual: str, expected: str) -> None:
        actual_mol = Chem.MolFromSmiles(actual)
        expected_mol = Chem.MolFromSmiles(expected)
        self.assertIsNotNone(actual_mol, actual)
        self.assertIsNotNone(expected_mol, expected)
        self.assertEqual(
            Chem.MolToSmiles(actual_mol, canonical=True, isomericSmiles=True),
            Chem.MolToSmiles(expected_mol, canonical=True, isomericSmiles=True),
        )

    def test_whole_sru_expands(self) -> None:
        value = mutils.substitute_markush(
            "*CC*<sep><d>0:<dum></d><d>3:<dum></d>|Sg:3|", {}
        )
        self.assertEqual(value, "*CCCCCC*")
        self.assert_plain_smiles(value)

    def test_local_sgroup_expands(self) -> None:
        value = mutils.substitute_markush(
            "CCOCC<sep><g>[1:0]:[3:4]:|Sg:n|</g>", {"n": 3}
        )
        self.assertEqual(value, "CCOCCOCCOCC")
        self.assert_plain_smiles(value)

    def test_peg_twelve_hydrogen_caps(self) -> None:
        value = mutils.substitute_markush(
            "*OCCO*<sep><d>0:<dum></d><d>5:<dum></d>"
            "<g>[1:0]:[3:4]:|Sg:12|</g>",
            {},
            terminal_policy="hydrogen",
        )
        self.assertEqual(value, "O" + "CCO" * 12)
        mol = Chem.MolFromSmiles(value)
        self.assertIsNotNone(mol)
        self.assertEqual(mol.GetNumHeavyAtoms(), 37)

    def test_ch2_range_enumerates_valid_smiles(self) -> None:
        values = mutils.substitute_markush(
            "*CCO<sep><a>0:CH2?1-3</a>", {}
        )
        self.assertEqual(set(values), {"CCCO", "CCCCO", "CCCCCO"})
        for value in values:
            self.assert_plain_smiles(value)

    def test_ring_symbolic_count_uses_definition(self) -> None:
        values = mutils.substitute_markush(
            "c1ccccc1<sep><r>0:R[1]?n</r>", {"R1": "Me", "n": 3}
        )
        self.assertEqual(len(values), 3)
        for value in values:
            self.assert_plain_smiles(value)

    def test_best_effort_substitutes_and_preserves_unresolved(self) -> None:
        value = mutils.substitute_markush(
            "*CC*<sep><a>0:R[1]</a><a>3:R[2]</a>", {"R1": "Me"}
        )
        self.assertIn("<sep>", value)
        self.assertIn("R[2]", value)
        self.assertNotIn("R[1]", value)
        self.assertIsNotNone(Chem.MolFromSmiles(value.split("<sep>", 1)[0]))

    def test_nested_scope_resolves_count_but_stays_esmiles(self) -> None:
        value = mutils.substitute_markush(
            "C<sep><s>*C<sep><a>0:R[9]</a>|Sg:n|</s>", {"n": 3}
        )
        self.assertEqual(
            value, "C<sep><s>*C<sep><a>0:R[9]</a>|Sg:3|</s>"
        )

    def test_unsafe_ch2_target_is_not_silently_rewritten(self) -> None:
        value = mutils.substitute_markush(
            "C=*C<sep><a>1:CH2?2</a>", {}
        )
        self.assertIn("<sep>", value)

    def test_strict_rejects_residual_repeat(self) -> None:
        with self.assertRaises(ValueError):
            mutils.substitute_markush(
                "C<sep><s>*C<sep><a>0:R[9]</a>|Sg:3|</s>",
                {},
                repeat_policy="strict",
            )

    def test_max_outputs_is_a_hard_limit(self) -> None:
        with self.assertRaises(ValueError):
            mutils.substitute_markush(
                "*C<sep><a>0:R[1]</a>",
                {"R1": ["C", "CC", "CCC"]},
                max_outputs=2,
            )
        with self.assertRaises(ValueError):
            mutils.substitute_markush("CC", {}, max_outputs=0)

    def test_attachment_preserves_tetrahedral_configuration(self) -> None:
        target = mutils.substitute_markush(
            "F[C@H](*)Cl<sep><a>2:Me</a>", {}
        )
        source = mutils.substitute_markush(
            "*C<sep><a>0:R[1]</a>", {"R1": "F[C@H](*)Cl"}
        )
        self.assert_same_isomer(target, "C[C@@H](F)Cl")
        self.assert_same_isomer(source, "C[C@@H](F)Cl")

    def test_attachment_preserves_double_bond_stereo(self) -> None:
        target = mutils.substitute_markush(
            "*/C=C/C<sep><a>0:Me</a>", {}
        )
        source = mutils.substitute_markush(
            "*C<sep><a>0:R[1]</a>", {"R1": "*/C=C/C"}
        )
        self.assert_same_isomer(target, "C/C=C/C")
        self.assert_same_isomer(source, "C/C=C/C")

    def test_postprocess_repeat_and_abbreviation_preserve_ez(self) -> None:
        result = mutils.postprocess_caption(
            "*/C=C/C*<sep><a>0:Me</a><a>4:CH2?2</a>"
        )
        self.assert_same_isomer(result["smi"], "C/C=C/CCC")
        self.assertEqual(result["esmi"], f"{result['smi']}<sep>")
        self.assertEqual(result["cxsmiles"], result["smi"])

    def test_postprocess_drops_only_unrepairable_indices(self) -> None:
        atom = mutils.postprocess_caption("CC<sep><a>9:Me</a>")
        ring = mutils.postprocess_caption("c1ccccc1<sep><r>9:Me</r>")
        mixed = mutils.postprocess_caption(
            "CC<sep><a>9:Me</a><a>1:<id>[OH]</a>"
        )

        self.assertEqual(atom["esmi"], "CC<sep>")
        self.assertEqual(ring["esmi"], "c1ccccc1<sep>")
        self.assertFalse(atom["markush"])
        self.assertFalse(ring["markush"])
        self.assertNotIn("<a>9:", mixed["esmi"])
        self.assertIn("<id>[OH]", mixed["esmi"])
        self.assertTrue(mixed["markush"])


if __name__ == "__main__":
    unittest.main()
