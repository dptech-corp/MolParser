from __future__ import annotations

import unittest

from rdkit import Chem

from molparser import utils as mutils


class MarkushRepeatRegressionTests(unittest.TestCase):
    def assert_plain_smiles(self, value: str) -> None:
        self.assertNotIn("<sep>", value)
        self.assertIsNotNone(Chem.MolFromSmiles(value), value)

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


if __name__ == "__main__":
    unittest.main()
