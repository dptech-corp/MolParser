from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


VALIDATOR_PATH = (
    Path(__file__).parents[1]
    / "skills"
    / "molparser-extended-smiles"
    / "validate_esmiles.py"
)
SPEC = importlib.util.spec_from_file_location("molparser_skill_validator", VALIDATOR_PATH)
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


def warnings(caption: str) -> list[str]:
    return [
        message.text
        for message in VALIDATOR.validate(caption, strict=True)
        if message.level == "warning"
    ]


class ESmilesVersionValidatorTests(unittest.TestCase):
    def test_esmiles_1_atom_ring_and_legacy_dummy_remain_valid(self) -> None:
        captions = (
            "*C<sep><a>0:R[1]</a>",
            "c1ccccc1<sep><r>0:R[1]</r>",
            "*CC*<sep><a>0:<dum></a><a>3:<dum></a>|Sg:n|",
        )
        for caption in captions:
            with self.subTest(caption=caption):
                self.assertEqual(errors(caption), [])

    def test_esmiles_2_named_and_unnamed_arcs_are_valid(self) -> None:
        captions = (
            "C=CC<sep><v>0:A:[0:2]</v>",
            "C=CC<sep><v>0::[0:2]</v>",
            "C=CC<sep><v>0::[0:2]</v><r><v>0:R[1]</r>",
        )
        for caption in captions:
            with self.subTest(caption=caption):
                self.assertEqual(errors(caption), [])

    def test_multiple_arcs_and_whitespace_tolerant_refs_are_valid(self) -> None:
        caption = (
            "C=CCC<sep><v>0::[0:2]</v><v>1:B:[1:3]</v>"
            "<r> <v>0:R[1]</r><r><v>1:R[2]</r>"
        )
        self.assertEqual(errors(caption), [])
        self.assertEqual(
            errors("CCC<sep><v>0:A:[0:2]</v><r><v>0:<id>[DNA]</r>"),
            [],
        )

    def test_reversed_arc_is_accepted_with_canonicalization_warning(self) -> None:
        messages = VALIDATOR.validate("C=CC<sep><v>0:A:[2:0]</v>", strict=True)
        self.assertFalse(any(message.level == "error" for message in messages))
        self.assertTrue(any("legacy-compatible" in message.text for message in messages))

    def test_arc_ids_refs_and_endpoint_bounds_are_checked(self) -> None:
        self.assertTrue(errors("C=CC<sep><v>1:A:[0:2]</v><r><v>0:R[1]</r>"))
        self.assertTrue(errors("C=CC<sep><v>0:A:[0:3]</v>"))
        self.assertTrue(errors("C=CC<sep><v>0:A:[1:1]</v>"))
        self.assertTrue(errors("CCC<sep><v>0:A:[0:2]</v><r><v>9:<id>[DNA]</r>"))
        self.assertTrue(errors("CCC<sep><r><v>9:<id>[DNA]</r>"))

    def test_arc_identity_and_pair_uniqueness_are_checked(self) -> None:
        self.assertTrue(
            errors("CCCC<sep><v>0:A:[0:2]</v><v>0:B:[1:3]</v>")
        )
        self.assertTrue(
            errors("CCCC<sep><v>0:A:[0:2]</v><v>1:B:[2:0]</v>")
        )
        warning_only = "CCCC<sep><v>1:A:[0:2]</v>"
        self.assertEqual(errors(warning_only), [])
        self.assertTrue(any("consecutive ids" in item for item in warnings(warning_only)))

    def test_sgroup_ports_are_checked_against_base_graph(self) -> None:
        self.assertEqual(
            errors("CCOCC<sep><g>[1:0]:[3:4]:|Sg:n|</g>"),
            [],
        )
        self.assertTrue(errors("CCOCC<sep><g>[1:0]:[3:9]:|Sg:n|</g>"))

    def test_sgroup_allows_one_or_more_unique_nonself_ports(self) -> None:
        self.assertEqual(errors("CCO<sep><g>[1:0]:|Sg:n|</g>"), [])
        self.assertEqual(
            errors("CCOCC<sep><g>[1:0]:[2:3]:[3:4]:|Sg:n|</g>"),
            [],
        )
        self.assertTrue(errors("CCO<sep><g>[1:1]:|Sg:n|</g>"))
        self.assertTrue(errors("CCO<sep><g>[1:0]:[1:0]:|Sg:n|</g>"))
        self.assertTrue(errors("CCOCC<sep><g>[1:0][3:4]|Sg:n|</g>"))

    def test_atom_ring_and_dummy_record_indices_are_checked(self) -> None:
        captions = (
            "CC<sep><a>2:R[1]</a>",
            "*C<sep><d>2:<dum></d>",
            "CC<sep><c>2:A</c>",
            "C1CC1<sep><r>1:R[1]</r>",
        )
        for caption in captions:
            with self.subTest(caption=caption):
                self.assertTrue(errors(caption))

    def test_nested_substructure_uses_its_own_atom_and_ring_namespaces(self) -> None:
        captions = (
            "C<sep><s>C=CC<sep><v>0:A:[0:2]</v></s>",
            "C<sep><s>CCOCC<sep><g>[1:0]:[3:4]:|Sg:n|</g></s>",
            "C<sep><s>C1CC1<sep><r>0:R[1]</r></s>",
        )
        for caption in captions:
            with self.subTest(caption=caption):
                self.assertEqual(errors(caption), [])

    def test_substructure_records_cannot_nest_another_substructure_record(self) -> None:
        caption = "C<sep><s>C<sep><s>CC<sep>|Sg:n|</s></s>"
        self.assertTrue(errors(caption))

    def test_top_level_sgroup_counts_are_not_double_validated_or_stale(self) -> None:
        for count in ("n", "m", "20", "2-8", "m-n"):
            with self.subTest(count=count):
                caption = f"*CC*<sep><d>0:<dum></d><d>3:<dum></d>|Sg:{count}|"
                self.assertEqual(errors(caption), [])
                self.assertFalse(
                    any("only flags" in item or "can be resolved" in item for item in warnings(caption))
                )
        self.assertTrue(errors("*CC*<sep>|Sg:n||Sg:m|"))
        self.assertTrue(errors("*CC*<sep>|Sg:0|"))
        self.assertTrue(errors("*CC*<sep>|Sg:4-2|"))


if __name__ == "__main__":
    unittest.main()
