from __future__ import annotations

import unittest

from molparser.utils import postprocess_caption, substitute_markush
from molparser.utils.translator import Translator


class TranslatorPublicSubstitutionTests(unittest.TestCase):
    def test_translator_exposes_markush_substitution(self) -> None:
        self.assertEqual(
            Translator.substitute_markush(
                "*C<sep><a>0:R[1]</a>",
                {"R1": "Me"},
            ),
            "CC",
        )

    def test_existing_v1_return_shapes_stay_unchanged(self) -> None:
        atom_result = Translator.substitute_markush(
            "*C<sep><a>0:R[1]</a>",
            {"R1": "Me"},
        )
        ring_result = Translator.substitute_markush(
            "c1ccccc1<sep><r>0:R[1]?1-2</r>",
            {"R1": "Me"},
        )
        ch2_result = Translator.substitute_markush(
            "C*CO<sep><a>1:CH2?3</a>",
            {},
        )

        self.assertIsInstance(atom_result, str)
        self.assertEqual(atom_result, substitute_markush(
            "*C<sep><a>0:R[1]</a>", {"R1": "Me"}
        ))
        self.assertIsInstance(ring_result, list)
        self.assertEqual(
            ring_result,
            [
                "Cc1ccc(C)cc1",
                "Cc1cccc(C)c1",
                "Cc1ccccc1",
                "Cc1ccccc1C",
            ],
        )
        self.assertEqual(ch2_result, "CCCCCO")

    def test_existing_v1_abbreviation_substitution_is_unchanged(self) -> None:
        self.assertEqual(
            Translator.substitute_markush("*C<sep><a>0:CF3</a>", {}),
            "CC(F)(F)F",
        )


class ESmiles2CountSubstitutionTests(unittest.TestCase):
    def test_nested_substructure_label_and_count_substitute_recursively(self) -> None:
        caption = "C<sep><s>*C<sep><a>0:R[1]</a>|Sg:n|</s>"
        self.assertEqual(
            substitute_markush(caption, {"R1": "Me", "n": 3}),
            "C<sep><s>CC<sep>|Sg:3|</s>",
        )

    def test_nested_substructure_preserves_v2_dummy_endpoints(self) -> None:
        caption = (
            "C<sep><s>*C*<sep><d>0:<dum></d><d>2:<dum></d>"
            "|Sg:n|</s>"
        )
        self.assertEqual(
            substitute_markush(caption, {"n": 4}),
            "C<sep><s>*C*<sep><d>0:<dum></d><d>2:<dum></d>|Sg:4|</s>",
        )

    def test_sgroup_count_substitutes_without_graph_expansion(self) -> None:
        caption = "CCOCC<sep><g>[1:0]:[3:4]:|Sg:n|</g>"
        self.assertEqual(
            substitute_markush(caption, {"n": 3}),
            "CCOCC<sep><g>[1:0]:[3:4]:|Sg:3|</g>",
        )

    def test_nested_sgroup_count_substitutes_without_graph_expansion(self) -> None:
        caption = (
            "C<sep><s>CCOCC<sep>"
            "<g>[1:0]:[3:4]:|Sg:m|</g></s>"
        )
        self.assertEqual(
            substitute_markush(caption, {"m": "12"}),
            "C<sep><s>CCOCC<sep><g>[1:0]:[3:4]:|Sg:12|</g></s>",
        )

    def test_missing_sgroup_count_definition_preserves_symbol(self) -> None:
        caption = "CCOCC<sep><g>[1:0]:[3:4]:|Sg:n|</g>"
        self.assertEqual(substitute_markush(caption, {}), caption)

    def test_top_level_sru_count_preserves_v2_dummy_endpoints(self) -> None:
        caption = "*CC*<sep><d>0:<dum></d><d>3:<dum></d>|Sg:n|"
        self.assertEqual(
            substitute_markush(caption, {"n": 3}),
            "*CC*<sep><d>0:<dum></d><d>3:<dum></d>|Sg:3|",
        )

    def test_top_level_sru_count_preserves_v1_dummy_syntax(self) -> None:
        caption = "*CC*<sep><a>0:<dum></a><a>3:<dum></a>|Sg:n|"
        self.assertEqual(
            substitute_markush(caption, {"n": 3}),
            "*CC*<sep><a>0:<dum></a><a>3:<dum></a>|Sg:3|",
        )

    def test_multiple_nested_label_values_keep_list_return_type(self) -> None:
        caption = "C<sep><s>*C<sep><a>0:R[1]</a>|Sg:n|</s>"
        result = substitute_markush(
            caption,
            {"R1": ["Me", "Cl"], "n": 2},
        )
        self.assertIsInstance(result, list)
        self.assertEqual(
            set(result),
            {
                "C<sep><s>CC<sep>|Sg:2|</s>",
                "C<sep><s>CCl<sep>|Sg:2|</s>",
            },
        )

    def test_sgroup_count_definition_must_be_one_positive_integer(self) -> None:
        caption = "CCOCC<sep><g>[1:0]:[3:4]:|Sg:n|</g>"
        for invalid in (0, -1, "Me", [2, 3]):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    substitute_markush(caption, {"n": invalid})

    def test_explicit_and_range_sgroup_counts_are_not_definition_keys(self) -> None:
        explicit = "CCOCC<sep><g>[1:0]:[3:4]:|Sg:20|</g>"
        ranged = "CCOCC<sep><g>[1:0]:[3:4]:|Sg:2-8|</g>"
        self.assertEqual(substitute_markush(explicit, {"20": 7}), explicit)
        self.assertEqual(substitute_markush(ranged, {"2-8": 7}), ranged)

    def test_special_id_preservation_is_unchanged(self) -> None:
        caption = "*C<sep><a>0:<id>[DNA]</a>"
        self.assertEqual(substitute_markush(caption, {}), caption)


class PrecompatibleClassificationTests(unittest.TestCase):
    def test_precompatible_records_are_reported_as_markush(self) -> None:
        captions = (
            "C<sep><s>*C<sep><a>0:R[1]</a>|Sg:n|</s>",
            "CCOCC<sep><g>[1:0]:[3:4]:|Sg:n|</g>",
            "C=CC<sep><v>0:A:[0:2]</v>",
        )
        for caption in captions:
            with self.subTest(caption=caption):
                result = postprocess_caption(caption)
                self.assertTrue(result["markush"])
                self.assertEqual(result["groups"], caption.split("<sep>", 1)[1])

    def test_explicit_top_level_sru_count_is_recognized(self) -> None:
        caption = "*CC*<sep><d>0:<dum></d><d>3:<dum></d>|Sg:3|"
        self.assertTrue(postprocess_caption(caption)["sru"])

    def test_supported_symbol_and_range_sru_counts_are_recognized(self) -> None:
        for count in ("n", "m", "20", "2-8", "m-n"):
            with self.subTest(count=count):
                caption = (
                    "*CC*<sep><d>0:<dum></d><d>3:<dum></d>"
                    f"|Sg:{count}|"
                )
                self.assertTrue(postprocess_caption(caption)["sru"])

    def test_unusual_top_level_sru_count_is_not_recognized(self) -> None:
        caption = "*CC*<sep><d>0:<dum></d><d>3:<dum></d>|Sg:!|"
        self.assertFalse(postprocess_caption(caption)["sru"])


if __name__ == "__main__":
    unittest.main()
