from __future__ import annotations

import unittest

from molparser import utils as mutils


class PublicUtilsApiTests(unittest.TestCase):
    def test_batch_draw_and_translator_substitution_are_public(self) -> None:
        self.assertEqual(mutils.draw_many(["CCO"]), [mutils.draw("CCO")])
        self.assertEqual(
            mutils.Translator.substitute_markush(
                "CCOCC<sep><g>[1:0]:[3:4]:|Sg:n|</g>",
                {"n": 3},
            ),
            "CCOCC<sep><g>[1:0]:[3:4]:|Sg:3|</g>",
        )


if __name__ == "__main__":
    unittest.main()
