from __future__ import annotations

import hashlib
import json
import math
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
from unittest.mock import patch

import rdkit
from rdkit.Chem import Draw

from molparser.utils import drawer


class DrawerLegacyCompatibilityTests(unittest.TestCase):
    def test_rdkit_2025_09_6_legacy_svgs_are_byte_exact(self) -> None:
        fixture_path = (
            Path(__file__).parent
            / "fixtures"
            / "legacy_svg_rdkit_2025_09_6.json"
        )
        fixture = json.loads(fixture_path.read_text(encoding="utf-8"))
        if rdkit.__version__ != fixture["rdkit"]:
            self.skipTest(
                f"golden requires RDKit {fixture['rdkit']}, got {rdkit.__version__}"
            )

        for name, record in fixture["cases"].items():
            with self.subTest(case=name):
                svg = drawer.draw(record["caption"], output_format="svg")
                encoded = svg.encode("utf-8")
                self.assertEqual(len(encoded), record["bytes"])
                self.assertEqual(
                    hashlib.sha256(encoded).hexdigest(),
                    record["sha256"],
                )

    def test_font_none_preserves_default_bytes(self) -> None:
        caption = "CC(=O)NCCO<sep>"
        self.assertEqual(
            drawer.draw(caption),
            drawer.draw(caption, config={"visual": {"fontFile": None}}),
        )

    def test_explicit_font_is_deterministic_and_opt_in(self) -> None:
        font = Path(Draw.__file__).resolve().parent / "FreeSans.ttf"
        self.assertTrue(font.is_file())
        caption = "CC(=O)NCCO<sep>"
        configured = {"visual": {"fontFile": str(font)}}
        first = drawer.draw(caption, config=configured)
        self.assertEqual(first, drawer.draw(caption, config=configured))
        self.assertNotEqual(first, drawer.draw(caption))


class DrawerRepeatOverlayTests(unittest.TestCase):
    def test_new_repeat_syntax_auto_enables_overlays(self) -> None:
        cases = (
            "*CCOCC*<sep><d>0:<dum></d><d>6:<dum></d>|Sg:n|",
            "COCCNCCOC<sep><g>[4:3]:[5:6]:|Sg:n|</g>",
            "CC*CC<sep><a>2:CH2?3</a>",
        )
        for caption in cases:
            with self.subTest(caption=caption):
                automatic = drawer.draw(caption)
                disabled = drawer.draw(
                    caption,
                    config={"features": {"repeat_units": False}},
                )
                self.assertNotEqual(automatic, disabled)
                ET.fromstring(automatic)

    def test_single_ch2_is_an_unlabelled_vertex_with_two_incident_bonds(self) -> None:
        svg = drawer.draw("CC*CC<sep><a>2:CH2?3</a>")
        root = ET.fromstring(svg)
        paths = [
            node
            for node in root.iter()
            if node.tag.rsplit("}", 1)[-1] == "path"
        ]
        target_glyphs = [
            node
            for node in paths
            if "atom-2" in node.attrib.get("class", "")
            and "bond-" not in node.attrib.get("class", "")
        ]
        target_bonds = {
            node.attrib.get("class", "").split()[0]
            for node in paths
            if "atom-2" in node.attrib.get("class", "")
            and "bond-" in node.attrib.get("class", "")
        }
        self.assertEqual(target_glyphs, [])
        self.assertEqual(len(target_bonds), 2)

    def test_round_and_square_local_brackets_are_deterministic(self) -> None:
        caption = "COCCNCCOC<sep><g>[4:3]:[5:6]:|Sg:n|</g>"
        round_svg = drawer.draw(
            caption,
            config={"styling": {"sgroup_bracket_style": "round"}},
        )
        square_svg = drawer.draw(
            caption,
            config={"styling": {"sgroup_bracket_style": "square"}},
        )
        self.assertEqual(round_svg, drawer.draw(caption))
        self.assertNotEqual(round_svg, square_svg)
        ET.fromstring(square_svg)

    def test_paired_near_vertical_axes_snap_together(self) -> None:
        def axis(degrees: float) -> tuple[float, float]:
            radians = math.radians(degrees)
            return math.sin(radians), math.cos(radians)

        snapped = drawer._snap_paired_bracket_axes(
            [axis(44.9), axis(-44.9)]
        )
        self.assertEqual([value[0] for value in snapped], [0.0, 0.0])
        mixed = [axis(44.9), axis(45.1)]
        self.assertEqual(drawer._snap_paired_bracket_axes(mixed), mixed)


class DrawerVirtualArcTests(unittest.TestCase):
    old_named = (
        "C=CCC(C(C)*)*<sep><a>6:R[2]</a><a>7:R[1]</a>"
        "<v>0:A:[0:2]</v><r><v>0:R[3]</r>"
    )

    def test_old_single_named_arc_stays_legacy_unless_forced(self) -> None:
        automatic = drawer.draw(self.old_named)
        legacy = drawer.draw(
            self.old_named,
            config={"features": {"academic_virtual_arcs": False}},
        )
        academic = drawer.draw(
            self.old_named,
            config={"features": {"academic_virtual_arcs": True}},
        )
        self.assertEqual(automatic, legacy)
        self.assertNotEqual(automatic, academic)

    def test_empty_name_and_multiple_arcs_auto_use_academic_layout(self) -> None:
        captions = (
            "C=CCC(C(C)*)*<sep><a>6:R[2]</a><a>7:R[1]</a>"
            "<v>0::[0:2]</v>",
            "CCCCCC<sep><v>0:A:[0:3]</v><v>1:B:[0:5]</v>",
        )
        for caption in captions:
            with self.subTest(caption=caption):
                original = drawer._virtual_arc_layout
                with patch.object(
                    drawer,
                    "_virtual_arc_layout",
                    wraps=original,
                ) as layout:
                    automatic = drawer.draw(caption)
                self.assertGreater(layout.call_count, 0)
                with patch.object(
                    drawer,
                    "_virtual_arc_layout",
                    wraps=original,
                ) as disabled_layout:
                    drawer.draw(
                        caption,
                        config={"features": {"academic_virtual_arcs": False}},
                    )
                self.assertEqual(disabled_layout.call_count, 0)
                ET.fromstring(automatic)

    def test_empty_name_and_zero_refs_parse(self) -> None:
        self.assertEqual(
            drawer._virtual_arcs("<v>0::[0:2]</v>"),
            [(0, "", 0, 2)],
        )
        self.assertEqual(drawer._virtual_arc_substituents(""), {})


class DrawerEndpointBallTests(unittest.TestCase):
    def test_allow_list_tokens_auto_draw_flat_standard_colours(self) -> None:
        for token, (fill, _outline) in drawer._ENDPOINT_BALL_STYLES.items():
            with self.subTest(token=token):
                caption = f"CC*<sep><a>2:<id>[{token}]</a>"
                svg = drawer.draw(caption)
                fill_hex = "#" + "".join(
                    f"{int(channel * 255):02X}" for channel in fill
                )
                self.assertIn(f"fill:{fill_hex}", svg)
                self.assertNotEqual(
                    svg,
                    drawer.draw(
                        caption,
                        config={"features": {"endpoint_balls": False}},
                    ),
                )

    def test_unapproved_special_id_remains_text(self) -> None:
        caption = "CC*<sep><a>2:<id>[DNA]</a>"
        self.assertEqual(
            drawer.draw(caption),
            drawer.draw(
                caption,
                config={"features": {"endpoint_balls": False}},
            ),
        )


class DrawManyTests(unittest.TestCase):
    captions = ["CCO", "CCN", "CCC", "COC"]

    def test_reuses_one_validated_config_and_preserves_order(self) -> None:
        expected = [drawer.draw(caption) for caption in self.captions]
        original = drawer._validate_config
        with patch.object(drawer, "_validate_config", wraps=original) as validate:
            actual = drawer.draw_many(self.captions)
        self.assertEqual(validate.call_count, 1)
        self.assertEqual(actual, expected)

    def test_process_results_match_sequential_results(self) -> None:
        sequential = drawer.draw_many(self.captions, workers=1)
        self.assertEqual(
            drawer.draw_many(self.captions, workers=2),
            sequential,
        )

    def test_rejects_string_and_invalid_worker_count(self) -> None:
        with self.assertRaises(TypeError):
            drawer.draw_many("CCO")
        with self.assertRaises(ValueError):
            drawer.draw_many([], workers=0)


if __name__ == "__main__":
    unittest.main()
