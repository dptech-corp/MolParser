from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SKILL_ROOT = REPO_ROOT / "skills" / "molparser-extended-smiles"
MANIFEST_PATH = SKILL_ROOT / "references" / "synthetic-examples.json"
REFERENCE_PATH = SKILL_ROOT / "references" / "synthetic-examples.md"
SCRIPT_PATH = SKILL_ROOT / "scripts" / "generate_synthetic_examples.py"
VALIDATOR_PATH = SKILL_ROOT / "validate_esmiles.py"
SKILL_PATH = SKILL_ROOT / "SKILL.md"
AGENT_PATH = SKILL_ROOT / "agents" / "openai.yaml"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_generator():
    spec = importlib.util.spec_from_file_location(
        "molparser_synthetic_example_generator", SCRIPT_PATH
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


GENERATOR = _load_generator()


def _load_validator():
    spec = importlib.util.spec_from_file_location(
        "molparser_synthetic_example_validator", VALIDATOR_PATH
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


VALIDATOR = _load_validator()


class SyntheticSkillExampleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
        cls.examples = cls.manifest["examples"]

    def test_manifest_contract_and_coverage(self) -> None:
        report = GENERATOR.validate_manifest(self.manifest)
        self.assertTrue(report["all_valid"])
        self.assertEqual(report["fixture_count"], 10)
        self.assertEqual(
            report["category_counts"],
            {
                "endpoint_ball": 2,
                "sgroup": 3,
                "virtual_arc": 4,
                "whole_sru": 1,
            },
        )
        self.assertEqual(
            report["sgroup_subtype_counts"],
            {"complex_local": 1, "simple_short": 1, "single_ch2": 1},
        )
        self.assertEqual(
            report["arc_appearances"],
            [
                "named_no_refs",
                "named_with_refs",
                "unnamed_no_refs",
                "unnamed_with_refs",
            ],
        )

    def test_reference_contains_every_fixture_verbatim(self) -> None:
        reference = REFERENCE_PATH.read_text(encoding="utf-8")
        for example in self.examples:
            with self.subTest(fixture=example["id"]):
                self.assertIn(example["esmiles"], reference)
        self.assertIn("Raw drawing versus substitution", reference)
        self.assertIn("E-SMILES 1.0 and 2.0 compatibility", reference)

    def test_every_fixture_passes_the_skill_validator_in_strict_mode(self) -> None:
        for example in self.examples:
            with self.subTest(fixture=example["id"]):
                messages = VALIDATOR.validate(example["esmiles"], strict=True)
                errors = [message.text for message in messages if message.level == "error"]
                self.assertEqual(errors, [])

    def test_pending_or_frozen_asset_state_is_self_consistent(self) -> None:
        status = self.manifest["status"]
        self.assertIn(status, {"pending-renderer-freeze", "frozen"})
        for example in self.examples:
            asset = SKILL_ROOT / example["asset"]["path"]
            expected_sha = example["asset"]["sha256"]
            with self.subTest(fixture=example["id"]):
                if status == "pending-renderer-freeze":
                    self.assertIsNone(expected_sha)
                    self.assertFalse(asset.exists())
                else:
                    self.assertRegex(expected_sha, r"^[0-9a-f]{64}$")
                    self.assertTrue(asset.is_file())
                    self.assertEqual(_sha256(asset), expected_sha)

    def test_check_only_cli_is_read_only_and_machine_readable(self) -> None:
        before = sorted(
            (path.relative_to(SKILL_ROOT).as_posix(), path.stat().st_size)
            for path in SKILL_ROOT.rglob("*")
            if path.is_file()
        )
        completed = subprocess.run(
            [sys.executable, "-B", str(SCRIPT_PATH), "--check-only"],
            cwd=REPO_ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr or completed.stdout)
        report = json.loads(completed.stdout)
        self.assertTrue(report["all_valid"])
        self.assertEqual(report["mode"], "check-only")
        self.assertEqual(report["fixture_count"], 10)
        self.assertEqual(report["asset_state"], self.manifest["status"])
        expected_verified = 0 if self.manifest["status"] == "pending-renderer-freeze" else 10
        self.assertEqual(report["verified_asset_count"], expected_verified)
        after = sorted(
            (path.relative_to(SKILL_ROOT).as_posix(), path.stat().st_size)
            for path in SKILL_ROOT.rglob("*")
            if path.is_file()
        )
        self.assertEqual(after, before)

    def test_render_rejects_an_unfrozen_drawer_sha_without_writing(self) -> None:
        asset_root = SKILL_ROOT / "assets" / "synthetic-examples"
        before = {
            path.name: _sha256(path)
            for path in asset_root.glob("*.svg")
        }
        self.assertEqual(len(before), 10)
        completed = subprocess.run(
            [
                sys.executable,
                "-B",
                str(SCRIPT_PATH),
                "--render",
                "--renderer-sha256",
                "0" * 64,
            ],
            cwd=REPO_ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 1)
        report = json.loads(completed.stdout)
        self.assertFalse(report["all_valid"])
        self.assertIn("drawer SHA mismatch", report["error"])
        after = {
            path.name: _sha256(path)
            for path in asset_root.glob("*.svg")
        }
        self.assertEqual(after, before)

    def test_skill_links_reference_and_agent_metadata_stays_aligned(self) -> None:
        skill_text = SKILL_PATH.read_text(encoding="utf-8")
        agent_text = AGENT_PATH.read_text(encoding="utf-8")
        self.assertIn("references/synthetic-examples.md", skill_text)
        self.assertIn('display_name: "MolParser E-SMILES"', agent_text)
        self.assertIn("$molparser-extended-smiles", agent_text)


if __name__ == "__main__":
    unittest.main()
