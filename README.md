# MolParser

MolParser is a toolkit for working with **E-SMILES** (extended SMILES) in OCSR and Markush workflows. The notation follows the formulation introduced in the [MolParser paper](https://arxiv.org/abs/2411.11098).


| Path                                | Role                                                                               |
| ----------------------------------- | ---------------------------------------------------------------------------------- |
| `utils/`                            | MolParser utils (substitute known abbreviations, and render E-SMILES to structure) |
| `skills/molparser-extended-smiles/` | E-SMILES skills (concise rules and examples for LLM / OCSR agents)                 |


## Installation

```bash
pip install -r requirements.txt
```

Run the examples below from the **repository root** so that `from utils import ...` resolves correctly.

## E-SMILES overview

E-SMILES combines a base SMILES with an optional extension:

```text
SMILES<sep>EXTENSION
```

Common extension records:

- `<a>0:R[1]</a>` — atom-indexed substituent or Markush placeholder
- `<r>0:R[1]</r>` — ring-indexed substituent (regio-uncertain attachment)
- `<c>9:B</c>` — abstract-ring or superatom placeholder
- `<a>0:<dum></a>` — explicit dummy attachment point
- `|Sg:n|` — structural repeating unit (SRU) marker

Full specification: `[skills/molparser-extended-smiles/extended-smiles-spec.md](skills/molparser-extended-smiles/extended-smiles-spec.md)`

## Quick start

### Post-process E-SMILES

Normalize an E-SMILES string and substitute known abbreviations. In this example, `CF3` is read from `utils/abbrevs_example.csv`, attached to the dummy atom `*`, and folded into ordinary SMILES. No unresolved Markush group remains, so `markush` is `False`.

```python
from utils import postprocess_caption

raw = "*c1ccccc1<sep><a>0:CF3</a>"
result = postprocess_caption(raw)

# caption: original input
# smi: abbreviation-substituted RDKit SMILES
# esmi: abbreviation-substituted E-SMILES
# markush: whether unresolved Markush groups remain
# sru: whether an SRU marker was detected
# groups: unresolved E-SMILES extension records
for key in ("caption", "smi", "esmi", "markush", "sru", "groups"):
    print(f"{key}: {result[key]}")
```

Expected output:

```text
caption: *c1ccccc1<sep><a>0:CF3</a>
smi: FC(F)(F)c1ccccc1
esmi: FC(F)(F)c1ccccc1<sep>
markush: False
sru: False
groups:
```

### Render E-SMILES

Render the E-SMILES as SVG and save it locally:

```python
from pathlib import Path
from utils import draw

svg_text = draw("*C(O)c1cc(C(=O)N(*)*)cc(-c2*ccc*2)c1<sep><a>0:CF3</a><a>9:R[3]</a><a>10:R[2]</a><a>14:X</a><a>18:Y</a><r>1:R[1]?1-3</r>", output_format="svg")

svg_path = Path("molecule.svg")
svg_path.write_text(svg_text, encoding="utf-8")
```

`molecule.svg` is a local render artifact.

To obtain a PNG from that SVG (requires `cairosvg` from `requirements.txt`):

```python
import cairosvg

png_path = Path("molecule.png")
cairosvg.svg2png(url=str(svg_path), write_to=str(png_path))
```

## LLM / OCSR workflow

### Skill context

Load these files for the agent:

- `skills/molparser-extended-smiles/SKILL.md`
- `skills/molparser-extended-smiles/extended-smiles-spec.md`
- `skills/molparser-extended-smiles/figure-index.md`

### Expected model output

```text
1. Base SMILES
2. E-SMILES in SMILES<sep>EXTENSION format
3. Markush status
4. Unsupported or ambiguous chemistry
```

### Validate and normalize

```bash
python skills/molparser-extended-smiles/validate_esmiles.py "<your_esmiles>"
```

Then normalize and render with `postprocess_caption` and `draw`.

## Related resources

- [Uni-Parser](https://arxiv.org/abs/2512.15098) — agent-oriented scientific document parsing with the latest MolParser. [Demo](https://uniparser.dp.tech/)
- [MolParser](https://arxiv.org/abs/2411.11098) — end-to-end molecular recognition. [Demo](https://ocsr.dp.tech/)
- [MolDetv2 weights](https://huggingface.co/UniParser/MolDetv2) — lightweight molecule detector. [Demo](https://huggingface.co/spaces/AI4Industry/MolDet)

