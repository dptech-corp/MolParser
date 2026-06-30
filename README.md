# MolParser

MolParser toolkit for working with **E-SMILES** (extended SMILES) in OCSR and Markush workflows. The notation follows the formulation introduced in the [MolParser paper](https://arxiv.org/abs/2411.11098).


| Path                                | Role                                                                               |
| ----------------------------------- | ---------------------------------------------------------------------------------- |
| `molparser/utils/`                  | MolParser utils (normalize E-SMILES, substitute abbreviations, convert to CXSMILES, and render structures) |
| `skills/molparser-extended-smiles/` | E-SMILES skills (concise rules and examples for LLM / OCSR agents)                 |


## Installation

```bash
pip install -e .
```

After installation, use the utilities through the `molparser` package namespace:

```python
from molparser import utils as mutils
```

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
- `?n` — local substructure multiplicity suffixes

Full specification: [skills/molparser-extended-smiles/extended-smiles-spec.md](skills/molparser-extended-smiles/extended-smiles-spec.md)

## Quick start

### Post-process E-SMILES

Convert E-SMILES to SMILES with abbreviation substitution, and convert E-SMILES to CXSMILES on a best-effort basis.

```python
from molparser import utils as mutils

raw = "*c1ccccc1<sep><a>0:CF3</a>"
result = mutils.postprocess_caption(raw)

# caption: original input string
# smi: normalized RDKit SMILES after substituting known abbreviations
# esmi: normalized E-SMILES after substitution and index repair
# cxsmiles: CXSMILES generated from the normalized E-SMILES
# markush: True if unresolved Markush labels remain
# sru: True if a structural repeating unit marker was detected
# groups: unresolved E-SMILES extension records kept after normalization
for key in ("caption", "smi", "esmi", "cxsmiles", "markush", "sru", "groups"):
    print(f"{key}: {result[key]}")
```

Expected output:

```text
caption: *c1ccccc1<sep><a>0:CF3</a>
smi: FC(F)(F)c1ccccc1
esmi: FC(F)(F)c1ccccc1<sep>
cxsmiles: FC(F)(F)c1ccccc1
markush: False
sru: False
groups:
```

### Substitute Markush Definitions

Substitute Markush labels with a definition dictionary. Definition keys can use
either `R1` or `R[1]`; values can be known abbreviations or SMILES fragments.
When a fragment contains `*`, that atom is treated as the attachment point.

```python
from molparser import utils as mutils

result = mutils.substitute_markush(
    "*c1ccccc1<sep><a>0:R[1]</a>",
    {"R1": "Me", "R2": "*CCO"},
)
print(result)
# Cc1ccccc1
```

Ring-indexed Markush records expand regio-uncertain attachments into a SMILES
list. Multiplicity suffixes such as `?3`, `?1-3`, and `?n` encode local
substructure replication over possible ring sites; `?n` reads the replication
count from the definition dictionary.

```python
result = mutils.substitute_markush(
    "c1ccccc1<sep><r>0:R[1]?1-3</r>",
    {"R1": "Me"},
)
print(result)

# ['Cc1cc(C)cc(C)c1', 'Cc1ccc(C)c(C)c1', 'Cc1ccc(C)cc1', 'Cc1cccc(C)c1', 'Cc1cccc(C)c1C', 'Cc1ccccc1', 'Cc1ccccc1C']
```

### Render E-SMILES

Render the E-SMILES as SVG and save it locally:

```python
from pathlib import Path
from molparser import utils as mutils

svg_text = mutils.draw("*C(O)c1cc(C(=O)N(*)*)cc(-c2*ccc*2)c1<sep><a>0:CF3</a><a>9:R[3]</a><a>10:R[2]</a><a>14:X</a><a>18:Y</a><r>1:R[1]?1-3</r>", output_format="svg")

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

Then normalize and render with `mutils.postprocess_caption` and `mutils.draw`.

## Related resources

- [Uni-Parser](https://arxiv.org/abs/2512.15098) — agent-oriented scientific document parsing with the latest MolParser. [Demo](https://uniparser.dp.tech/)
- [MolParser](https://arxiv.org/abs/2411.11098) — end-to-end molecular recognition. [Demo](https://ocsr.dp.tech/)
- [MolDetv2 weights](https://huggingface.co/UniParser/MolDetv2) — lightweight molecule detector. [Demo](https://huggingface.co/spaces/AI4Industry/MolDet)

## 📖 Citation

```bibtex
@inproceedings{fang2025molparser,
  title={Molparser: End-to-end visual recognition of molecule structures in the wild},
  author={Fang, Xi and Wang, Jiankun and Cai, Xiaochen and Chen, Shangqian and Yang, Shuwen and Tao, Haoyi and Wang, Nan and Yao, Lin and Zhang, Linfeng and Ke, Guolin},
  booktitle={Proceedings of the IEEE/CVF International Conference on Computer Vision},
  pages={24528--24538},
  year={2025}
}
```

```bibtex
@article{fang2025uniparser,
  title={Uni-Parser Technical Report},
  author={Fang, Xi and Tao, Haoyi and Yang, Shuwen and Zhong, Suyang and Lu, Haocheng and Lyu, Han and Huang, Chaozheng and Li, Xinyu and Zhang, Linfeng and Ke, Guolin},
  journal={arXiv preprint arXiv:2512.15098},
  year={2025}
}
```

## License

The project code is licensed under the Apache License 2.0. See [LICENSE](LICENSE) for details. Apache 2.0 permits **commercial use**, modification, and distribution, **provided that the license and copyright notices are retained**.

Model weights, datasets, and third-party dependencies used with this project are subject to their respective licenses. Please review and comply with those licenses when using them.
