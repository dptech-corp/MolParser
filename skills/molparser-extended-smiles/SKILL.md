---
name: molparser-extended-smiles
description: "Use for MolParser E-SMILES generation, validation, normalization, substituent substitution, rendering, and repair in OCSR/Markush workflows, including atom-indexed substituents, regio-uncertain ring attachments, abstract-ring superatoms, dummy attachment points, and SRU repeat markers."
---

# MolParser E-SMILES Skill

Use this skill when reading, writing, validating, normalizing, or rendering MolParser E-SMILES for OCSR and Markush tasks.

## Fast Workflow

1. Build an RDKit-compatible base SMILES from the molecular graph.
2. Emit `SMILES<sep>EXTENSION`; keep `<sep>` even when `EXTENSION` is empty.
3. Add extension records only when the base SMILES cannot carry the annotation:
   - `<a>[ATOM_INDEX]:[GROUP_LABEL]</a>`: atom-indexed substituent or abbreviation.
   - `<r>[RING_INDEX]:[GROUP_LABEL]</r>`: ring-indexed substituent with unspecified attachment atom.
   - `<c>[ATOM_INDEX]:[RING_LABEL]</c>`: abstract-ring or superatom placeholder at a dummy atom.
   - `<a>[ATOM_INDEX]:<dum></a>`: explicit dummy attachment point.
   - `?n`, `?1-3`, `?3`: multiplicity suffix on a group label.
   - `|Sg:n|`: structural repeating unit (SRU) repeat marker.
4. Keep indexes zero-based. Atom indexes and ring indexes are separate namespaces.
5. Return the E-SMILES first; add only concise notes for unsupported chemistry or ambiguity.

## Normalization And Substitution

Use `postprocess_caption` to normalize E-SMILES, substitute known abbreviations,
and get CXSMILES:

```python
from utils import postprocess_caption

result = postprocess_caption(raw_esmiles)

# result["smi"]: abbreviation-substituted RDKit SMILES
# result["esmi"]: normalized E-SMILES
# result["cxsmiles"]: CXSMILES converted from normalized E-SMILES
# result["markush"]: whether unresolved Markush groups remain
# result["sru"]: whether an SRU marker was detected
# result["groups"]: unresolved E-SMILES extension records
```

- Use `result["esmi"]` for normalized E-SMILES output.
- Use `result["cxsmiles"]` when CXSMILES output is needed.

## Rendering

Use rendering for visual QA, not as chemical validation:

```python
from utils import draw

svg_text = draw(result["esmi"], output_format="svg")
png_bytes = draw(result["esmi"], output_format="png")
```

The drawer displays atom substituents, dummy attachment points, abstract rings, and ring-level annotations from the E-SMILES extension.

## Validation Priorities

- Exactly one top-level `<sep>`.
- `<a>` indexes atoms; `<r>` indexes rings; `<c>` indexes the dummy atom carrying the abstract-ring label.
- `GROUP_LABEL` may be a common abbreviation (`Me`, `OMe`, `CF3`), a Markush label (`R[1]`), or `<dum>`.
- Use group-level multiplicity suffixes (`?n`, `?1-3`, `?3`) separately from SRU-level `|Sg:n|`.
- After canonicalization or abbreviation substitution, regenerate all affected extension indexes.

## Boundary Policy

- Do not invent tokens for unsupported coordination bonds, electron-transfer arrows, uncertain bond styles, or uncertain chirality.
- Preserve the encodable molecular backbone and report unencoded chemistry explicitly.

## Reference Order

1. `extended-smiles-spec.md`
2. `figure-index.md`
3. `validate_esmiles.py`
4. `source-provenance.md`
