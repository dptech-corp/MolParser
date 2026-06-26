---
name: molparser-extended-smiles
description: "Use for MolParser E-SMILES generation, validation, normalization, abbreviation substitution, Markush definition substitution, rendering, and repair in OCSR/Markush workflows, including atom-indexed substituents, regio-uncertain ring attachments, abstract-ring superatoms, dummy attachment points, local substructure multiplicity suffixes, and SRU repeat markers."
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
   - `?n`, `?1-3`, `?3`: local substructure multiplicity suffix on a group label.
   - `|Sg:n|`: structural repeating unit (SRU) repeat marker.
4. Keep indexes zero-based. Atom indexes and ring indexes are separate namespaces.
5. Return the E-SMILES first; add only concise notes for unsupported chemistry or ambiguity.

## Normalization And Substitution

Use `postprocess_caption` to convert E-SMILES to SMILES with abbreviation
substitution and get best-effort CXSMILES:

```python
from molparser import utils as mutils

result = mutils.postprocess_caption(raw_esmiles)

# result["smi"]: abbreviation-substituted RDKit SMILES
# result["esmi"]: normalized E-SMILES
# result["cxsmiles"]: CXSMILES converted from normalized E-SMILES
# result["markush"]: whether unresolved Markush groups remain
# result["sru"]: whether an SRU marker was detected
# result["groups"]: unresolved E-SMILES extension records
```

- Use `result["esmi"]` for normalized E-SMILES output.
- Use `result["cxsmiles"]` when CXSMILES output is needed.

## Markush Definition Substitution

Use `substitute_markush` when unresolved Markush labels should be expanded into
concrete SMILES. Definition keys can use `R1` or `R[1]`; values can be known
abbreviations (`Me`, `Cl`, `Ph`) or SMILES fragments. If a fragment contains
`*`, that atom is used as the attachment point. The return value is a single
SMILES string when the substitution is fully determined, or a list of SMILES
when regio-uncertain attachments or multiplicity ranges produce multiple
deduplicated structures.

```python
from molparser import utils as mutils

result = mutils.substitute_markush(
    "*c1ccccc1<sep><a>0:R[1]</a>",
    {"R1": "Me"},
)
# "Cc1ccccc1"
```

- Ring-indexed records such as `<r>0:R[1]</r>` enumerate regio-uncertain
  attachments and return a SMILES list after canonical de-duplication:
  `mutils.substitute_markush("c1ccccc1<sep><r>0:R[1]</r>", {"R1": "Me"})`
  returns `"Cc1ccccc1"` after symmetry de-duplication, while
  `<r>0:R[1]?1-3</r>` returns a list for 1-3 methyl substitutions on benzene.
- Multiplicity suffixes `?3`, `?1-3`, and `?n` encode local substructure
  replication; `?n` reads the replication count from the definition dictionary, e.g.
  `<a>2:CH2?n</a>` with `{"n": 10}`.

## Rendering

Use rendering for visual QA, not as chemical validation:

```python
from molparser import utils as mutils

svg_text = mutils.draw(result["esmi"], output_format="svg")
png_bytes = mutils.draw(result["esmi"], output_format="png")
```

The drawer displays atom substituents, dummy attachment points, abstract rings, and ring-level annotations from the E-SMILES extension.

## Validation Priorities

- Exactly one top-level `<sep>`.
- `<a>` indexes atoms; `<r>` indexes rings; `<c>` indexes the dummy atom carrying the abstract-ring label.
- `GROUP_LABEL` may be a common abbreviation (`Me`, `OMe`, `CF3`), a Markush label (`R[1]`), or `<dum>`.
- Use local substructure multiplicity suffixes (`?n`, `?1-3`, `?3`) separately from SRU-level `|Sg:n|`.
- After canonicalization or abbreviation substitution, regenerate all affected extension indexes.

## Boundary Policy

- Do not invent tokens for unsupported coordination bonds, electron-transfer arrows, uncertain bond styles, or uncertain chirality.
- Preserve the encodable molecular backbone and report unencoded chemistry explicitly.

## Reference Order

1. `extended-smiles-spec.md`
2. `figure-index.md`
3. `validate_esmiles.py`
4. `source-provenance.md`
