---
name: molparser-extended-smiles
description: "Use for MolParser E-SMILES generation, validation, normalization, abbreviation substitution, Markush definition substitution, rendering, and repair in OCSR/Markush workflows, including atom-indexed substituents, regio-uncertain ring attachments, abstract-ring superatoms, dummy attachment points, local substructure multiplicity suffixes, virtual arcs, SRU repeat markers, and MolParser colored endpoint-ball labels."
---

# MolParser E-SMILES Skill

Use this skill when reading, writing, validating, normalizing, or rendering MolParser E-SMILES for OCSR and Markush tasks.

## Fast Workflow

1. Build an RDKit-compatible base SMILES from the molecular graph.
2. Emit `SMILES<sep>EXTENSION`; keep `<sep>` even when `EXTENSION` is empty.
3. Add extension records only when the base SMILES cannot carry the annotation:
   - `<a>[ATOM_INDEX]:[GROUP_LABEL]</a>`: atom-indexed substituent or abbreviation.
   - `<a>[ATOM_INDEX]:<id>[NOTE]</a>`: atom-indexed special Markush label. The drawer shows `NOTE` as group text unless the indexed atom is `*` and `NOTE` is an approved endpoint-ball token.
   - `<r>[RING_INDEX]:[GROUP_LABEL]</r>`: ring-indexed substituent with unspecified attachment atom.
   - `<c>[ATOM_INDEX]:[RING_LABEL]</c>`: abstract-ring or superatom placeholder at a dummy atom.
   - `<d>[ATOM_INDEX]:<dum></d>`: explicit dummy attachment point. This is the E-SMILES 2.0 form; legacy `<a>[ATOM_INDEX]:<dum></a>` is still accepted.
   - `<v>[VIRTUALARC_INDEX]:[VIRTUALARC_NAME]:[FROM_ATOM:TO_ATOM]</v>`: pre-compatible virtualArc annotation for a special abstract ring; the name field may be empty.
   - `<r><v>[VIRTUALARC_INDEX]:[GROUP_LABEL]</r>`: pre-compatible substituent attached to a virtualArc.
   - `<x>[ATOM_1:ATOM_2]:[STATE]</x>`: biphenyl axial chirality, where the directly bonded aromatic axis atoms satisfy `ATOM_1 < ATOM_2` and `STATE` is `Ra` or `Sa`.
   - `<s>[SUBSTRUCTURE_ESMILES]</s>`: pre-compatible nested substructure record for ring-external repeat fragments.
   - `<g>[INNER_PORT:OUTER_PORT]:...:|Sg:n|</g>`: pre-compatible s-group repeat record.
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
- Axial chirality `<x>` is retained in `result["esmi"]`, with both axis atom
  indexes remapped after abbreviation substitution and canonicalization.
  `result["smi"]` ignores this annotation and remains ordinary SMILES.

## Markush Definition Substitution

Use `substitute_markush` when unresolved Markush labels or repeat counts should
be resolved. Definition keys can use `R1` or `R[1]`; values can be known
abbreviations (`Me`, `Cl`, `Ph`) or SMILES fragments. If a fragment contains
`*`, that atom is used as the attachment point. The return value is a single
RDKit-parseable SMILES string when the branch is fully determined, a list for
multiple enumerated structures, or residual E-SMILES when best-effort expansion
cannot infer a connection without guessing.

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
- Treat `<id>[NOTE]` as a preserved special label, not as a definition key.
  `substitute_markush` may resolve other labels in the same branch while keeping
  the `<id>` record in residual E-SMILES.
- Labels inside a single-level pre-compatible `<s>` record are substituted by `substitute_markush`, and the record is returned as a preserved E-SMILES annotation, e.g. `<s>**<sep><a>0:L[3]</a><a>1:R[3]</a>|Sg:n|</s>` can become `<s>ClBr<sep>|Sg:n|</s>`. An `<s>` nested inside another `<s>` is rejected. This does not change the existing `?n` copy behavior.
- With the default `repeat_policy="best_effort"`, physically expand Whole-SRU and two-port `<g>` repeats only when their boundary bonds and scope are unambiguous. Preserve under-specified `<s>`, `<g>`, `<v>`, and `<r><v>...` portions as E-SMILES, while substituting independent resolvable labels and remapping retained outer-graph indices. Use `strict` to reject residual repeats and `preserve` for annotation-only repeat handling.
- For the PEG-like example `*OCCO*<sep><d>0:<dum></d><d>5:<dum></d><g>[1:0]:[3:4]:|Sg:12|</g>`, pass `terminal_policy="hydrogen"`; the result is `"O" + "CCO" * 12` (C24H50O13), a plain RDKit-parseable SMILES.

## Rendering

Use rendering for visual QA, not as chemical validation:

```python
from molparser import utils as mutils

svg_text = mutils.draw(result["esmi"], output_format="svg")
png_bytes = mutils.draw(result["esmi"], output_format="png")
```

The drawer displays atom substituents, dummy attachment points, abstract rings,
and ring-level annotations from the E-SMILES extension. For atom-indexed
`<id>[NOTE]`, it draws `NOTE` as group text by default. It draws an endpoint
ball only when the indexed base atom is `*` and `NOTE` is one of `ball`,
`grey`, `black`, `green`, `blue`, `yellow`, `purple`, `orange`, `pink`, or
`brown`. Matching is exact and case-sensitive. Setting
`features.endpoint_balls=False` disables the ball treatment and keeps the text
label.

For representative fixed-color endpoint balls, Whole SRU, local `<g>` repeats,
and a large stereochemical atom-Markush structure, read
[`figure-index.md`](figure-index.md). Draw each raw E-SMILES fixture. Do not call
`substitute_markush` first unless the user has explicitly requested concrete
Markush expansion.

## Validation Priorities

- Exactly one top-level `<sep>`.
- Nested `<s>` records may contain an additional `<sep>` for the substructure E-SMILES.
- `<a>` indexes atoms; `<d>` indexes explicit dummy attachment points; `<r>` indexes rings; `<c>` indexes the dummy atom carrying the abstract-ring label.
- `<v>` indexes virtualArc annotations in a separate namespace; `<r><v>0:R[3]</r>` attaches an unresolved group to virtualArc `0`.
- `<x>[ATOM_1:ATOM_2]:Ra</x>` and `<x>[ATOM_1:ATOM_2]:Sa</x>` identify the directly bonded aromatic atoms defining a biphenyl stereogenic axis. Both indexes use the base-SMILES atom namespace and must be written in ascending order.
- A virtualArc name may be empty. For new datasets, emit canonical endpoint pairs with `start < end`, sort multiple pairs lexicographically, assign consecutive ids, and rebind `<r><v>` references; continue accepting legacy reversed endpoints when reading.
- `<g>` port pairs use `[INNER_PORT:OUTER_PORT]`; repeat count defaults to `n` and may be explicit, e.g. `|Sg:20|`.
- `GROUP_LABEL` may be a common abbreviation (`Me`, `OMe`, `CF3`) or a Markush label (`R[1]`). Use `<id>[NOTE]` only as an atom-indexed `<a>` payload, with a non-empty, whitespace-free `NOTE` that contains no `]` and has no multiplicity suffix. Endpoint-ball rendering on `*` is only a visual treatment of approved `<id>` values; it is not a separate token or chemical identity. For dummy attachment points, prefer `<d>[ATOM_INDEX]:<dum></d>` and accept legacy `<a>[ATOM_INDEX]:<dum></a>`.
- Use local substructure multiplicity suffixes (`?n`, `?1-3`, `?3`) separately from SRU-level `|Sg:n|`.
- After canonicalization or abbreviation substitution, remap supported retained `<g>`, `<v>`, `<x>`, and atom/ring annotations. If a resolved substitution removes an endpoint needed by a retained annotation, best-effort mode rolls that branch back instead of emitting a silently corrupted record.

## Boundary Policy

- Do not invent tokens for unsupported coordination bonds, electron-transfer arrows, uncertain bond styles, or uncertain chirality.
- Preserve the encodable molecular backbone and report unencoded chemistry explicitly.

## Reference Order

1. `extended-smiles-spec.md`
2. `figure-index.md`
3. `validate_esmiles.py`
4. `source-provenance.md`
