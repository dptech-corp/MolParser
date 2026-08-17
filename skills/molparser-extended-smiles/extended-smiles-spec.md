# MolParser E-SMILES Spec

This spec follows the current `utils` implementation and the documented MolParser figure set.

## 1. Top-Level Format

```text
SMILES<sep>EXTENSION
```

- `SMILES`: RDKit-compatible molecular backbone.
- `<sep>`: required separator.
- `EXTENSION`: optional annotation string; keep `<sep>` even when empty.

For ordinary molecules, use `SMILES<sep>`.

## 2. Core Tokens

### Atom-Indexed Substituent

```text
<a>[ATOM_INDEX]:[GROUP_LABEL]</a>
```

- `ATOM_INDEX`: zero-based atom index in the base SMILES.
- `GROUP_LABEL`: substituent, abbreviation, or Markush placeholder. Legacy inputs may also use `<dum>` here for dummy attachment points.
- Atom-indexed special Markush labels use `<id>[NOTE]`, where `NOTE` is a
  non-empty, whitespace-free custom remark. Use this payload only inside
  `<a>...</a>`;

Example:

```text
*c1ccccc1<sep><a>0:R[1]</a>
```

```text
N(*)(*)C1=NC(NC(*)C(=O)N*)=NC(N(*)***C(=O)N(***)*)=N1<sep><a>12:<id>[DNA]</a>
```

### Ring-Indexed Substituent

```text
<r>[RING_INDEX]:[GROUP_LABEL]</r>
```

- `RING_INDEX`: zero-based ring index, independent from atom indexes.
- Use this when the substituent is attached to a ring but the exact attachment atom is regio-uncertain.

Example:

```text
c1ccccc1<sep><r>0:R[1]</r><r>0:R[2]</r>
```

### Abstract-Ring / Superatom Placeholder

```text
<c>[ATOM_INDEX]:[RING_LABEL]</c>
```

- `ATOM_INDEX`: zero-based dummy atom index carrying the abstract ring or superatom.
- `RING_LABEL`: abstract-ring label such as `B` or `Ar`.

Example:

```text
*C(NC(*)(*)C(*)(*)*)C(=O)N(*)*<sep><c>9:B</c>
```

### Dummy Attachment Point

`<d>[ATOM_INDEX]:<dum></d>` marks an explicit dummy atom attachment point. This is the E-SMILES 2.0 representation. The legacy `<a>[ATOM_INDEX]:<dum></a>` form remains valid for backward compatibility.

```text
*C(O)=O<sep><d>0:<dum></d>
```

### VirtualArc

`<v>[VIRTUALARC_INDEX]:[VIRTUALARC_NAME]:[FROM_ATOM:TO_ATOM]</v>` records a special abstract ring connection.

- `VIRTUALARC_INDEX`: zero-based index in a namespace separate from atoms and rings.
- `VIRTUALARC_NAME`: source label such as `A`, `Ar`, or `M`; it may be empty when the source arc is unnamed, for example `<v>0::[0:2]</v>`.
- `FROM_ATOM` and `TO_ATOM`: zero-based base-SMILES atom indexes for the endpoints. General readers may accept `[0:2]` and legacy `[2:0]` as equivalent. New datasets should emit the canonical order described below.
- A substituent attached to a virtualArc uses `<r><v>[VIRTUALARC_INDEX]:[GROUP_LABEL]</r>`.
- Current utilities preserve these records as pre-compatible annotations. The drawer renders their virtual-arc geometry, while normalization and Markush substitution do not turn the annotation into a chemical bond.

```text
C=CCC(C(C)*)*<sep><a>6:R[2]</a><a>7:R[1]</a><v>0:A:[0:2]</v><r><v>0:R[3]</r>
```

## 3. Multiplicity And Structural Repetition

Use a suffix on a group label for local substructure multiplicity:

```text
<r>1:R[1]?1-3</r>
<r>1:R[5]?n</r>
<a>0:CH2?3</a>
```

Use `|Sg:n|` for structural repeating unit (SRU) repetition:

```text
*CC*<sep><d>0:<dum></d><d>2:<dum></d>|Sg:n|
```

Current `molparser.utils.postprocess_caption` recognizes syntactically valid
`|Sg:COUNT|` forms, including symbolic, numeric, and simple range counts, as SRU
markers.

Use `<s>...</s>` for a nested substructure record, typically a ring-external repeat fragment. The body is another `SMILES<sep>EXTENSION` fragment.

```text
N1=C(N)C2=C(C(=C(N2)*)***)N=C1**<sep><a>8:R[2]</a><a>9:L[2]</a><a>10:B</a><s>**<sep><a>0:L[3]</a><a>1:R[3]</a>|Sg:n|</s><a>14:L[1]</a><a>15:R[1]</a>
```

Use `<g>[INNER_PORT:OUTER_PORT]:...:|Sg:n|</g>` for an s-group repeat record.

- `INNER_PORT`: atom index inside the repeat structure.
- `OUTER_PORT`: atom index outside the repeat structure.
- Multiple port pairs are written by repeating `[INNER_PORT:OUTER_PORT]` and separating fields with `:`.
- The repeat count defaults to `n` when abstract; explicit numbers, letters, and ranges such as `20`, `m`, or `m-n` are syntactically accepted.
- `substitute_markush` defaults to best-effort repeat handling. It physically expands a concrete Whole-SRU or two-port `<g>` repeat only when the boundary bonds, cut set, and scope are unambiguous; a fully resolved result is RDKit-parseable SMILES. Otherwise it preserves the unresolved record as E-SMILES and still substitutes independent resolvable labels. Use `repeat_policy="strict"` to reject such residuals or `"preserve"` for annotation-only repeat handling.
- Round parentheses and square brackets are equivalent drawing styles for the same `<g>` record; bracket shape is not encoded in E-SMILES.

```text
C=CCC(C(C)*)*<sep><a>6:R[2]</a><a>7:R[1]</a><g>[3:2]:[4:5]:|Sg:20|</g>
```

## 4. Labels

- Common abbreviations: `Me`, `OMe`, `Ph`, `CF3`.
- Markush placeholders: `R[1]`, `R[2]`, `X[1]`.
- Atom-indexed special Markush labels use one grammar:
  `<a>[ATOM_INDEX]:<id>[NOTE]</a>`. Do not use `<id>` as the payload of `<d>`,
  `<r>`, `<c>`, or `<r><v>`. The note is opaque payload text rather than a
  Markush definition key, so `substitute_markush` preserves it while resolving
  independent labels where possible.
- The drawer renders `NOTE` as ordinary group text by default. The approved
  values `ball`, `grey`, `black`, `green`, `blue`, `yellow`, `purple`,
  `orange`, `pink`, and `brown` select endpoint-ball rendering only when
  `ATOM_INDEX` points to `*`. Matching is exact and case-sensitive. Otherwise
  the value remains group text. Setting `features.endpoint_balls` to `false`
  also keeps the text form. This is a MolParser visual extension of the same
  `<id>[NOTE]` token, not a separate E-SMILES grammar or chemical identity.
- Dataset-specific labels may appear as payload text when they cannot be reduced to a standard abbreviation.

## 5. Utility Behavior

- `molparser.utils.postprocess_caption` / `Translator.refactor` canonicalize SMILES and substitute known atom-indexed abbreviations from `molparser/utils/abbrevs_example.csv` when the attachment is chemically valid.
- Resolved substituents are folded into the base `smi`; unresolved Markush or ring-level annotations remain in `groups`.
- Atom-indexed `<id>[NOTE]` records remain in `groups`; they are not looked up
  in the Markush definition dictionary. Non-ball notes render as text, while
  only approved values on `*` in Section 4 render as endpoint balls.
- Pre-compatible `<s>` records are preserved in `groups`; `substitute_markush` can substitute labels inside one single-level substructure record and return the result as an E-SMILES annotation, without attaching or expanding it into the main molecule. An `<s>` nested inside another `<s>` is outside the supported grammar and is rejected.
- Pre-compatible `<s>` and under-specified `<g>`, `<v>`, or `<r><v>...` portions remain in `groups`. When other substitutions reorder the outer graph, retained `<g>` ports, `<v>` endpoints, and `<c>` atom indices are remapped; if an endpoint was removed, best-effort mode rolls back that branch instead of emitting a stale index. A symbolic `<g>` count may resolve to one positive integer even when its graph remains residual.
- `draw` renders SMILES or E-SMILES to SVG/PNG for visual QA, including SRU and local-repeat brackets, virtual arcs, and approved MolParser endpoint-ball labels.

## 6. Unsupported Chemistry

The current token set does not encode:

- coordination or dative-bond semantics beyond the representational scope of standard SMILES (e.g., metal complexes);
- electron-transfer arrows;
- uncertain bond styles;
- uncertain chirality.

Preserve the encodable backbone, do not invent tokens, and report unencoded chemistry explicitly.

## 7. Validation Checklist

- at least one top-level `<sep>`; nested `<s>` records may contain their own `<sep>`;
- balanced `<a>`, `<d>`, `<r>`, `<c>`, `<s>`, `<g>`, and `<v>` tags, allowing inline `<r><v>...</r>`;
- non-negative indexes in the correct namespace;
- no whitespace inside group labels;
- special labels use exactly `<a>[ATOM_INDEX]:<id>[NOTE]</a>` with a non-empty,
  whitespace-free `NOTE` that contains no `]` and has no `?` multiplicity
  suffix; endpoint-ball values use the same syntax as text labels and render as
  balls only when the indexed base atom is `*` and ball rendering is enabled;
- local substructure multiplicity uses `?n`, `?1-3`, or `?3`;
- SRU repetition uses `|Sg:n|` or an explicit count such as `|Sg:20|`;
- supported retained `<g>`, `<v>`, and atom/ring indexes are remapped after canonicalization or substitution; if an endpoint is removed and cannot be remapped, best-effort mode preserves the original branch instead of silently changing the annotated scope.
- an unnamed virtualArc uses an empty name field, `<v>0::[START:END]</v>`;
