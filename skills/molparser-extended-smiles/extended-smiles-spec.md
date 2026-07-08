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
- Special Markush labels use `<id>[NOTE]`, where `NOTE` is a custom remark such as `DNA`, `RNA`, `protein`, `red`, or `ball`.

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

`<d>[ATOM_INDEX]:<dum></d>` marks an explicit dummy atom attachment point. This is the new SMILES 2.0 representation. The legacy `<a>[ATOM_INDEX]:<dum></a>` form remains valid for backward compatibility.

```text
*C(O)=O<sep><d>0:<dum></d>
```

### VirtualArc

`<v>[VIRTUALARC_INDEX]:[VIRTUALARC_NAME]:[FROM_ATOM:TO_ATOM]</v>` records a special abstract ring connection.

- `VIRTUALARC_INDEX`: zero-based index in a namespace separate from atoms and rings.
- `VIRTUALARC_NAME`: source label such as `A`, `Ar`, or `M`.
- `FROM_ATOM` and `TO_ATOM`: zero-based base-SMILES atom indexes for the unordered endpoints; `[0:2]` and `[2:0]` are equivalent.
- A substituent attached to a virtualArc uses `<r><v>[VIRTUALARC_INDEX]:[GROUP_LABEL]</r>`.
- Current utilities preserve these records as pre-compatible annotations; they do not expand or render the virtualArc geometry.

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

Current `molparser.utils.postprocess_caption` recognizes `|Sg:n|` as the SRU marker.

Use `<s>...</s>` for a nested substructure record, typically a ring-external repeat fragment. The body is another `SMILES<sep>EXTENSION` fragment.

```text
N1=C(N)C2=C(C(=C(N2)*)***)N=C1**<sep><a>8:R[2]</a><a>9:L[2]</a><a>10:B</a><s>**<sep><a>0:L[3]</a><a>1:R[3]</a>|Sg:n|</s><a>14:L[1]</a><a>15:R[1]</a>
```

Use `<g>[INNER_PORT:OUTER_PORT]:...:|Sg:n|</g>` for an s-group repeat record.

- `INNER_PORT`: atom index inside the repeat structure.
- `OUTER_PORT`: atom index outside the repeat structure.
- Multiple port pairs are written by repeating `[INNER_PORT:OUTER_PORT]` and separating fields with `:`.
- The repeat count defaults to `n` when abstract; explicit numbers, letters, and ranges such as `20`, `m`, or `m-n` are syntactically accepted.
- Current utilities preserve `<s>` and `<g>` records as pre-compatible annotations; they do not expand repeat structures.

```text
C=CCC(C(C)*)*<sep><a>6:R[2]</a><a>7:R[1]</a><g>[3:2]:[4:5]:|Sg:20|</g>
```

## 4. Labels

- Common abbreviations: `Me`, `OMe`, `Ph`, `CF3`.
- Markush placeholders: `R[1]`, `R[2]`, `X[1]`.
- Special Markush labels: `<id>[DNA]`, `<id>[RNA]`, `<id>[protein]`, or other task-specific notes.
- Dataset-specific labels may appear as payload text when they cannot be reduced to a standard abbreviation.

## 5. Utility Behavior

- `molparser.utils.postprocess_caption` / `Translator.refactor` canonicalize SMILES and substitute known atom-indexed abbreviations from `molparser/utils/abbrevs_example.csv` when the attachment is chemically valid.
- Resolved substituents are folded into the base `smi`; unresolved Markush or ring-level annotations remain in `groups`.
- Pre-compatible `<s>` records are preserved in `groups`; `substitute_markush` can recursively substitute labels inside the nested substructure and return the result as an E-SMILES annotation, without attaching or expanding it into the main molecule.
- Pre-compatible `<g>`, `<v>`, and `<r><v>...` records are preserved in `groups` but are not chemically expanded, remapped, or specially rendered yet.
- `draw` renders SMILES or E-SMILES to SVG/PNG for visual QA.

## 6. Unsupported Chemistry

The current token set does not encode:

- coordination or dative bond semantics beyond ordinary SMILES support;
- electron-transfer arrows;
- uncertain bond styles;
- uncertain chirality.

Preserve the encodable backbone, do not invent tokens, and report unencoded chemistry explicitly.

## 7. Validation Checklist

- at least one top-level `<sep>`; nested `<s>` records may contain their own `<sep>`;
- balanced `<a>`, `<d>`, `<r>`, `<c>`, `<s>`, `<g>`, and `<v>` tags, allowing inline `<r><v>...</r>`;
- non-negative indexes in the correct namespace;
- no whitespace inside group labels;
- local substructure multiplicity uses `?n`, `?1-3`, or `?3`;
- SRU repetition uses `|Sg:n|` or an explicit count such as `|Sg:20|`;
- supported extension indexes are regenerated after canonicalization or substitution; pre-compatible `<s>`, `<g>`, and `<v>` annotations are preserved without semantic remapping.
