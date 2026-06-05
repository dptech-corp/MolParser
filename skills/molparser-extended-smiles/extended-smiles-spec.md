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
- `GROUP_LABEL`: substituent, abbreviation, Markush placeholder, or `<dum>`.

Example:

```text
*c1ccccc1<sep><a>0:R[1]</a>
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

`<dum>` marks an explicit dummy atom attachment point.

```text
*C(O)=O<sep><a>0:<dum></a>
```

## 3. Repetition

Use a suffix on a group label for local multiplicity:

```text
<r>1:R[1]?1-3</r>
<r>1:R[5]?n</r>
<a>0:CH2?3</a>
```

Use `|Sg:n|` for structural repeating unit (SRU) repetition:

```text
*CC*<sep><a>0:<dum></a><a>2:<dum></a>|Sg:n|
```

Current `utils.postprocess_caption` recognizes `|Sg:n|` as the SRU marker.

## 4. Labels

- Common abbreviations: `Me`, `OMe`, `Ph`, `CF3`.
- Markush placeholders: `R[1]`, `R[2]`, `X[1]`.
- Dataset-specific labels may appear as payload text when they cannot be reduced to a standard abbreviation.

## 5. Utility Behavior

- `postprocess_caption` / `Translator.refactor` canonicalize SMILES and substitute known atom-indexed abbreviations from `utils/abbrevs_example.csv` when the attachment is chemically valid.
- Resolved substituents are folded into the base `smi`; unresolved Markush or ring-level annotations remain in `groups`.
- `draw` renders SMILES or E-SMILES to SVG/PNG for visual QA.

## 6. Unsupported Chemistry

The current token set does not encode:

- coordination or dative bond semantics beyond ordinary SMILES support;
- electron-transfer arrows;
- uncertain bond styles;
- uncertain chirality.

Preserve the encodable backbone, do not invent tokens, and report unencoded chemistry explicitly.

## 7. Validation Checklist

- exactly one top-level `<sep>`;
- balanced `<a>`, `<r>`, and `<c>` tags;
- non-negative indexes in the correct namespace;
- no whitespace inside group labels;
- local multiplicity uses `?n`, `?1-3`, or `?3`;
- SRU repetition uses `|Sg:n|`;
- extension indexes are regenerated after canonicalization or substitution.
