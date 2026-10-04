# MolParser E-SMILES 2.0 Spec

This document describes E-SMILES 2.0 notation and the current `utils` behavior.
Repository compatibility extensions and normalization conventions are identified
separately from the core notation.

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
  `<a>...</a>`. The literal form is `<a>0:<id>[DNA]</a>`; `<id>` is an
  inline marker and has no closing `</id>` token.

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
<r>[RING_INDEX_1],[RING_INDEX_2],...:[GROUP_LABEL]</r>
```

- `RING_INDEX`: zero-based ring index, independent from atom indexes.
- Use this when the substituent is attached to a ring but the exact attachment atom is regio-uncertain.
- Comma-separated ring indexes name one substituent that may attach at any substitutable atom in any of the listed rings. The site set is the union of those rings. Duplicate indexes name the same ring once. This form is only for plain `<r>` records, not `<r><v>`.

Example:

```text
c1ccccc1<sep><r>0:R[1]</r><r>0:R[2]</r>
```

```text
c1ccc2ccccc2c1<sep><r>0,1:R[1]</r>
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

The index must point to `*` in the base SMILES. `<dum>` is an inline marker;
do not add `</dum>`.

```text
*C(O)=O<sep><d>0:<dum></d>
```

### VirtualArc

`<v>[VIRTUALARC_INDEX]:[VIRTUALARC_NAME]:[FROM_ATOM:TO_ATOM]</v>` records a special abstract ring connection.

- `VIRTUALARC_INDEX`: zero-based index in a namespace separate from atoms and rings.
- `VIRTUALARC_NAME`: source label such as `A`, `Ar`, or `M`; it may be empty when the source arc is unnamed, for example `<v>0::[0:2]</v>`.
- `FROM_ATOM` and `TO_ATOM`: zero-based base-SMILES atom indexes for the endpoints. The repository accepts `[0:2]` and `[2:0]` as equivalent. Its normalization convention is to emit ascending endpoint pairs, sort the pairs lexicographically, assign consecutive arc ids, and update `<r><v>` references. This is a repository convention, distinct from the required ascending axis atom indexes for `<x>`.
- A substituent attached to a virtualArc uses `<r><v>[VIRTUALARC_INDEX]:[GROUP_LABEL]</r>`.
- Virtual arcs are E-SMILES 2.0 annotations for abstract connections, not chemical bonds. The drawer renders their geometry, while normalization and Markush substitution preserve their annotation semantics without adding a chemical bond.

```text
C=CCC(C(C)*)*<sep><a>6:R[2]</a><a>7:R[1]</a><v>0:A:[0:2]</v><r><v>0:R[3]</r>
```

### Axial Chirality (Stereogenic Axis)

A stereogenic axis uses a single `<x></x>` record. The former biaryl-bond
endpoints are generalized to the terminal reference points of the chiral-axis
vector:

```text
<x>[ATOM_1:ATOM_2]:[STATE]</x>
```

- `ATOM_1` and `ATOM_2` are the zero-based key atom indexes that span the
  topology of the chiral axis. They must satisfy `ATOM_1 < ATOM_2`.
- `[ATOM_1:ATOM_2]` uniquely locks the skeleton that carries the axial
  chirality: the two atoms of a hindered single bond, or the key atoms at
  both ends of the principal axis in a multi-atom cumulene or orthogonal
  system.
- `STATE` is the absolute axial configuration:
  - `Ra`: R axial
  - `Sa`: S axial
- Atom order is not reversed to change the state.
- Normalized E-SMILES preserves this record and remaps both atom indexes after
  abbreviation substitution or canonicalization. Plain SMILES output ignores
  the annotation; current utilities do not assign RDKit axial stereochemistry.
- Definite axial configuration alone is not Markush. The following example has
  `markush=False` and `sru=False`; an independent unresolved label or virtual
  arc can still make a structure Markush.

```text
Nc1ccc2ccccc2c1-c1c(O)ccc2ccccc12<sep><x>[10:11]:Sa</x>
```

| Axial chirality type | Topology / atropisomerism | Axis-end atom rule | Example |
| --- | --- | --- | --- |
| Biaryl systems | Hindered single-bond rotation (C–C) | Atom indexes of the C–C single bond joining the two aryl rings | `<x>[5:12]:Ra</x>` |
| Aryl–alkenyl / diene systems | Hindered single-bond rotation (C–C) | Aryl carbon and alkenyl α-carbon indexes | `<x>[4:9]:Sa</x>` |
| Aryl amide / imide systems | Hindered single-bond rotation (C–N) | Aryl carbon and amide nitrogen indexes | `<x>[6:7]:Ra</x>` |
| Allene / cumulene systems | Orthogonal cumulative double-bond π-orbitals | Terminal sp² carbon indexes | `<x>[2:4]:Sa</x>` |
| Alkylidenecycloalkanes | Ring orthogonal to the double bond | Ring attachment atom and terminal alkene carbon indexes | `<x>[3:8]:Sa</x>` |
| Spiro compounds | Orthogonal geometry at the spiro atom | Distal atoms at both ends of the principal axis through the spiro atom (para vertices) | `<x>[1:10]:Ra</x>` |

## 3. Multiplicity And Structural Repetition

Use a suffix on a group label for local substructure multiplicity:

```text
<r>1:R[1]?1-3</r>
<r>1:R[5]?n</r>
<a>0:CH2?3</a>
```

Use top-level `|Sg:COUNT|` for whole-molecule structural repetition. When the
source identifies an unspecified whole repeat, write `|Sg:n|`:

```text
*CC*<sep><d>0:<dum></d><d>3:<dum></d>|Sg:n|
```

Symbolic, numeric, and simple range counts are valid repeat annotations. The
public `sru` classification is narrower: it is true only for a symbolic
whole-molecule count, not integers, numeric ranges, or local/nested repeats.

The repository additionally accepts `<s>...</s>` as a compatibility extension
for nested substructure records, typically ring-external repeat fragments. It
is not one of the core tokens described in the E-SMILES 2.0 announcement. The
body is another `SMILES<sep>EXTENSION` fragment with its own atom indexes.

```text
N1=C(N)C2=C(C(=C(N2)*)***)N=C1**<sep><a>8:R[2]</a><a>9:L[2]</a><a>10:B</a><s>**<sep><a>0:L[3]</a><a>1:R[3]</a>|Sg:n|</s><a>14:L[1]</a><a>15:R[1]</a>
```

Use `<g>[INNER_PORT:OUTER_PORT]:...:|Sg:n|</g>` for an E-SMILES 2.0 local
s-group repeat record. It does not itself set the whole-molecule `sru` flag.

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
- A comma-separated `<r>` record enumerates attachment atoms from the union of the listed rings. `<r>0,1:R[1]</r>` with `{"R1": "Me"}` returns the distinct canonical methyl products on either ring. Multiplicity such as `?1-2` counts copies on that same union. An unresolved list is preserved, for example `<r>0,1:R[1]</r>`, and each listed ring index is remapped after canonicalization. If any listed ring cannot be remapped, that record is dropped rather than narrowed to the rings that remain.
- `substitute_markush` checks complete user definitions first, then complete
  abbreviation matches, before attempting composite splitting. Supported
  linear prefixes include `SO2`, `CO2`, `CO`, `CH2`, `CF2`, `NH`, `N`, `O`,
  and `S`; multiple prefixes can precede a resolvable terminal group. In
  best-effort mode, partial atom-group expansion can retain labels such as
  `R`, `X`, `Y`, `Ar`, `Het`, or supported indexed forms on a new dummy atom,
  with annotation indexes remapped. Unknown forms that cannot be split remain
  whole labels. These are implementation heuristics, not new E-SMILES tokens.
- The automatic composite pass in `postprocess_caption` currently triggers on
  `SO2`, `CO2`, `CO`, `CH2`, `CF2`, or `NH`. For standalone composite forms
  such as `NR`, `OR`, and `SR`, call `substitute_markush` explicitly. For
  example, `*CCO<sep><a>0:CONHAr</a>` with no definitions becomes
  `*NC(=O)CCO<sep><a>0:Ar</a>` in best-effort mode. A fully resolved branch is
  returned as SMILES; unresolved branches may remain E-SMILES.
- Resolved substituents are folded into the base `smi`; unresolved Markush or ring-level annotations remain in `groups`.
- Atom-indexed `<id>[NOTE]` records remain in `groups`; they are not looked up
  in the Markush definition dictionary. Non-ball notes render as text, while
  only approved values on `*` in Section 4 render as endpoint balls.
- Compatibility `<s>` records are preserved in `groups`; `substitute_markush` can substitute labels inside one single-level substructure record and return the result as an E-SMILES annotation, without attaching or expanding it into the main molecule. An `<s>` nested inside another `<s>` is outside the supported grammar and is rejected.
- Retained `<s>`, unresolved `<g>`, virtual-arc records, and definite `<x>` annotations remain in `groups`. When other substitutions reorder the outer graph, retained `<g>` ports, `<v>` endpoints, `<x>` axis atoms, and `<c>` atom indices are remapped; if an endpoint was removed, best-effort mode rolls back that branch instead of emitting a stale index. A symbolic `<g>` count may resolve to one positive integer even when its graph remains residual.
- `draw` renders SMILES or E-SMILES to SVG/PNG for visual QA, including SRU and local-repeat brackets, virtual arcs, and approved MolParser endpoint-ball labels.

### Classification Flags

`markush` and `sru` describe different properties and may both be true. They do
not control annotation preservation or the scope of a rendered repeat.

| Situation, with no other annotations | `markush` | `sru` |
| --- | --- | --- |
| Fully specified molecule with only `<x>` | `False` | `False` |
| Unresolved atom/ring group or virtual arc | `True` | `False` |
| Whole repeat with top-level `|Sg:n|` or `|Sg:m|` | `False` | `True` |
| Whole repeat with top-level `|Sg:3|` or `|Sg:1-3|` | `False` | `False` |
| Symbolic local `<g>`, `?n`, or a repeat inside `<s>` | `True` | `False` |

The whole-repeat rows assume explicit terminal dummies and no unresolved
substituents. Numeric repeat annotations can remain valid even with
`sru=False`; post-processing expands a numeric repeat only when it can return
one deterministic structure. Use `esmi` to retain extensions such as `<x>`
that plain SMILES and current CXSMILES conversion do not encode.

## 6. Unsupported Chemistry

The current token set does not encode:

- coordination or dative-bond semantics beyond the representational scope of standard SMILES (e.g., metal complexes);
- electron-transfer arrows;
- uncertain bond styles;
- uncertain chirality, including axial chirality whose absolute configuration is not a definite `Ra` or `Sa`.

Preserve the encodable backbone, do not invent tokens, and report unencoded chemistry explicitly.

## 7. Validation Checklist

- exactly one top-level `<sep>`; nested `<s>` records may contain their own `<sep>`;
- balanced `<a>`, `<d>`, `<r>`, `<c>`, `<s>`, `<g>`, `<v>`, and `<x>` tags, allowing inline `<r><v>...</r>`;
- non-negative indexes in the correct namespace; a plain `<r>` record may list several ring indexes separated by commas;
- `<d>...:<dum></d>` and legacy `<a>...:<dum></a>` must index a dummy `*` atom;
- `<x>` uses two terminal reference atom indexes of a stereogenic axis in ascending order (`ATOM_1 < ATOM_2`) and a `Ra` or `Sa` state;
- no whitespace inside group labels;
- special labels use exactly `<a>[ATOM_INDEX]:<id>[NOTE]</a>` with a non-empty,
  whitespace-free `NOTE` that contains no `]` and has no `?` multiplicity
  suffix; endpoint-ball values use the same syntax as text labels and render as
  balls only when the indexed base atom is `*` and ball rendering is enabled;
- local substructure multiplicity uses `?n`, `?1-3`, or `?3`;
- SRU repetition uses `|Sg:n|` or an explicit count such as `|Sg:20|`;
- supported retained `<g>`, `<v>`, and atom/ring indexes are remapped after canonicalization or substitution; if an endpoint is removed and cannot be remapped, best-effort mode preserves the original branch instead of silently changing the annotated scope.
- an unnamed virtualArc uses an empty name field, `<v>0::[START:END]</v>`;
