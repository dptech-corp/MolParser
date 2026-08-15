# Audited Synthetic E-SMILES Examples

Use this reference when generating, explaining, or regression-testing Whole SRU,
local repetition, virtual arcs, or colored endpoint balls. The compact fixture
manifest is [`synthetic-examples.json`](synthetic-examples.json).

## Contents

- [Status and provenance](#status-and-provenance)
- [E-SMILES 1.0 and 2.0 compatibility](#e-smiles-10-and-20-compatibility)
- [Raw drawing versus substitution](#raw-drawing-versus-substitution)
- [Ten fixture cases](#ten-fixture-cases)
- [Integration gates](#integration-gates)

## Status and provenance

The ten records are a deliberately small documentation set, not a copy of the
150,000-record training dataset. Eight records were selected from the audited
`arc_sru_2_150k_v19` source and two from the fixed-color endpoint-ball examples.
The manifest records the source record, candidate id where applicable, and
source SVG SHA-256 for every fixture.

Generate the SVG documentation assets only through the pinned workflow. Run
`scripts/generate_synthetic_examples.py --check-only` before rendering, then
supply the exact drawer SHA and write a machine report:

```text
python -B scripts/generate_synthetic_examples.py --render --renderer-sha256 DRAWER_SHA256 --report QA_OUTPUT/render-report.json
```

The report binds the manifest and drawer hashes to each asset SHA, byte size,
and viewBox. Keep it as a CI or review artifact rather than an additional skill
reference. After visual review, copy the final asset hashes into the manifest,
change its status to `frozen`, and rerun the document tests.

## E-SMILES 1.0 and 2.0 compatibility

| Feature | E-SMILES 1.0-compatible handling | E-SMILES 2.0 form |
| --- | --- | --- |
| Whole SRU endpoints | Keep the terminal `*` atoms and `|Sg:n|`; endpoint roles are implicit. Legacy `<a>INDEX:<dum></a>` may also be encountered. | Prefer explicit `<d>INDEX:<dum></d>` records plus `|Sg:COUNT|`. |
| One-atom local repeat | No lossless token; preserve the base graph and report the repeat separately. | `<a>INDEX:CH2?COUNT</a>`; `COUNT` may be a number or symbol such as `n`, `m`, `x`, or `y`. |
| Multi-atom local repeat | No lossless token. | `<g>[INNER:OUTER]:[INNER:OUTER]:|Sg:COUNT|</g>`. Bracket shape is a drawing option, not a different E-SMILES token. |
| Virtual arc | No lossless token. | `<v>ID:NAME:[START:END]</v>` plus optional `<r><v>ID:GROUP</r>`. `NAME` may be empty. |
| Colored endpoint ball | A legacy ordinal `<id>[N]` does not preserve the color semantics without an external note. | MolParser project extension used in 2.0 workflows: an atom-indexed semantic label such as `<a>3:<id>[blue]</a>`. |

Do not silently downgrade a 2.0-only feature. If 1.0 output is required, return
the base graph and an explicit note describing the omitted repeat, arc, or color.

For newly generated virtual-arc data, normalize every endpoint pair to
`START < END`, sort multiple pairs lexicographically, and assign consecutive
arc ids `0..n-1` after sorting. Keep `<r><v>` references bound to the reassigned
ids. General readers may continue to accept legacy reversed endpoint order.

## Raw drawing versus substitution

Draw these fixtures from the raw E-SMILES:

```python
from molparser import utils as mutils

svg = mutils.draw(raw_esmiles, config=render_config, output_format="svg")
```

This preserves the intended annotation identity and visual contract. In
particular, it preserves dummy endpoints, repeat brackets, virtual-arc names and
references, and semantic ball labels.

Use `postprocess_caption` when normalized E-SMILES, SMILES, or CXSMILES is the
requested product. Use `substitute_markush` only when the user supplies concrete
definitions for unresolved atom- or ring-level Markush groups. Substitution is
not a preprocessing step for these drawing fixtures:

- `<g>`, `<v>`, and `<r><v>` are annotations, not instructions to build a new
  concrete molecular graph; `substitute_markush` may resolve a symbolic `<g>`
  count to one positive integer, but it does not physically repeat the graph;
- `<id>[blue]` and `<id>[yellow]` are rendering semantics, not Markush
  definition keys;
- local multiplicity expansion can change the structure and must be requested
  explicitly rather than inferred during drawing.

## Ten fixture cases

### 1. Whole SRU

```text
*CSCC(*)CO<sep><d>0:<dum></d><d>5:<dum></d>|Sg:n|
```

The two `<d>` records identify the dummy endpoints. The top-level `|Sg:n|`
marks the whole backbone between them as the structural repeating unit.

### 2. One unlabeled CH2 repeat

```text
CCCO*CC(CCl)OCCO<sep><a>4:CH2?16</a>
```

Atom `4` is the repeated CH2 vertex and `16` is the copy count. The publication
drawing should show the carbon as an unlabeled skeletal vertex inside round
parentheses; it should not print the literal text `CH2`.

### 3. Short local repeat with round brackets

```text
CCC(CO)CCCCCC(CC)C(C)O<sep><g>[6:5]:[9:10]:|Sg:3|</g>
```

The first number in each port pair is inside the repeated substructure and the
second is outside. Both outside boundary atoms are unlabeled skeletal carbons,
so the lower-right repeat count `3` remains visually distinct from atom labels.
Render this fixture with `styling.sgroup_bracket_style="round"`.

### 4. Complex local repeat with square brackets

```text
O=C(CCC1CCC1)OCCC1COCCN1c1c2c(nc3ccc(Cl)cc13)CCCC2<sep><g>[2:3]:[9:10]:|Sg:10|</g>
```

This fixture covers a longer heteroatom-containing local repeat. Render it with
`styling.sgroup_bracket_style="square"`.

### 5. Named arc without arc substituents

```text
C#CCC(C)OC(=O)C(N)CC#C<sep><v>0:Y:[0:12]</v>
```

Arc `0` is named `Y`; no `<r><v>` record is present.

### 6. Named arcs with arc substituents

```text
C#CCCc1cc(OCCCC=C)cc(C)c1C(=O)OCCN1CCC(CC=C)C1<sep><v>0:D:[0:12]</v><v>1:M:[12:28]</v><r><v>0:R[2]</r><r><v>1:R[2]</r>
```

This case exercises two sorted arcs with consecutive ids. Each `R[2]` record is
bound to its own arc id.

### 7. Unnamed arc with an arc substituent

```text
C=CCOCC(O)C(=O)OCC=C<sep><v>0::[0:12]</v><r><v>0:R[1]</r>
```

The empty field between the two colons means that the arc has no visible name.
The `R[1]` label remains attached to arc `0`.

### 8. Unnamed arc without arc substituents

```text
C=CCc1n[nH]c2c(O)nc(C=C)nc12<sep><v>0::[0:12]</v>
```

This is the minimal unnamed/no-reference appearance. Do not invent an arc name
or add a placeholder label while rendering it.

### 9. Blue endpoint ball

```text
CC(=O)*<sep><a>3:<id>[blue]</a>
```

The semantic note replaces dummy atom `3` with the approved fixed blue ball.
The source palette is fixed; documentation generation must disable color
variation.

### 10. Yellow endpoint ball

```text
c1ccccc1*<sep><a>6:<id>[yellow]</a>
```

This case verifies the same endpoint-ball rule on an aromatic attachment.

## Integration gates

Before committing the ten SVG assets, require all of the following:

1. The manifest passes `scripts/generate_synthetic_examples.py --check-only`.
2. The validator accepts named and unnamed virtual arcs, including either
   presence or absence of `<r><v>` records.
3. The drawer renders `<g>` as molecular brackets rather than a text note and
   honors both `round` and `square` bracket styles.
4. The single-CH2 fixture has no visible `CH2` glyph and places the repeat count
   at the lower right of the closing parenthesis.
5. Arc endpoints attach only to the intended unsaturated atoms; no arc, label,
   substituent connector, or molecular structure intersects another element.
6. The drawer renders approved endpoint-ball tokens in their fixed colors.
7. Every SVG parses as XML, contains no scripts or external references, and has
   a non-empty viewBox. Record each final asset SHA-256 in the manifest.
