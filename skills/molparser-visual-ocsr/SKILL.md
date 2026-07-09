---
name: molparser-visual-ocsr
description: "Use for MolParser visual OCSR workflows that parse molecule images or PDFs into E-SMILES/SMILES with MolDet detection and MolParser recognition models."
---

# MolParser Visual OCSR Skill

Use this skill when extracting molecular structures from images, PDFs, URLs, or mixed input lists with `molparser.MolParser`.

## Workflow

1. Use `molparser/models/config.yaml` defaults or pass a custom YAML path to `MolParser`.
2. Prefer `device="auto"` unless CPU/GPU must be forced.
3. For PIL images, image paths, or image URLs:
   - Use `rec_only=True` when the image is already a single molecule crop.
   - Use `rec_only=False` when detection is needed; this uses `moldet_v2_yolo11n_640_general.pt`.
4. For PDF paths or PDF URLs:
   - `rec_only` is ignored; render pages with `pypdfium2` at `pdf_dpi=300`.
   - Detect molecules with `moldet_v2_yolo11n_960_doc.pt`.
5. Recognize crops with the MolParser recognition model, then normalize with `postprocess_caption`.
6. Use batch defaults unless memory requires tuning: MolDet batch size is 8, MolParser OCSR batch size is 16.
7. Multi-image inputs are batched together; PDF page detections are cropped first, then molecule crops are batched together for recognition.
8. Use `moldet.expand_px` for detector box expansion and `molparser.padding_px` for recognition input padding; defaults are `expand_px=2` and `padding_px=0`, and both can be overridden in `parse(...)`.

## Minimal Example

```python
from molparser import MolParser

parser = MolParser(device="auto", pdf_dpi=300)

results = parser.parse(["mol.png", "paper.pdf"], rec_only=False)
for item in results:
    print(item.input_index, item.page_index, item.bbox, item.esmi)
```

For high-throughput molecule crops:

```python
parser = MolParser(device="cuda", molparser_batch_size=1024)
results = parser.parse(["mol_001.png", "mol_002.png"], rec_only=True)
```

## Model Download Policy

- Local paths win when configured.
- Hugging Face is tried first when the HF model ID is set.
- If Hugging Face fails and the ModelScope ID is set, automatically fall back to ModelScope.
- Default MolDet files are `moldet_v2_yolo11n_640_general.pt` for image detection and `moldet_v2_yolo11n_960_doc.pt` for PDF page detection.
- Use environment variables for proxy and private tokens:

```bash
export HTTP_PROXY=100.68.165.249:3128
export HTTPS_PROXY=100.68.165.249:3128
export HF_TOKEN=hf_xxx
```

## Output Fields

Each result includes `source`, `input_index`, `page_index`, `bbox`, `confidence`, `raw_caption`, `caption`, `smi`, `esmi`, `cxsmiles`, `markush`, `sru`, and `groups`.

- All parse methods return `list[MolParserResult]`.
- Single image with `rec_only=True`: one whole-image result; `bbox`, `confidence`, and `page_index` are `None`.
- Image with `rec_only=False`: one result per detected molecule; no detections returns an empty list.
- PDF input: one result per detected molecule on rendered PDF pages; `page_index`, `bbox`, and `confidence` are populated.
- List input: flat result list; use `input_index` to map each result back to the original list item.

## License Notes

MolDetv2 weights are non-commercial according to their model card. The PyTorch detection path uses Ultralytics YOLO; closed-source, private, SaaS, internal business, or commercial usage may require AGPL-3.0 compliance or an Ultralytics Enterprise license.
