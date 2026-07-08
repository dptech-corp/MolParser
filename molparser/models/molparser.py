from __future__ import annotations

import os
from contextlib import contextmanager
from dataclasses import asdict, dataclass, fields, replace
from importlib.resources import files
from pathlib import Path
from typing import Any, Optional

from .io import load_image, normalize_inputs, render_pdf
from .runtime import resolve_model
from ..utils import postprocess_caption


DEFAULT_CONFIG = "config.yaml"


@dataclass
class MolParserConfig:
    moldet_hf_repo: str = "UniParser/MolDetv2"
    moldet_modelscope_repo: str = "UniParser/MolDetv2"
    moldet_image_model_path: str = ""
    moldet_pdf_model_path: str = ""
    moldet_image_modelname: str = "moldet_v2_yolo11n_640_general.pt"
    moldet_pdf_modelname: str = "moldet_v2_yolo11n_960_doc.pt"
    molparser_hf_repo: str = "AI4Industry/MolParser-Mobile"
    molparser_modelscope_repo: str = ""
    molparser_model_path: str = ""
    cache_dir: str = ""
    device: str = "auto"
    detector_conf: float = 0.5
    image_imgsz: int = 640
    pdf_imgsz: int = 960
    moldet_batch_size: int = 8
    pdf_dpi: int = 300
    max_length: int = 256
    molparser_batch_size: int = 16
    expand_px: int = 2
    padding_px: int = 0
    hf_token: str | bool | None = None

    @classmethod
    def from_default(cls, **overrides) -> "MolParserConfig":
        return cls.from_yaml(files("molparser.models").joinpath(DEFAULT_CONFIG), **overrides)

    @classmethod
    def from_yaml(cls, path: str | Path, **overrides) -> "MolParserConfig":
        try:
            import yaml
        except ImportError as exc:
            raise ImportError("YAML config loading requires PyYAML.") from exc

        with Path(path).open("r", encoding="utf-8") as f:
            data = yaml.safe_load(f) or {}
        return cls.from_mapping(data, **overrides)

    @classmethod
    def from_mapping(cls, data: dict[str, Any], **overrides) -> "MolParserConfig":
        flat = {
            "moldet_hf_repo": data.get("moldet", {}).get("hf_repo", cls.moldet_hf_repo),
            "moldet_modelscope_repo": data.get("moldet", {}).get("modelscope_repo", cls.moldet_modelscope_repo),
            "moldet_image_model_path": data.get("moldet", {}).get("image_model_path", ""),
            "moldet_pdf_model_path": data.get("moldet", {}).get("pdf_model_path", ""),
            "moldet_image_modelname": data.get("moldet", {}).get("image_modelname", cls.moldet_image_modelname),
            "moldet_pdf_modelname": data.get("moldet", {}).get("pdf_modelname", cls.moldet_pdf_modelname),
            "detector_conf": data.get("moldet", {}).get("confidence", cls.detector_conf),
            "image_imgsz": data.get("moldet", {}).get("image_imgsz", cls.image_imgsz),
            "pdf_imgsz": data.get("moldet", {}).get("pdf_imgsz", cls.pdf_imgsz),
            "moldet_batch_size": data.get("moldet", {}).get("batch_size", cls.moldet_batch_size),
            "expand_px": data.get("moldet", {}).get("expand_px", cls.expand_px),
            "molparser_hf_repo": data.get("molparser", {}).get("hf_repo", cls.molparser_hf_repo),
            "molparser_modelscope_repo": data.get("molparser", {}).get("modelscope_repo", ""),
            "molparser_model_path": data.get("molparser", {}).get("model_path", ""),
            "max_length": data.get("molparser", {}).get("max_length", cls.max_length),
            "molparser_batch_size": data.get("molparser", {}).get("batch_size", cls.molparser_batch_size),
            "padding_px": data.get("molparser", {}).get("padding_px", cls.padding_px),
            "pdf_dpi": data.get("pdf", {}).get("dpi", cls.pdf_dpi),
            "device": data.get("runtime", {}).get("device", cls.device),
            "cache_dir": data.get("runtime", {}).get("cache_dir", ""),
            "hf_token": data.get("runtime", {}).get("hf_token"),
        }
        flat.update(_normalize_overrides(overrides))
        allowed = {field.name for field in fields(cls)}
        return cls(**{key: value for key, value in flat.items() if key in allowed})


@dataclass
class MolParserResult:
    source: str
    input_index: int
    page_index: Optional[int]
    bbox: Optional[tuple[float, float, float, float]]
    confidence: Optional[float]
    raw_caption: str
    caption: str
    smi: str
    esmi: str
    cxsmiles: str
    markush: bool
    sru: bool
    groups: Any

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


class MolParser:
    def __init__(self, config: MolParserConfig | str | Path | None = None, **overrides):
        if config is None:
            self.config = MolParserConfig.from_default(**overrides)
        elif isinstance(config, (str, Path)):
            self.config = MolParserConfig.from_yaml(config, **overrides)
        elif overrides:
            self.config = replace(config, **_normalize_overrides(overrides))
        else:
            self.config = config
        self._molparser = None
        self._image_detector = None
        self._pdf_detector = None

    def parse(
        self,
        inputs,
        *,
        rec_only: bool = False,
        pages=None,
        expand_px: int | None = None,
        padding_px: int | None = None,
    ) -> list[MolParserResult]:
        with self._temporary_margins(expand_px=expand_px, padding_px=padding_px):
            image_records: list[tuple[Any, str, int]] = []
            pdf_records: list[tuple[Path, str, int]] = []
            for input_index, item in enumerate(normalize_inputs(inputs)):
                if item.kind == "pdf":
                    if item.path is None:
                        raise ValueError("PDF input requires a file path.")
                    pdf_records.append((item.path, item.source, input_index))
                else:
                    image_records.append((load_image(item), item.source, input_index))
            results: list[MolParserResult] = []
            if image_records:
                results.extend(self._parse_image_batch(image_records, rec_only=rec_only))
            if pdf_records:
                results.extend(self._parse_pdf_batch(pdf_records, pages=pages))
            return sorted(results, key=lambda result: result.input_index)

    def parse_image(
        self,
        image,
        *,
        rec_only: bool = True,
        expand_px: int | None = None,
        padding_px: int | None = None,
    ) -> list[MolParserResult]:
        with self._temporary_margins(expand_px=expand_px, padding_px=padding_px):
            item = normalize_inputs(image)[0]
            return self._parse_image(load_image(item), item.source, 0, rec_only=rec_only)

    def parse_pdf(
        self,
        pdf,
        *,
        pages=None,
        expand_px: int | None = None,
        padding_px: int | None = None,
    ) -> list[MolParserResult]:
        with self._temporary_margins(expand_px=expand_px, padding_px=padding_px):
            item = normalize_inputs(pdf)[0]
            if item.kind != "pdf" or item.path is None:
                raise ValueError("parse_pdf expects a PDF path or URL.")
            return self._parse_pdf_batch([(item.path, item.source, 0)], pages=pages)

    def _token(self):
        return self.config.hf_token if self.config.hf_token is not None else os.environ.get("HF_TOKEN")

    def _cache_dir(self):
        return self.config.cache_dir or None

    @contextmanager
    def _temporary_margins(self, *, expand_px: int | None, padding_px: int | None):
        old_expand = self.config.expand_px
        old_padding = self.config.padding_px
        if expand_px is not None:
            self.config.expand_px = int(expand_px)
        if padding_px is not None:
            self.config.padding_px = int(padding_px)
        try:
            yield
        finally:
            self.config.expand_px = old_expand
            self.config.padding_px = old_padding

    def _recognize(self, images: list[Any]) -> list[str]:
        if self._molparser is None:
            from .runtime import MolParserRecognizer

            path = resolve_model(
                local_path=self.config.molparser_model_path,
                hf_model_id=self.config.molparser_hf_repo,
                modelscope_model_id=self.config.molparser_modelscope_repo,
                cache_dir=self._cache_dir(),
                token=self._token(),
            )
            self._molparser = MolParserRecognizer(
                str(path),
                device=self.config.device,
                token=self._token(),
                max_length=self.config.max_length,
            )
        captions: list[str] = []
        for batch in _batched(images, self.config.molparser_batch_size):
            captions.extend(self._molparser.recognize(batch))
        return captions

    def _detector(self, kind: str):
        from .runtime import MolDetDetector

        if kind == "pdf":
            if self._pdf_detector is None:
                path = self._resolve_detector(self.config.moldet_pdf_model_path, self.config.moldet_pdf_modelname)
                self._pdf_detector = MolDetDetector(
                    str(path), device=self.config.device, imgsz=self.config.pdf_imgsz, conf=self.config.detector_conf
                )
            return self._pdf_detector
        if self._image_detector is None:
            path = self._resolve_detector(self.config.moldet_image_model_path, self.config.moldet_image_modelname)
            self._image_detector = MolDetDetector(
                str(path), device=self.config.device, imgsz=self.config.image_imgsz, conf=self.config.detector_conf
            )
        return self._image_detector

    def _resolve_detector(self, local_path: str, filename: str) -> Path:
        return resolve_model(
            local_path=local_path,
            hf_model_id=self.config.moldet_hf_repo,
            modelscope_model_id=self.config.moldet_modelscope_repo,
            filename=filename,
            cache_dir=self._cache_dir(),
            token=self._token(),
        )

    def _parse_image(self, image, source: str, input_index: int, *, rec_only: bool) -> list[MolParserResult]:
        return self._parse_image_batch([(image, source, input_index)], rec_only=rec_only)

    def _parse_image_batch(self, image_records: list[tuple[Any, str, int]], *, rec_only: bool) -> list[MolParserResult]:
        images = [record[0] for record in image_records]
        sources = [record[1] for record in image_records]
        input_indexes = [record[2] for record in image_records]
        if rec_only:
            padded_images = [self._pad_image(image) for image in images]
            return self._results_many(padded_images, sources, input_indexes, [None] * len(images), [None] * len(images), [None] * len(images))

        detections_by_image = self._detect("image", images)
        crops: list[Any] = []
        crop_sources: list[str] = []
        crop_input_indexes: list[int] = []
        bboxes: list[tuple[float, float, float, float] | None] = []
        confidences: list[float | None] = []
        for image, source, input_index, detections in zip(images, sources, input_indexes, detections_by_image):
            for detection in detections:
                crops.append(self._crop(image, detection))
                crop_sources.append(source)
                crop_input_indexes.append(input_index)
                bboxes.append(detection.bbox)
                confidences.append(detection.confidence)
        return self._results_many(
            crops,
            crop_sources,
            crop_input_indexes,
            [None] * len(crops),
            bboxes,
            confidences,
        )

    def _parse_pdf(self, path: Path, source: str, input_index: int, *, pages=None) -> list[MolParserResult]:
        return self._parse_pdf_batch([(path, source, input_index)], pages=pages)

    def _parse_pdf_batch(self, pdf_records: list[tuple[Path, str, int]], *, pages=None) -> list[MolParserResult]:
        page_records = []
        for path, source, input_index in pdf_records:
            for page in render_pdf(path, dpi=self.config.pdf_dpi, pages=pages):
                page_records.append((page.image, source, input_index, page.page_index))
        if not page_records:
            return []

        page_images = [record[0] for record in page_records]
        detections_by_page = self._detect("pdf", page_images)
        crops: list[Any] = []
        crop_sources: list[str] = []
        crop_input_indexes: list[int] = []
        page_indexes: list[int | None] = []
        bboxes: list[tuple[float, float, float, float] | None] = []
        confidences: list[float | None] = []
        for (image, source, input_index, page_index), detections in zip(page_records, detections_by_page):
            for detection in detections:
                crops.append(self._crop(image, detection))
                crop_sources.append(source)
                crop_input_indexes.append(input_index)
                page_indexes.append(page_index)
                bboxes.append(detection.bbox)
                confidences.append(detection.confidence)
        return self._results_many(crops, crop_sources, crop_input_indexes, page_indexes, bboxes, confidences)

    def _crop(self, image, detection):
        from .runtime import crop_detection

        return crop_detection(image, detection, expand_px=self.config.expand_px, pad_px=self.config.padding_px)

    def _pad_image(self, image):
        from .runtime import pad_image

        return pad_image(image, pad_px=self.config.padding_px)

    def _detect(self, kind: str, images: list[Any]):
        detector = self._detector(kind)
        detections: list[Any] = []
        for batch in _batched(images, self.config.moldet_batch_size):
            detections.extend(detector.detect(batch))
        return detections

    def _results(
        self,
        images: list[Any],
        source: str,
        input_index: int,
        page_indexes: list[int | None],
        bboxes: list[tuple[float, float, float, float] | None],
        confidences: list[float | None],
    ) -> list[MolParserResult]:
        return self._results_many(
            images,
            [source] * len(images),
            [input_index] * len(images),
            page_indexes,
            bboxes,
            confidences,
        )

    def _results_many(
        self,
        images: list[Any],
        sources: list[str],
        input_indexes: list[int],
        page_indexes: list[int | None],
        bboxes: list[tuple[float, float, float, float] | None],
        confidences: list[float | None],
    ) -> list[MolParserResult]:
        if not images:
            return []
        raws = self._recognize(images)
        results: list[MolParserResult] = []
        for raw, source, input_index, page_index, bbox, confidence in zip(raws, sources, input_indexes, page_indexes, bboxes, confidences):
            post = postprocess_caption(raw)
            results.append(
                MolParserResult(
                    source=source,
                    input_index=input_index,
                    page_index=page_index,
                    bbox=bbox,
                    confidence=confidence,
                    raw_caption=str(raw),
                    caption=str(post.get("caption", raw)),
                    smi=str(post.get("smi", "")),
                    esmi=str(post.get("esmi", "")),
                    cxsmiles=str(post.get("cxsmiles", "")),
                    markush=bool(post.get("markush", False)),
                    sru=bool(post.get("sru", False)),
                    groups=post.get("groups", ""),
                )
            )
        return results


def _normalize_overrides(overrides: dict[str, Any]) -> dict[str, Any]:
    aliases = {
        "moldet_hf_model_id": "moldet_hf_repo",
        "moldet_modelscope_model_id": "moldet_modelscope_repo",
        "recognizer_hf_model_id": "molparser_hf_repo",
        "recognizer_modelscope_model_id": "molparser_modelscope_repo",
        "recognizer_model_path": "molparser_model_path",
        "det_batch_size": "moldet_batch_size",
        "ocsr_batch_size": "molparser_batch_size",
    }
    return {aliases.get(key, key): value for key, value in overrides.items()}


def _batched(items: list[Any], batch_size: int):
    size = max(1, int(batch_size or 1))
    for start in range(0, len(items), size):
        yield items[start : start + size]


__all__ = ["MolParser", "MolParserConfig", "MolParserResult"]
