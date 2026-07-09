from __future__ import annotations

import contextlib
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

from PIL import Image


@dataclass
class ModelResolveError(RuntimeError):
    message: str

    def __str__(self) -> str:
        return self.message


@dataclass
class Detection:
    bbox: tuple[float, float, float, float]
    confidence: float


def resolve_model(
    *,
    local_path: str | Path | None = None,
    hf_model_id: str = "",
    modelscope_model_id: str = "",
    filename: str | None = None,
    cache_dir: str | Path | None = None,
    token: str | bool | None = None,
) -> Path:
    existing = _as_existing_path(local_path)
    if existing is not None:
        return _find_file(existing, filename) if filename and existing.is_dir() else existing

    errors: list[str] = []
    if hf_model_id:
        try:
            if filename:
                from huggingface_hub import hf_hub_download

                return Path(
                    hf_hub_download(
                        repo_id=hf_model_id,
                        filename=filename,
                        cache_dir=str(cache_dir) if cache_dir else None,
                        token=token,
                    )
                )
            from huggingface_hub import snapshot_download

            return Path(snapshot_download(repo_id=hf_model_id, cache_dir=str(cache_dir) if cache_dir else None, token=token))
        except Exception as exc:
            errors.append(f"Hugging Face failed for {hf_model_id}: {exc!r}")

    if modelscope_model_id:
        try:
            try:
                from modelscope import snapshot_download
            except ImportError:
                from modelscope.hub.snapshot_download import snapshot_download

            root = Path(snapshot_download(modelscope_model_id, cache_dir=str(cache_dir) if cache_dir else None))
            return _find_file(root, filename) if filename else root
        except Exception as exc:
            errors.append(f"ModelScope failed for {modelscope_model_id}: {exc!r}")

    details = " ".join(errors) if errors else "No model id was configured."
    raise ModelResolveError(f"Unable to resolve model. {details} Set a local path or configure HF/ModelScope.")


class MolDetDetector:
    def __init__(self, model_path: str, *, device: str = "auto", imgsz: int = 640, conf: float = 0.5):
        try:
            from ultralytics import YOLO
        except ImportError as exc:
            raise ImportError("MolDet requires ultralytics. Install molparser[moldet] or molparser[inference].") from exc

        self.model = YOLO(str(model_path))
        self.device = None if device == "auto" else device
        self.imgsz = int(imgsz)
        self.conf = float(conf)

    def detect(self, images: Sequence[Image.Image]) -> list[list[Detection]]:
        results = self.model.predict(list(images), imgsz=self.imgsz, conf=self.conf, device=self.device, verbose=False)
        batches: list[list[Detection]] = []
        for result in results:
            detections: list[Detection] = []
            boxes = getattr(result, "boxes", None)
            if boxes is not None:
                for box in boxes:
                    xyxy = tuple(float(v) for v in box.xyxy[0].tolist())
                    detections.append(Detection(bbox=xyxy, confidence=float(box.conf)))
            batches.append(detections)
        return batches


class MolParserRecognizer:
    def __init__(
        self,
        model_path: str,
        *,
        device: str = "auto",
        token=None,
        max_length: int = 256,
    ):
        try:
            import torch
            from transformers import AutoModelForImageTextToText, AutoProcessor
        except ImportError as exc:
            raise ImportError("MolParser recognition requires torch and transformers. Install molparser[inference].") from exc

        self.torch = torch
        self.device = torch.device("cuda" if device == "auto" and torch.cuda.is_available() else ("cpu" if device == "auto" else device))
        self.dtype = torch.float16 if self.device.type == "cuda" else torch.float32
        self.max_length = int(max_length)
        self.processor = AutoProcessor.from_pretrained(model_path, trust_remote_code=True, token=token)
        self.model = AutoModelForImageTextToText.from_pretrained(
            model_path,
            trust_remote_code=True,
            token=token,
            dtype=self.dtype,
        ).to(self.device)
        self.model.eval()

    def recognize(self, images: Sequence[Image.Image]) -> list[str]:
        inputs = self.processor(images=[image.convert("RGB") for image in images], return_tensors="pt")
        pixel_values = inputs["pixel_values"].to(device=self.device, dtype=self.dtype)
        autocast = (
            self.torch.autocast(device_type="cuda", dtype=self.dtype)
            if self.device.type == "cuda"
            else contextlib.nullcontext()
        )
        with self.torch.inference_mode(), autocast:
            outputs = self.model.generate(
                pixel_values=pixel_values,
                max_length=self.max_length,
                num_beams=1,
                do_sample=False,
                use_cache=True,
            )
        return self.processor.batch_decode(outputs.detach().cpu().tolist(), skip_special_tokens=True)


def crop_detection(image: Image.Image, detection: Detection, *, expand_px: int = 0, pad_px: int = 0) -> Image.Image:
    width, height = image.size
    x1, y1, x2, y2 = detection.bbox
    box = (
        max(0, int(round(x1)) - expand_px),
        max(0, int(round(y1)) - expand_px),
        min(width, int(round(x2)) + expand_px),
        min(height, int(round(y2)) + expand_px),
    )
    crop = image.crop(box).convert("RGB")
    return pad_image(crop, pad_px=pad_px)


def pad_image(image: Image.Image, *, pad_px: int = 0) -> Image.Image:
    image = image.convert("RGB")
    if pad_px <= 0:
        return image
    padded = Image.new("RGB", (image.width + pad_px * 2, image.height + pad_px * 2), "white")
    padded.paste(image, (pad_px, pad_px))
    return padded


def _as_existing_path(path: str | Path | None) -> Optional[Path]:
    if not path:
        return None
    resolved = Path(path).expanduser()
    return resolved if resolved.exists() else None


def _find_file(root: Path, filename: str | None) -> Path:
    if not filename:
        return root
    candidate = root / filename
    if candidate.exists():
        return candidate
    matches = list(root.rglob(filename))
    if matches:
        return matches[0]
    raise ModelResolveError(f"{filename} not found under model directory: {root}")


MobileRecognizer = MolParserRecognizer


__all__ = [
    "Detection",
    "MobileRecognizer",
    "ModelResolveError",
    "MolDetDetector",
    "MolParserRecognizer",
    "crop_detection",
    "pad_image",
    "resolve_model",
]
