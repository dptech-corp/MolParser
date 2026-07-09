from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
from urllib.parse import urlparse
from urllib.request import urlopen

from PIL import Image


IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".bmp", ".gif", ".webp", ".tif", ".tiff"}


@dataclass
class ParserInput:
    kind: str
    source: str
    path: Path | None = None
    image: Image.Image | None = None


@dataclass
class RenderedPage:
    page_index: int
    image: Image.Image


def is_url(value: str) -> bool:
    return urlparse(value).scheme in {"http", "https"}


def _suffix_from_url(url: str, content_type: str = "") -> str:
    suffix = Path(urlparse(url).path).suffix.lower()
    if suffix:
        return suffix
    if "pdf" in content_type:
        return ".pdf"
    if "png" in content_type:
        return ".png"
    if "jpeg" in content_type or "jpg" in content_type:
        return ".jpg"
    return ".bin"


def download_url(url: str) -> Path:
    with urlopen(url) as response:
        content_type = response.headers.get("content-type", "")
        suffix = _suffix_from_url(url, content_type)
        with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
            tmp.write(response.read())
            return Path(tmp.name)


def _kind_from_path(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".pdf":
        return "pdf"
    if suffix in IMAGE_SUFFIXES:
        return "image"
    raise ValueError(f"Unsupported input suffix: {path.suffix}")


def normalize_inputs(inputs) -> list[ParserInput]:
    if isinstance(inputs, (list, tuple)):
        items: Iterable = inputs
    else:
        items = [inputs]

    normalized: list[ParserInput] = []
    for item in items:
        if isinstance(item, Image.Image):
            normalized.append(ParserInput(kind="image", source="PIL.Image", image=item.convert("RGB")))
            continue
        if isinstance(item, (str, Path)):
            raw = str(item)
            path = download_url(raw) if is_url(raw) else Path(raw).expanduser()
            normalized.append(ParserInput(kind=_kind_from_path(path), source=raw, path=path))
            continue
        raise TypeError(f"Unsupported input type: {type(item)!r}")
    return normalized


def load_image(item: ParserInput) -> Image.Image:
    if item.image is not None:
        return item.image.convert("RGB")
    if item.path is None:
        raise ValueError("Image input has no image or path.")
    return Image.open(item.path).convert("RGB")


def render_pdf(path: str | Path, *, dpi: int = 300, pages: Iterable[int] | None = None) -> list[RenderedPage]:
    try:
        import pypdfium2 as pdfium
    except ImportError as exc:
        raise ImportError("PDF parsing requires pypdfium2. Install molparser[pdf] or molparser[inference].") from exc

    document = pdfium.PdfDocument(str(path))
    page_indexes = list(pages) if pages is not None else list(range(len(document)))
    scale = float(dpi) / 72.0
    return [
        RenderedPage(page_index=page_index, image=document[page_index].render(scale=scale).to_pil().convert("RGB"))
        for page_index in page_indexes
    ]
