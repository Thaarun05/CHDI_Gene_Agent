"""UCSC conservation figure validation, portable storage, and safe browser URLs."""

from __future__ import annotations

import hashlib
import io
import re
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import parse_qsl, urlencode, urlparse, urlunparse

from gene_dossier.config import PROJECT_ROOT, Settings, get_settings

_PNG_MAGIC = b"\x89PNG\r\n\x1a\n"
_JPEG_MAGIC = b"\xff\xd8\xff"
_GIF_MAGIC = b"GIF8"

_BLOCKED_HTML_MARKERS = (
    "cf-turnstile",
    "turnstile",
    "protecting itself from bots",
    "cloudflare",
    "captcha",
    "challenge-platform",
    "ucsc genome browser is temporarily unavailable",
    "hgsid",
    "hgtracks.js",
)

_ALLOWED_UCSC_HOSTS = frozenset({"genome.ucsc.edu", "hgdownload.soe.ucsc.edu"})


@dataclass
class ValidatedImage:
    content: bytes
    media_type: str
    width: int
    height: int
    sha256: str
    byte_size: int


@dataclass
class FigureValidationError:
    code: str
    message: str


def sha256_hex(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def redact_api_key(text: str | None) -> str:
    """Remove apiKey query values case-insensitively from text/URLs/errors."""
    if not text:
        return ""
    # query param forms
    out = re.sub(
        r"([?&]api[_-]?key=)[^&\s\"']+",
        r"\1REDACTED",
        str(text),
        flags=re.IGNORECASE,
    )
    return out


def sanitize_params(params: dict[str, Any] | None) -> dict[str, Any]:
    """Return a copy of params with apiKey removed (any case)."""
    if not params:
        return {}
    out: dict[str, Any] = {}
    for key, value in params.items():
        if str(key).lower().replace("-", "_") in {"apikey", "api_key"}:
            continue
        out[key] = value
    return out


def _image_size_png(content: bytes) -> tuple[int, int] | None:
    if len(content) < 24 or not content.startswith(_PNG_MAGIC):
        return None
    import struct

    width, height = struct.unpack(">II", content[16:24])
    return int(width), int(height)


def _image_size_via_pillow(content: bytes) -> tuple[str, int, int] | None:
    try:
        from PIL import Image
    except ImportError:
        return None
    try:
        with Image.open(io.BytesIO(content)) as img:
            img.load()
            fmt = (img.format or "").upper()
            media = {
                "PNG": "image/png",
                "JPEG": "image/jpeg",
                "GIF": "image/gif",
            }.get(fmt)
            if not media:
                return None
            w, h = img.size
            return media, int(w), int(h)
    except Exception:  # noqa: BLE001
        return None


def validate_image_bytes(content: bytes) -> tuple[ValidatedImage | None, FigureValidationError | None]:
    """Validate PNG/JPEG/GIF bytes; reject HTML/CAPTCHA/bootstrap pages."""
    if not content:
        return None, FigureValidationError("empty_figure", "Figure body is empty")
    head = content[:4096]
    # HTML / challenge detection
    lower = head.lower()
    if b"<html" in lower or b"<!doctype html" in lower or b"<script" in lower:
        text = head.decode("utf-8", errors="replace").lower()
        for marker in _BLOCKED_HTML_MARKERS:
            if marker in text:
                return None, FigureValidationError(
                    "blocked_browser_render",
                    f"Rejected HTML response containing {marker!r}",
                )
        return None, FigureValidationError(
            "html_figure_wrapper",
            "Response is HTML rather than raw image bytes",
        )

    media: str | None = None
    width = height = 0
    if content.startswith(_PNG_MAGIC):
        media = "image/png"
        size = _image_size_png(content)
        if size:
            width, height = size
    elif content.startswith(_JPEG_MAGIC):
        media = "image/jpeg"
    elif content.startswith(_GIF_MAGIC):
        media = "image/gif"

    pillow = _image_size_via_pillow(content)
    if pillow:
        media, width, height = pillow
    elif media == "image/png" and width and height:
        pass
    elif media is None:
        return None, FigureValidationError("invalid_figure", "Unsupported or undecodable image")

    if width < 8 or height < 8:
        return None, FigureValidationError("trivial_figure", "Image dimensions are trivial")

    digest = sha256_hex(content)
    return (
        ValidatedImage(
            content=content,
            media_type=media or "image/png",
            width=width,
            height=height,
            sha256=digest,
            byte_size=len(content),
        ),
        None,
    )


def relative_to_artifact_root(path: Path, *, root: Path | None = None) -> str:
    """Return a portable relative path under the artifact root."""
    settings = get_settings()
    base = (root or settings.raw_data_path).resolve()
    resolved = path.resolve()
    try:
        rel = resolved.relative_to(base)
    except ValueError:
        # Fall back to project-root relative when possible.
        try:
            rel = resolved.relative_to(PROJECT_ROOT.resolve())
        except ValueError as exc:
            raise ValueError(f"figure path escapes artifact root: {path}") from exc
    return rel.as_posix()


def resolve_artifact_path(relative_or_absolute: str, *, root: Path | None = None) -> Path:
    """Resolve a stored figure path against the configured artifact root."""
    settings = get_settings()
    base = (root or settings.raw_data_path).resolve()
    raw = Path(relative_or_absolute)
    if raw.is_absolute():
        # Legacy absolute paths: accept only if under artifact root or project root.
        resolved = raw.resolve()
        for allowed in (base, PROJECT_ROOT.resolve()):
            try:
                resolved.relative_to(allowed)
                return resolved
            except ValueError:
                continue
        raise ValueError("absolute figure path is outside managed storage")
    return (base / raw).resolve()


def stage_and_commit_figure(
    *,
    dossier_run_id: str,
    content: bytes,
    extension: str = "png",
    settings: Settings | None = None,
) -> tuple[Path, str, str]:
    """Copy figure bytes via temp path into managed storage.

    Returns ``(final_absolute_path, relative_path, sha256)``.
    Reuses an existing matching checksum file when present.
    """
    cfg = settings or get_settings()
    validated, err = validate_image_bytes(content)
    if err or validated is None:
        raise ValueError(err.message if err else "invalid figure")

    digest = validated.sha256
    rel_dir = Path(dossier_run_id) / "ucsc" / "figures"
    final_rel = rel_dir / f"{digest}.{extension.lstrip('.')}"
    final_abs = (cfg.raw_data_path / final_rel).resolve()
    final_abs.parent.mkdir(parents=True, exist_ok=True)

    if final_abs.is_file() and sha256_hex(final_abs.read_bytes()) == digest:
        return final_abs, final_rel.as_posix(), digest

    # Stage to temporary file in the same directory for atomic replace.
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{digest}.",
        suffix=f".{extension.lstrip('.')}.tmp",
        dir=str(final_abs.parent),
    )
    tmp_path = Path(tmp_name)
    try:
        with open(fd, "wb") as handle:
            handle.write(validated.content)
        tmp_path.replace(final_abs)
    except Exception:
        if tmp_path.exists():
            tmp_path.unlink(missing_ok=True)
        raise
    return final_abs, final_rel.as_posix(), digest


def build_safe_hgtracks_url(
    *,
    genome: str,
    display_position: str,
    transcript_id: str | None,
) -> str | None:
    """Build a credential-free UCSC hgTracks URL with an allowlisted query."""
    if not genome or not display_position:
        return None
    if not re.match(r"^chr[\w.]+:\d+-\d+$", display_position, re.IGNORECASE):
        return None
    if transcript_id is not None and not _ENST_OR_SAFE.match(transcript_id):
        return None
    params: dict[str, str] = {
        "db": genome,
        "position": display_position,
    }
    if transcript_id:
        params["hgFind.matches"] = transcript_id
    query = urlencode(params)
    return f"https://genome.ucsc.edu/cgi-bin/hgTracks?{query}"


_ENST_OR_SAFE = re.compile(r"^[A-Za-z0-9._-]+$")


def is_safe_ucsc_browser_url(url: str) -> bool:
    """Validate a UCSC browser URL against the allowlist."""
    try:
        parsed = urlparse(url)
    except Exception:  # noqa: BLE001
        return False
    if parsed.scheme != "https":
        return False
    if parsed.hostname not in {"genome.ucsc.edu"}:
        return False
    if parsed.path != "/cgi-bin/hgTracks":
        return False
    allowed = {"db", "position", "hgfind.matches"}
    for key, _value in parse_qsl(parsed.query, keep_blank_values=True):
        if key.lower() == "apikey" or key.lower().replace("-", "_") == "api_key":
            return False
        if key.lower() not in allowed:
            return False
    return True


def extract_ucsc_image_url_from_html(html_text: str) -> str | None:
    """Locate an approved-host image reference inside an hgRenderTracks HTML wrapper."""
    # src="..." or src='...'
    for match in re.finditer(r"""src=["']([^"']+)["']""", html_text, re.IGNORECASE):
        candidate = match.group(1)
        try:
            parsed = urlparse(candidate)
        except Exception:  # noqa: BLE001
            continue
        if parsed.scheme and parsed.scheme not in {"http", "https"}:
            continue
        host = parsed.hostname
        if host is None and candidate.startswith("/"):
            return f"https://genome.ucsc.edu{candidate}"
        if host in _ALLOWED_UCSC_HOSTS:
            if parsed.scheme == "http":
                return urlunparse(("https", parsed.netloc, parsed.path, "", parsed.query, ""))
            return candidate
    return None


__all__ = [
    "ValidatedImage",
    "FigureValidationError",
    "sha256_hex",
    "redact_api_key",
    "sanitize_params",
    "validate_image_bytes",
    "relative_to_artifact_root",
    "resolve_artifact_path",
    "stage_and_commit_figure",
    "build_safe_hgtracks_url",
    "is_safe_ucsc_browser_url",
    "extract_ucsc_image_url_from_html",
]
