"""UCSC conservation figure validation, presets, storage, and safe browser URLs."""

from __future__ import annotations

import hashlib
import html
import io
import logging
import re
import tempfile
import threading
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
UCSC_SECTION_1B_TRACK_PRESET_ID = "ucsc_section_1b_comprehensive_v1"
UCSC_SECTION_1B_TRACK_PRESET_VERSION = 1
UCSC_SECTION_1B_TRACK_PRESET: dict[str, str] = {
    # Reset browser state so rendering does not depend on cookies/session defaults.
    "hideTracks": "1",
    "pix": "1400",
    "knownGene": "pack",
    "ncbiRefSeqCurated": "pack",
    "mane": "pack",
    "omimGene2": "pack",
    "snp155": "pack",
    "gtexGeneV8": "pack",
    "wgEncodeReg": "pack",
    "wgEncodeRegMarkH3k27ac": "pack",
    "cons100way": "full",
    "multiz100way": "pack",
    "rmsk": "pack",
    "guidelines": "off",
    "textSize": "10",
}


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


@dataclass
class StagedFigure:
    """Temporary or reusable managed figure path information."""

    final_absolute_path: Path
    relative_path: str
    sha256: str
    media_type: str
    width: int
    height: int
    byte_size: int
    extension: str
    temp_path: Path | None
    existed_already: bool


def sha256_hex(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


_API_KEY_QUERY_RE = re.compile(
    r"([?&]api[_-]?key=)[^&\s\"']+",
    re.IGNORECASE,
)
_API_KEY_KV_RE = re.compile(
    r"(api[_-]?key\s*[=:]\s*)([^\s,\"']+)",
    re.IGNORECASE,
)


def redact_api_key(text: str | None) -> str:
    """Remove apiKey query values case-insensitively from text/URLs/errors."""
    if not text:
        return ""
    out = _API_KEY_QUERY_RE.sub(r"\1REDACTED", str(text))
    out = _API_KEY_KV_RE.sub(r"\1REDACTED", out)
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


class ApiKeyRedactionFilter(logging.Filter):
    """Defense-in-depth filter; prefer the LogRecordFactory for hierarchy coverage."""

    def filter(self, record: logging.LogRecord) -> bool:  # noqa: A003
        _sanitize_log_record(record)
        return True


def _sanitize_log_record(record: logging.LogRecord) -> logging.LogRecord:
    """Redact apiKey material on a LogRecord in place."""
    try:
        try:
            rendered = record.getMessage()
        except Exception:  # noqa: BLE001 — malformed msg/args
            rendered = f"{record.msg!s} {record.args!s}"
        record.msg = redact_api_key(rendered)
        record.args = ()
        if record.exc_info and record.exc_info[1] is not None:
            exc = record.exc_info[1]
            try:
                exc.args = tuple(
                    redact_api_key(a) if isinstance(a, str) else a for a in exc.args
                )
            except Exception:  # noqa: BLE001
                pass
    except Exception:  # noqa: BLE001
        pass
    return record


_REDACTION_LOCK = threading.Lock()
_REDACTION_INSTALLED = False
_previous_log_record_factory = logging.getLogRecordFactory()


def _redacting_record_factory(*args: Any, **kwargs: Any) -> logging.LogRecord:
    """Sanitize every LogRecord at creation, including httpcore.* descendants."""
    record = _previous_log_record_factory(*args, **kwargs)
    return _sanitize_log_record(record)


def install_ucsc_api_key_log_redaction() -> None:
    """Install permanent apiKey redaction for all logging hierarchies.

    A chained :func:`logging.setLogRecordFactory` sanitizes every new record
    regardless of logger name (including ``httpcore.connection`` /
    ``httpcore.http11``). Logger-level :class:`ApiKeyRedactionFilter` instances
    remain as defense in depth.
    """
    global _REDACTION_INSTALLED, _previous_log_record_factory
    with _REDACTION_LOCK:
        if _REDACTION_INSTALLED:
            return
        _previous_log_record_factory = logging.getLogRecordFactory()
        logging.setLogRecordFactory(_redacting_record_factory)
        filt = ApiKeyRedactionFilter()
        for name in (
            "httpx",
            "httpcore",
            "gene_dossier.tools.ucsc",
            "gene_dossier.workflow",
        ):
            logging.getLogger(name).addFilter(filt)
        logging.getLogger().addFilter(filt)
        _REDACTION_INSTALLED = True


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


def validate_image_bytes(
    content: bytes,
    *,
    min_width: int = 8,
    min_height: int = 8,
) -> tuple[ValidatedImage | None, FigureValidationError | None]:
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

    if width < min_width or height < min_height:
        return None, FigureValidationError(
            "trivial_figure",
            f"Image dimensions {width}x{height} below minimum {min_width}x{min_height}",
        )

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


def validate_live_render_image_bytes(
    content: bytes,
    *,
    requested_pix: int = 1400,
) -> tuple[ValidatedImage | None, FigureValidationError | None]:
    """Validate a live hgRenderTracks image against a stronger size threshold."""
    # Require a meaningful fraction of the requested pixel width; reject logos.
    min_width = max(200, int(requested_pix * 0.5))
    min_height = 80
    return validate_image_bytes(content, min_width=min_width, min_height=min_height)


def stage_figure_tempfile(
    *,
    dossier_run_id: str,
    content: bytes,
    extension: str = "png",
    settings: Settings | None = None,
) -> StagedFigure:
    """Stage figure bytes to a temp file in the final managed directory.

    Does not move the file into place or commit any DB state. Callers are responsible
    for moving `temp_path` to `final_absolute_path`, verifying the checksum, and
    deleting temporary/new files on failure.
    """
    cfg = settings or get_settings()
    validated, err = validate_image_bytes(content)
    if err or validated is None:
        raise ValueError(err.message if err else "invalid figure")

    digest = validated.sha256
    rel_dir = Path(dossier_run_id) / "ucsc" / "figures"
    ext = extension.lstrip(".") or "png"
    final_rel = rel_dir / f"{digest}.{ext}"
    final_abs = (cfg.raw_data_path / final_rel).resolve()
    final_abs.parent.mkdir(parents=True, exist_ok=True)

    if final_abs.is_file():
        if sha256_hex(final_abs.read_bytes()) != digest:
            raise ValueError(f"existing figure at {final_abs} does not match checksum")
        return StagedFigure(
            final_absolute_path=final_abs,
            relative_path=final_rel.as_posix(),
            sha256=digest,
            media_type=validated.media_type,
            width=validated.width,
            height=validated.height,
            byte_size=validated.byte_size,
            extension=ext,
            temp_path=None,
            existed_already=True,
        )

    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{digest}.",
        suffix=f".{ext}.tmp",
        dir=str(final_abs.parent),
    )
    tmp_path = Path(tmp_name)
    try:
        with open(fd, "wb") as handle:
            handle.write(validated.content)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise

    return StagedFigure(
        final_absolute_path=final_abs,
        relative_path=final_rel.as_posix(),
        sha256=digest,
        media_type=validated.media_type,
        width=validated.width,
        height=validated.height,
        byte_size=validated.byte_size,
        extension=ext,
        temp_path=tmp_path,
        existed_already=False,
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
    staged = stage_figure_tempfile(
        dossier_run_id=dossier_run_id,
        content=content,
        extension=extension,
        settings=settings,
    )
    if staged.temp_path is not None:
        staged.temp_path.replace(staged.final_absolute_path)
        if sha256_hex(staged.final_absolute_path.read_bytes()) != staged.sha256:
            staged.final_absolute_path.unlink(missing_ok=True)
            raise ValueError("figure checksum mismatch after final move")
    return staged.final_absolute_path, staged.relative_path, staged.sha256


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
    """Locate an approved UCSC generated-track image inside an HTML wrapper.

    Prefers ``/trash/`` image assets over logos or unrelated ``<img>`` tags.
    HTML-decodes ``src`` values before parsing.
    """
    candidates: list[str] = []
    preferred: list[str] = []
    for match in re.finditer(r"""src=["']([^"']+)["']""", html_text, re.IGNORECASE):
        candidate = html.unescape(match.group(1).strip())
        try:
            parsed = urlparse(candidate)
        except Exception:  # noqa: BLE001
            continue
        if parsed.scheme and parsed.scheme not in {"http", "https"}:
            continue
        host = parsed.hostname
        path = (parsed.path or "").lower()
        # Skip obvious UI chrome.
        if any(token in path for token in ("logo", "favicon", "icon", "sprite")):
            continue
        resolved: str | None = None
        if host is None and candidate.startswith("/"):
            resolved = f"https://genome.ucsc.edu{candidate}"
        elif host in _ALLOWED_UCSC_HOSTS:
            if parsed.scheme == "http":
                resolved = urlunparse(
                    ("https", parsed.netloc, parsed.path, "", parsed.query, "")
                )
            else:
                resolved = candidate
        if not resolved:
            continue
        candidates.append(resolved)
        if "/trash/" in path or path.endswith((".png", ".jpg", ".jpeg", ".gif")):
            preferred.append(resolved)
    if preferred:
        return preferred[0]
    return candidates[0] if candidates else None


def split_url_for_provenance(url: str) -> tuple[str, dict[str, str]]:
    """Return credential-free base URL and sanitized query params for persistence."""
    parsed = urlparse(url)
    base = f"{parsed.scheme}://{parsed.netloc}{parsed.path}" if parsed.scheme and parsed.netloc else url
    params = sanitize_params(dict(parse_qsl(parsed.query, keep_blank_values=True)))
    return base, {str(k): str(v) for k, v in params.items()}


__all__ = [
    "ValidatedImage",
    "FigureValidationError",
    "StagedFigure",
    "UCSC_SECTION_1B_TRACK_PRESET_ID",
    "UCSC_SECTION_1B_TRACK_PRESET_VERSION",
    "UCSC_SECTION_1B_TRACK_PRESET",
    "ApiKeyRedactionFilter",
    "sha256_hex",
    "redact_api_key",
    "sanitize_params",
    "install_ucsc_api_key_log_redaction",
    "validate_image_bytes",
    "validate_live_render_image_bytes",
    "stage_figure_tempfile",
    "relative_to_artifact_root",
    "resolve_artifact_path",
    "stage_and_commit_figure",
    "build_safe_hgtracks_url",
    "is_safe_ucsc_browser_url",
    "extract_ucsc_image_url_from_html",
    "split_url_for_provenance",
]
