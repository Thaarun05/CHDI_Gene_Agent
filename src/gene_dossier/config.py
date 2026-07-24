"""Central configuration for the Gene Dossier Platform.

Settings are loaded from environment variables and an optional ``.env`` file (see
``.env.example``). All API keys are optional: missing keys degrade gracefully and the
retrieval / normalization / reporting pipeline still runs. LLM usage is optional too.

Access settings via :func:`get_settings`, which returns a cached singleton.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from pydantic import Field, computed_field
from pydantic_settings import BaseSettings, SettingsConfigDict

# Repo root = two levels up from this file (src/gene_dossier/config.py -> repo root).
PROJECT_ROOT: Path = Path(__file__).resolve().parents[2]


def _resolve(path: str | Path) -> Path:
    """Resolve ``path`` against the project root if it is relative."""
    p = Path(path)
    return p if p.is_absolute() else (PROJECT_ROOT / p)


class Settings(BaseSettings):
    """Environment-backed configuration.

    Field names map to upper-cased environment variables (e.g. ``ncbi_api_key`` ->
    ``NCBI_API_KEY``). Unknown environment variables are ignored.
    """

    model_config = SettingsConfigDict(
        env_file=".env",
        env_file_encoding="utf-8",
        case_sensitive=False,
        extra="ignore",
    )

    # --- API keys (all optional) ---
    ncbi_api_key: str | None = None
    biogrid_accesskey: str | None = None
    omim_api_key: str | None = None
    serpapi_api_key: str | None = None

    # --- LLM providers (optional) ---
    openai_api_key: str | None = None
    openai_base_url: str | None = None
    anthropic_api_key: str | None = None
    nvidia_nim_api_key: str | None = None
    nvidia_nim_base_url: str | None = None
    nvidia_nim_model: str | None = None
    default_llm_model: str | None = None
    # openai | anthropic | nvidia_nim | unset (auto: openai → nim → anthropic)
    default_llm_provider: str | None = None

    # --- Storage / database ---
    database_url: str = "sqlite:///data/gene_dossier.db"
    raw_data_dir: Path = Field(default=Path("data/raw"))
    output_dir: Path = Field(default=Path("data/outputs"))
    index_dir: Path = Field(default=Path("data/indexes"))

    # HTTP defaults shared by API clients.
    http_timeout_seconds: float = 30.0
    caller_identity: str = "gene_dossier_platform"

    # --- Resolved absolute paths ---
    @computed_field  # type: ignore[prop-decorator]
    @property
    def raw_data_path(self) -> Path:
        """Absolute path to the raw artifact directory."""
        return _resolve(self.raw_data_dir)

    @computed_field  # type: ignore[prop-decorator]
    @property
    def output_path(self) -> Path:
        """Absolute path to the report output directory."""
        return _resolve(self.output_dir)

    @computed_field  # type: ignore[prop-decorator]
    @property
    def index_path(self) -> Path:
        """Absolute path to the vector/index directory (reserved for future RAG)."""
        return _resolve(self.index_dir)

    # --- Helpers ---
    def has_key(self, name: str) -> bool:
        """Return True if the named key setting is present and non-empty.

        ``name`` may be a field name (``"ncbi_api_key"``) or an env-style name
        (``"NCBI_API_KEY"``).
        """
        value = getattr(self, name.lower(), None)
        return bool(value and str(value).strip())

    def has_llm(self) -> bool:
        """Return True if any LLM provider key is configured."""
        return (
            self.has_key("openai_api_key")
            or self.has_key("anthropic_api_key")
            or self.has_key("nvidia_nim_api_key")
        )

    def ensure_dirs(self) -> None:
        """Create the raw, output, and index directories if they do not exist."""
        for path in (self.raw_data_path, self.output_path, self.index_path):
            path.mkdir(parents=True, exist_ok=True)


@lru_cache(maxsize=1)
def get_settings() -> Settings:
    """Return the cached :class:`Settings` singleton."""
    return Settings()
