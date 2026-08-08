"""Unit tests for BioGRID Section-5b client helpers."""

from __future__ import annotations

from unittest.mock import patch

from gene_dossier.config import Settings
from gene_dossier.tools import biogrid as bg
from gene_dossier.workflow import WorkflowTransientContext


def test_fetch_version_plain_text_no_json_parse() -> None:
    settings = Settings(biogrid_accesskey="SECRETKEY")

    class _Resp:
        status_code = 200
        is_success = True
        url = "https://webservice.thebiogrid.org/version"
        history: list = []
        headers = {"content-type": "text/plain"}
        content = b"5.0.259\n"

        def json(self):
            raise AssertionError("version must not be JSON-parsed")

    class _Client:
        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def get(self, url, params=None):
            assert "version" in url
            return _Resp()

    with patch("gene_dossier.tools.biogrid.httpx.Client", _Client):
        result = bg.fetch_version(gene_symbol="SREBF2", settings=settings)
    assert result.success is True
    data = bg.unwrap_biogrid_payload(result.data)
    assert data["version"] == "5.0.259"
    assert result.request_params.get("accesskey") == "***"


def test_section_page_uses_include_self_and_cross() -> None:
    captured = {}

    class _Resp:
        status_code = 200
        is_success = True
        url = "https://webservice.thebiogrid.org/interactions/"
        history: list = []
        headers = {"content-type": "application/json"}
        content = b"{}"

        def json(self):
            return {}

    class _Client:
        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def get(self, url, params=None):
            captured["params"] = dict(params or {})
            return _Resp()

    with patch("gene_dossier.tools.biogrid.httpx.Client", _Client):
        result = bg.fetch_interactions_page(
            "SREBF2",
            start=0,
            settings=Settings(biogrid_accesskey="SECRETKEY"),
        )
    assert result.success is True
    assert captured["params"]["selfInteractionsExcluded"] == "false"
    assert captured["params"]["interSpeciesExcluded"] == "false"
    assert captured["params"]["start"] == "0"


def test_pagination_fetches_second_page() -> None:
    calls: list[int] = []

    def _payload(start: int, n: int):
        return {
            str(i): {
                "BIOGRID_INTERACTION_ID": start + i,
                "BIOGRID_ID_A": 1,
                "BIOGRID_ID_B": 2 + i,
                "OFFICIAL_SYMBOL_A": "GENE",
                "OFFICIAL_SYMBOL_B": f"P{i}",
                "ENTREZ_GENE_A": "10",
                "ENTREZ_GENE_B": str(100 + i),
                "ORGANISM_A": 9606,
                "ORGANISM_B": 9606,
                "EXPERIMENTAL_SYSTEM_TYPE": "physical",
            }
            for i in range(n)
        }

    class _Resp:
        def __init__(self, payload):
            import json as _json

            self.status_code = 200
            self.is_success = True
            self.url = "https://webservice.thebiogrid.org/interactions/"
            self.history = []
            self.headers = {"content-type": "application/json"}
            self._payload = payload
            self.content = _json.dumps(payload).encode("utf-8")

        def json(self):
            return self._payload

    class _Client:
        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def get(self, url, params=None):
            start = int((params or {}).get("start") or 0)
            calls.append(start)
            n = 3 if start == 0 else 1
            return _Resp(_payload(start, n))

    with patch("gene_dossier.tools.biogrid.httpx.Client", _Client):
        pages, rows, errors = bg.fetch_all_interactions_section_5b(
            "GENE",
            max_results=3,
            settings=Settings(biogrid_accesskey="SECRETKEY"),
        )
    assert errors == []
    assert calls == [0, 3]
    assert len(pages) == 2
    assert len(rows) == 4
    assert rows[0]["_page_start"] == 0
    assert rows[3]["_page_start"] == 3


def test_request_cache_identity_reuses_page() -> None:
    hits = {"n": 0}

    class _Resp:
        status_code = 200
        is_success = True
        url = "https://webservice.thebiogrid.org/interactions/"
        history: list = []
        headers = {"content-type": "application/json"}
        content = b"{}"

        def json(self):
            return {}

    class _Client:
        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def get(self, url, params=None):
            hits["n"] += 1
            return _Resp()

    transient = WorkflowTransientContext()
    with patch("gene_dossier.tools.biogrid.httpx.Client", _Client):
        a = bg.fetch_interactions_page(
            "GENE", start=0, settings=Settings(biogrid_accesskey="SECRETKEY"), transient=transient
        )
        b = bg.fetch_interactions_page(
            "GENE", start=0, settings=Settings(biogrid_accesskey="SECRETKEY"), transient=transient
        )
    assert a.success and b.success
    assert hits["n"] == 1
