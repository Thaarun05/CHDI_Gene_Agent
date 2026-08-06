#!/usr/bin/env python3
"""Focused Section 3a render smoke wrapper."""

from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
_SCRIPTS = _REPO_ROOT / "scripts"
for path in (str(_SRC), str(_SCRIPTS)):
    if path not in sys.path:
        sys.path.insert(0, path)

from run_section_bundle import main as bundle_main  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    args = list(argv if argv is not None else sys.argv[1:])
    if "--sections" not in args:
        args = ["--sections", "3a", *args]
    return bundle_main(args)


if __name__ == "__main__":
    raise SystemExit(main())
