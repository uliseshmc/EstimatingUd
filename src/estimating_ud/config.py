"""Where the data and the results live.

The alignment download root is read from the ``DATA_HUMCHIMPORANGOR114`` environment
variable (or the ``--data-root`` CLI option). Results default to ``<repo>/results``,
located from this file when the package is installed in editable mode, and can be
overridden with ``ESTIMATING_UD_RESULTS``.
"""

from __future__ import annotations

import os
from pathlib import Path

DATA_ENV = "DATA_HUMCHIMPORANGOR114"
RESULTS_ENV = "ESTIMATING_UD_RESULTS"


def repo_root() -> Path | None:
    """Return the repository root (the directory holding pyproject.toml), or None."""
    for parent in Path(__file__).resolve().parents:
        if (parent / "pyproject.toml").exists():
            return parent
    return None


def data_root(override: str | os.PathLike | None = None) -> Path:
    """Return the alignment download root as a Path."""
    value = override or os.environ.get(DATA_ENV)
    if not value:
        msg = f"set the {DATA_ENV} environment variable to the alignment download location"
        raise RuntimeError(msg)
    return Path(value).expanduser()


def results_dir(subdir: str | None = None, override: str | os.PathLike | None = None) -> Path:
    """Return (and create) the results directory, optionally a sub-directory of it."""
    if override:
        base = Path(override).expanduser()
    elif os.environ.get(RESULTS_ENV):
        base = Path(os.environ[RESULTS_ENV]).expanduser()
    else:
        root = repo_root()
        base = (root / "results") if root else (Path.cwd() / "results")
    path = base / subdir if subdir else base
    path.mkdir(parents=True, exist_ok=True)
    return path
