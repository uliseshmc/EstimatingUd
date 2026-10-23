"""Where the data and the results live.

The alignment download root is read from the ``DATA_HUMCHIMPORANGOR114`` environment
variable (or the ``--data-root`` CLI option). When the variable is not exported, a
``.env`` file at the repository root is consulted (``KEY=value`` lines, ``export`` prefix
allowed), so notebooks and shells that did not source ``~/.bashrc`` still find the data.
Results default to ``<repo>/results``, located from this file when the package is
installed in editable mode, and can be overridden with ``ESTIMATING_UD_RESULTS``.
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


def dotenv_path() -> Path | None:
    """``<repo>/.env`` if the repository root can be located."""
    root = repo_root()
    return (root / ".env") if root else None


def read_dotenv(path: Path | None = None) -> dict[str, str]:
    """Parse a ``.env`` file: ``KEY=value`` per line, optional ``export``, ``#`` comments, quotes."""
    path = path or dotenv_path()
    values: dict[str, str] = {}
    if not path or not path.is_file():
        return values
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        if line.startswith("export "):
            line = line[len("export ") :].lstrip()
        key, _, value = line.partition("=")
        values[key.strip()] = value.strip().strip("'\"")
    return values


def data_root(override: str | os.PathLike | None = None) -> Path:
    """Return the alignment download root as a Path.

    Order: explicit ``override``, the ``DATA_HUMCHIMPORANGOR114`` environment variable,
    then ``DATA_HUMCHIMPORANGOR114`` in ``<repo>/.env``.
    """
    value = override or os.environ.get(DATA_ENV) or read_dotenv().get(DATA_ENV)
    if not value:
        msg = (
            f"set the {DATA_ENV} environment variable to the alignment download location, "
            f"or put '{DATA_ENV}=/path/to/data' in {dotenv_path() or '<repo>/.env'}"
        )
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
