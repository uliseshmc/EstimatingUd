"""One module per pipeline stage. Each exposes ``run(...)`` with plain arguments; the CLI
in :mod:`estimating_ud.cli` only parses options and calls these functions."""
