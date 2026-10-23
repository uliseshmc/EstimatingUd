"""Estimating Ud, the deleterious mutation rate, from Human/Chimpanzee/Orangutan alignments.

Importing this package registers the ``GNC_CpG_ss`` and ``GT_CpG_ss`` substitution
models with cogent3 (see :mod:`estimating_ud.models`), so
``get_app("model", "GT_CpG_ss", ...)`` works after ``import estimating_ud``.
"""

from estimating_ud import models as _models  # noqa: F401  (registration side effect)

__version__ = "0.1.0"
