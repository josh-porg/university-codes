"""Write tables as LaTeX ``tabular`` environments (``table2latex``, ``latex(sym(...))``)."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Sequence


def _cell(value, fmt: str) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, bool):
        return str(int(value))
    if isinstance(value, float) and math.isinf(value):
        return r"$\infty$" if value > 0 else r"$-\infty$"
    try:
        return format(value, fmt)
    except (TypeError, ValueError):
        return str(value)


def to_latex_table(
    rows,
    headers: Sequence[str] | None = None,
    row_names: Sequence[str] | None = None,
    fmt: str = ".6g",
    path: str | Path | None = None,
) -> str:
    """Render ``rows`` (2-D array or list of lists) as a LaTeX ``tabular``.

    Mirrors the third-party ``table2latex`` used in AE 546 and Math 796:
    left-aligned columns, header line, ``\\hline`` above and below the body.
    Writes the text to ``path`` too if given, and returns it.
    """
    rows = [list(r) for r in rows]
    n_col = len(rows[0]) if rows else len(headers or [])
    spec = "l" * (n_col + (1 if row_names is not None else 0))
    lines = [rf"\begin{{tabular}}{{{spec}}}"]
    if headers is not None:
        head = " & ".join(headers)
        lines.append(("& " if row_names is not None else "") + head + r" \\")
    lines.append(r"\hline")
    for i, r in enumerate(rows):
        cells = [_cell(v, fmt) for v in r]
        if row_names is not None:
            cells.insert(0, str(row_names[i]))
        lines.append(" & ".join(cells) + r" \\")
    lines += [r"\hline", r"\end{tabular}"]
    text = "\n".join(lines)
    if path is not None:
        Path(path).write_text(text)
    return text
