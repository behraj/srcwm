"""Assemble the SRCWM E0 calibration notebook from cells_src/.

Code cells (*.py) and markdown cells (*.md) are concatenated in lexical
filename order into one .ipynb. Run from the repository root:

    python assemble_notebook.py
"""
import sys
from pathlib import Path

import nbformat as nbf

CELLS_DIR = Path(__file__).parent / "cells_src"
OUT = Path(__file__).parent / "srcwm_e0_calibration.ipynb"

nb = nbf.v4.new_notebook()
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python",
                   "name": "python3"},
    "language_info": {"name": "python", "version": sys.version.split()[0]},
}

sources = sorted(CELLS_DIR.iterdir(), key=lambda p: p.name)
for path in sources:
    text = path.read_text().rstrip()
    if path.suffix == ".md":
        nb.cells.append(nbf.v4.new_markdown_cell(text))
    elif path.suffix == ".py":
        nb.cells.append(nbf.v4.new_code_cell(text))
    else:
        raise ValueError(f"unexpected file in cells_src: {path.name}")

nbf.validate(nb)
nbf.write(nb, OUT)
print(f"assembled {len(nb.cells)} cells "
      f"({sum(1 for c in nb.cells if c.cell_type == 'code')} code) -> {OUT.name}")
