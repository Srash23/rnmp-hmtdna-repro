#!/usr/bin/env python
"""Execute a notebook in-process and write the outputs back into the .ipynb.

    python scripts/execute_notebook.py notebooks/rnmp_hmtdna_reproduction.ipynb

WHY NOT nbclient/jupyter: a Jupyter kernel needs to bind a local socket, which is
blocked in the sandbox this analysis was developed in (PermissionError from
jupyter_client's port finder). Rather than ship a notebook with empty outputs,
this runs each cell in the current interpreter, captures stdout, any figures the
cell created, and the value of a trailing expression, and injects them as cell
outputs. Use `jupyter nbconvert --execute` instead wherever that works -- this is
a fallback, not a better tool.

Cells are executed in order in one shared namespace, exactly as a kernel would.
Any exception propagates, so a broken notebook fails loudly.
"""
import ast
import base64
import contextlib
import io
import json
import os
import sys
import uuid
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def run(path: Path) -> dict:
    nb = json.loads(path.read_text())
    ns, n = {}, 0
    for cell in nb["cells"]:
        cell.setdefault("id", uuid.uuid4().hex[:8])
        if cell["cell_type"] != "code":
            continue
        src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
        n += 1
        cell["execution_count"] = n
        outs, buf = [], io.StringIO()
        before = set(plt.get_fignums())
        tree = ast.parse(src)
        tail = None
        if tree.body and isinstance(tree.body[-1], ast.Expr):
            tail = ast.Expression(tree.body.pop().value)
        with contextlib.redirect_stdout(buf):
            exec(compile(tree, "<nb>", "exec"), ns, ns)
            val = eval(compile(tail, "<nb>", "eval"), ns, ns) if tail else None
        if buf.getvalue():
            outs.append({"output_type": "stream", "name": "stdout", "text": buf.getvalue()})
        for num in [f for f in plt.get_fignums() if f not in before]:
            fig = plt.figure(num)
            b = io.BytesIO()
            fig.savefig(b, format="png", dpi=110, bbox_inches="tight")
            outs.append({"output_type": "display_data", "metadata": {},
                         "data": {"image/png": base64.b64encode(b.getvalue()).decode()}})
            plt.close(fig)
        if val is not None:
            data = {"text/plain": repr(val)}
            if hasattr(val, "to_html"):
                data["text/html"] = val.to_html()
            outs.append({"output_type": "execute_result", "execution_count": n,
                         "data": data, "metadata": {}})
        cell["outputs"] = outs
    path.write_text(json.dumps(nb, indent=1))
    return nb


if __name__ == "__main__":
    p = Path(sys.argv[1]).resolve()
    # Run with the notebook's own directory as cwd, so the command works from
    # anywhere -- `python scripts/execute_notebook.py notebooks/x.ipynb` from the
    # repo root and `python ../scripts/execute_notebook.py x.ipynb` from
    # notebooks/ must behave identically.
    os.chdir(p.parent)
    nb = run(p)
    code = [c for c in nb["cells"] if c["cell_type"] == "code"]
    imgs = sum(1 for c in code for o in c["outputs"] if "image/png" in o.get("data", {}))
    tabs = sum(1 for c in code for o in c["outputs"] if "text/html" in o.get("data", {}))
    print(f"{p}: {len(code)} code cells executed, {imgs} inline figures, {tabs} rendered tables")
