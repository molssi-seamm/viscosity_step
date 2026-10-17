#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Check the results of the SPC/E regression (spce.flow) against the expected ranges.

Usage::

    python check.py <job directory or job.out> [...]

With several jobs, their average is also checked.
Exits with status 1 if any value is outside its range.
"""

from pathlib import Path
import re
import sys

#: The accepted ranges: about the 95% error of one 3 ns run around the references
VISCOSITY = (0.62, 0.82)  # mPa*s, Green-Kubo average of the five shear components
D0 = (2.75, 3.05)  # 1e-9 m^2/s, oxygen MSD with the Yeh-Hummer correction


def read(path):
    """The viscosity (mPa*s) and corrected diffusion coefficient (1e-9 m^2/s)."""
    path = Path(path)
    if path.is_dir():
        path = path / "job.out"
    text = path.read_text()

    # The last viscosity table: the row "η (Green-Kubo)", or in the tables before
    # phase 8 the first "average" row
    start = text.rfind("Shear viscosity")
    if start < 0:
        raise RuntimeError(f"No viscosity table in {path}")
    match = re.search(r"│\s*η \(Green-Kubo\)\s*│\s*([-0-9.]+)", text[start:])
    if match is None:
        match = re.search(r"│\s*average\s*│\s*([-0-9.]+)", text[start:])
    viscosity = float(match.group(1))

    # The last table of diffusion coefficients: the total, already corrected
    start = text.rfind("Diffusion Coefficients (*")
    if start < 0:
        raise RuntimeError(f"No diffusion coefficients in {path}")
    scale = float(re.search(r"\(\*\s*([0-9.e+-]+)\s*m\^2/s\)", text[start:]).group(1))
    match = re.search(r"│\s*total\s*│\s*([-0-9.]+)", text[start:])
    d0 = float(match.group(1)) * scale / 1.0e-9
    return viscosity, d0


def check(name, value, limits):
    ok = limits[0] <= value <= limits[1]
    result = "ok" if ok else "FAIL"
    print(f"    {name:10s} {value:7.3f}   [{limits[0]}, {limits[1]}]   {result}")
    return ok


def main(paths):
    ok = True
    values = []
    for path in paths:
        viscosity, d0 = read(path)
        values.append((viscosity, d0))
        print(path)
        ok &= check("viscosity", viscosity, VISCOSITY)
        ok &= check("D0", d0, D0)
    if len(values) > 1:
        print(f"Average of {len(values)} runs")
        ok &= check("viscosity", sum(v for v, _ in values) / len(values), VISCOSITY)
        ok &= check("D0", sum(d for _, d in values) / len(values), D0)
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1:]))
