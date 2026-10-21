"""Collect the headline quantities of several run_all.py outputs into one table (Table S5.1).

Usage:
  python collect_variants.py --unit nm --out tables/S5_1_robustness_variants.csv \
      --variant main=results/main \
      --variant trim_50=results/trim_50 \
      --variant omae_boundary=results/omae_boundary@2026 --variant omae_open=results/omae_open@2026

NAME=DIR@EPISODE keeps only that episode (for the OMAE codings, which leave 2025 unchanged).
"""
import argparse
import json
from pathlib import Path

import pandas as pd

import common as C

DESCRIPTION = {
    "main": "main specification",
    "trim_50": "trajectory within {50} of airports excluded",
    "trim_100": "trajectory within {100} of airports excluded",
    "omae_boundary": "2026: OMAE removed from the closure boundary",
    "omae_open": "2026: OMAE open (UAE-endpoint flights admitted)",
}
DIST = ["d_half_km", "d_half_lo_km", "d_half_hi_km", "median_clearance_km", "median_floor_km",
        "median_residual_km", "median_residual_lo_km", "median_residual_hi_km",
        "median_residual_q0_km", "median_residual_q10_km"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", action="append", required=True)
    ap.add_argument("--unit", choices=("nm", "km"), default="nm")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    U = C.Unit(a.unit)
    rows = []
    for spec in a.variant:
        name, rest = spec.split("=", 1)
        path, _, only = rest.partition("@")
        S = json.loads((Path(path) / "summary.json").read_text())
        desc = DESCRIPTION.get(name, name)
        desc = desc.replace("{50}", U.thr(50)).replace("{100}", U.thr(100))
        for ep in ("2025", "2026"):
            if only and ep != only:
                continue
            s = S[ep]
            r = dict(variant=name, description=desc, episode=ep,
                     exposed_od_pairs=s["exposed_od_pairs"], exposed_event_flights=s["exposed_event_flights"],
                     profile_monotone=s["profile_monotone"],
                     monotone_all_7_definitions=s["monotone_all_definitions"])
            for c in DIST:
                r[c.replace("_km", f"_{U.tag}")] = U(s[c])
            lo, hi = s["d_half_range_7_definitions_km"]
            r[f"d_half_range_7_definitions_{U.tag}"] = f"{U(lo)}-{U(hi)}"
            rows.append(r)
    pd.DataFrame(rows).to_csv(a.out, index=False)
    print(f"-> {a.out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
