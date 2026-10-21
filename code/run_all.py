"""Run steps 2 and 3 on the flight tables of both episodes and write the result tables S2-S3.

Step 4 (operator-relative standoff) needs the trajectories again and is run separately with
step4_operators.py (see code/README.md).

Usage (main specification):
  python run_all.py --flights-2025 flights_2025.parquet --flights-2026 flights_2026.parquet \
                    --firs firs.geojson --out results/main

Sensitivity variants (Section III-D; see code/README.md):
  airport trimming   build the flight tables with step1_clearance.py --trim-km 50 (or 100)
  OMAE codings       build the 2026 flight table with step1_clearance.py --drop-fir OMAE, then
                     pass --geometry-drop 2026:OMAE (closure boundary without OMAE) and, for the
                     "open" coding, also --endpoint-drop 2026:OMAE (flights to and from UAE
                     airports admitted)

Outputs in --out:
  tables/   S2_*, S3_* in the format of the repository's tables/ (distances in --unit)
  summary.json   headline quantities used by collect_variants.py
  private/  flight-level residuals (input of Fig. 7); do not share if your data licence
            forbids redistribution of flight-level data
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import common as C
import step2_profile as S2
import step3_floor as S3


def _drops(tokens, episode):
    """'FIR' applies to both episodes, 'EPISODE:FIR' (e.g. 2026:OMAE) to one episode only."""
    return {t.split(":")[-1] for t in tokens if ":" not in t or t.split(":")[0] == episode}


def load_sample(path, episode, geometry_drop, endpoint_drop):
    F = pd.read_parquet(path)
    _, ref_set = C.reference_closure(episode)
    excl = [f for f in ref_set if f not in _drops(endpoint_drop, episode)]
    geom = [f for f in ref_set if f not in _drops(geometry_drop, episode)]
    return C.analysis_sample(F, excl), geom


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flights-2025", required=True)
    ap.add_argument("--flights-2026", required=True)
    ap.add_argument("--firs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--unit", choices=("nm", "km"), default="nm")
    ap.add_argument("--geometry-drop", action="append", default=[])
    ap.add_argument("--endpoint-drop", action="append", default=[])
    a = ap.parse_args()
    U = C.Unit(a.unit)
    out = Path(a.out); (out / "tables").mkdir(parents=True, exist_ok=True)
    (out / "private").mkdir(exist_ok=True)
    firs = C.load_firs(a.firs)

    def save(df, name):
        df.to_csv(out / "tables" / name, index=False)
        print(f"  -> tables/{name}")

    samples, geoms, summary = {}, {}, {}
    prof, fits, cells_rows, fir_rows, dec, qrows, srows, quant = [], [], [], [], [], [], [], []
    resid_private = []
    for ep, path in (("2026", a.flights_2026), ("2025", a.flights_2025)):
        F, geom = load_sample(path, ep, a.geometry_drop, a.endpoint_drop)
        gxy = C.closure_geometry(firs, geom)
        samples[ep], geoms[ep] = F, gxy
        X = C.split_exposure(F)
        print(f"[run/{ep}] analysis sample {len(F):,} flight-periods; closure-exposed OD pairs "
              f"{X['exposed'].od.nunique():,}, event flights {(X['exposed'].period == 'event').sum():,}")

        # ---- step 2 (one random stream per episode: unexposed OD pairs first, then exposed)
        rng_profile = np.random.default_rng(C.SEED_PROFILE)
        P = {}
        for lab in ("unexposed", "exposed"):
            P[lab] = S2.profile_with_ci(X[lab], rng_profile)
        for lab in ("exposed", "unexposed"):
            p, _ = P[lab]
            prof.append(p.assign(episode=ep, od_pairs="closure-exposed" if lab == "exposed" else lab))
        pe = P["exposed"][0]
        bf = S2.benchmark_fits(pe)
        fits.append(dict(episode=ep, **bf, h_at_anchor=float(np.interp(C.D_HALF_ANCHOR, pe.d_km, pe.h_hat)),
                         d_half_km=S2.half_decay(pe.d_km.values, pe.h_hat.values)))
        mono = bool(np.all(np.diff(pe.h_hat.values[(pe.d_km >= 20) & (pe.d_km <= 200)]) <= 1e-9))

        # ---- step 3
        Fa = S3.anchor_cells(F, gxy)
        cells, fl = S3.cell_floors(Fa), S3.fir_floors(Fa)
        E, n_without = S3.residuals(Fa, cells, fl)
        lo, hi = S3.median_ci(E)
        val = S3.validation_cells(Fa, cells)
        ev = S3.exposed_event(Fa)
        g = ev.groupby("cell").clearance_km
        for c, r in cells.iterrows():
            i, j = map(int, c.split("_"))
            v = val.loc[c] if c in val.index else None
            cells_rows.append(dict(episode=ep, cell=c, lat_min=i * C.CELL_DEG, lat_max=(i + 1) * C.CELL_DEG,
                                   lon_min=j * C.CELL_DEG, lon_max=(j + 1) * C.CELL_DEG,
                                   baseline_positive_flights=int(r.n), baseline_flights_near=int(r.n_near),
                                   floor_p5_km=r.floor, eligible=bool(r.eligible),
                                   exposed_event_flights=int(g.size().get(c, 0)),
                                   event_p10_km=v.p10 if v is not None else np.nan,
                                   validation_cell=v is not None))
        for fir, r in fl.iterrows():
            fir_rows.append(dict(episode=ep, fir=fir, baseline_positive_flights=int(r.n), floor_p5_km=r.floor))
        dec.append(dict(episode=ep, exposed_event_flights=len(E) + n_without, flights_without_floor=n_without,
                        share_cell_floor=(E.floor_source == "cell").mean(),
                        median_clearance_km=E.clearance_km.median(), median_floor_km=E.floor.median(),
                        median_residual_km=E.residual.median(), median_residual_lo_km=lo,
                        median_residual_hi_km=hi, share_negative_residual=(E.residual < 0).mean()))
        qs = S3.quantile_sensitivity(Fa, E, fl)
        qrows += [dict(episode=ep, floor_quantile_q=q, main=q == 0.05, median_residual_km=v) for q, v in qs.items()]
        srows.append(S3.eligibility_sensitivity(Fa, gxy, fl).assign(episode=ep))
        quant.append(dict(episode=ep, **{f"p{int(q * 100):02d}_km": E.residual.quantile(q)
                                         for q in (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)}))
        resid_private.append(E[["residual"]].assign(episode=ep))
        val.reset_index().rename(columns={"index": "cell"}).assign(episode=ep).to_csv(
            out / "private" / f"validation_cells_{ep}.csv", index=False)

        summary[ep] = dict(exposed_od_pairs=int(X["exposed"].od.nunique()),
                           exposed_event_flights=int((X["exposed"].period == "event").sum()),
                           profile_monotone=mono,
                           median_clearance_km=float(E.clearance_km.median()), median_floor_km=float(E.floor.median()),
                           median_residual_km=float(E.residual.median()), median_residual_lo_km=lo,
                           median_residual_hi_km=hi, median_residual_q0_km=qs[0], median_residual_q10_km=qs[0.10])

    # ---- exposure definitions (d_1/2 with intervals), both episodes, one random stream
    D = S2.exposure_definitions({ep: C.split_exposure(samples[ep])["both"] for ep in samples},
                                np.random.default_rng(C.SEED_DHALF))
    for ep in ("2025", "2026"):
        m = D[(D.episode == ep) & D.definition.str.contains("main")].iloc[0]
        d = D[D.episode == ep]
        summary[ep].update(d_half_km=m.d_half_km, d_half_lo_km=m.d_half_lo_km, d_half_hi_km=m.d_half_hi_km,
                           d_half_range_7_definitions_km=[d.d_half_km.min(), d.d_half_km.max()],
                           monotone_all_definitions=bool(d.monotone.all()))
        for f in fits:
            if f["episode"] == ep:
                f.update(d_half_lo_km=m.d_half_lo_km, d_half_hi_km=m.d_half_hi_km)

    # ---- write tables (distances converted to the reporting unit)
    def conv(df, cols, nd=1):
        """Convert distance columns in place (same position, unit suffix renamed)."""
        df = df.copy()
        for c in cols:
            df[c] = U(df[c], nd)
        return df.rename(columns={c: c.replace("_km", f"_{U.tag}") for c in cols})

    Pf = pd.concat(prof)
    save(conv(Pf[["episode", "od_pairs", "d_km", "h_hat", "ci_lo", "ci_hi", "se", "event_flights_within_d",
                   "baseline_flights_within_d"]], ["d_km"]).round(4), "S2_1_avoidance_profile.csv")
    Fi = pd.DataFrame(fits)
    lo_, hi_ = U.thr(C.FIT_NODES.min()), U.thr(C.FIT_NODES.max())
    Fi.insert(2, "fit_range", f"{lo_.split()[0]}-{hi_}")
    save(conv(Fi, ["sigma_km", "xi_bar_lo_km", "xi_bar_hi_km", "d_half_km", "d_half_lo_km", "d_half_hi_km"]).round(3),
         "S2_2_benchmark_fits.csv")
    save(conv(D, ["d_half_km", "d_half_lo_km", "d_half_hi_km"]).round(3), "S2_3_exposure_definitions.csv")
    save(conv(pd.DataFrame(cells_rows), ["floor_p5_km", "event_p10_km"]).assign(
        near_boundary_threshold=U.thr(C.CELL_NEAR_KM)), "S3_1_boundary_cells.csv")
    save(conv(pd.DataFrame(fir_rows), ["floor_p5_km"]), "S3_2_fir_fallback_floors.csv")
    save(conv(pd.DataFrame(dec), ["median_clearance_km", "median_floor_km", "median_residual_km",
                                  "median_residual_lo_km", "median_residual_hi_km"]).round(3), "S3_3_decomposition.csv")
    save(conv(pd.DataFrame(qrows), ["median_residual_km"]), "S3_4_floor_quantile.csv")
    S5 = pd.concat(srows)
    S5["near_boundary_threshold"] = S5.pop("near_boundary_km").map(U.thr)
    save(conv(S5, ["median_residual_km"]).round(3), "S3_5_eligibility_cell_size.csv")
    Q = pd.DataFrame(quant)
    save(conv(Q, [c for c in Q.columns if c.endswith("_km")]), "S3_6_residual_quantiles.csv")
    pd.concat(resid_private).to_parquet(out / "private" / "residuals.parquet")
    with open(out / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=1, default=float)
    print(f"[run] done -> {out}")


if __name__ == "__main__":
    main()
