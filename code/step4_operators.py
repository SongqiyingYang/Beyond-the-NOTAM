"""Step 4: operator-associated relative standoff and cross-event consistency (Section III-C).

Detour arcs (per closure-exposed OD pair)
  * the closure set is replaced by a smoothed outline (morphological closing with radius
    OUTLINE_SMOOTH_KM, largest component), whose boundary is a ring with arc length s;
  * every baseline flight of the OD pair that entered the closure (clearance 0) gives an entry
    point A and an exit point B: the first and last intersections of its trajectory, joined
    across data gaps, with the ring. The flight whose (A, B) is closest to all others (sum of
    ring distances) is the reference of the OD pair;
  * A and B split the ring into two arcs, sampled every ARC_STEP_KM. An event-period flight is
    assigned to the arc with the smaller mean distance to its trajectory (detour direction);
  * d_f(s) is the distance from arc position s to the flight's trajectory (densified to
    TRACK_DENSIFY_KM; nearest-point queries with a k-d tree);
  * flights of one OD pair and detour direction share a common arc segment: from the first to
    the last s at which the lower envelope of their d_f(s) is at most ARC_END_KM.
Trajectory segments separated by gaps longer than GAP_S with fewer than MIN_SEGMENT_POINTS
points are dropped as position outliers.

Standoff and operator scores
  * m_f = mean of the smallest share (ARC_FRACTIONS; main 0.20) of d_f(s) on the common segment;
  * R_k (Eq. 15) = mean over operator k's comparable flights of m_f minus the mean m of the other
    operators' flights in the same OD-pair x detour-direction group; a flight is comparable when
    its group contains at least one other operator.

Cross-event consistency
  * Spearman rho of R_k between 2025 and 2026, operators with at least OPERATOR_MIN_FLIGHTS
    comparable flights in each episode (main 10; sensitivity 15, 20);
  * permutation test: operator labels are permuted within each OD-pair x detour-direction
    group (preserving every operator's number of flights per group), B_PERMUTATION times with
    one generator seeded SEED_PERMUTATION; the same permutations serve all arc fractions and
    thresholds. Two-sided p = P(|rho*| >= |rho|);
    family p = max-statistic adjustment across the four arc fractions (max |rho*| over fractions);
  * sensitivity of the reference arcs: references restricted to baseline crossings without a data
    gap longer than REFERENCE_MAX_GAP_S between A and B (OD pairs without such a crossing are
    dropped); permutations from the same generator, after those of the main references.

Usage:
  python step4_operators.py --flights-2025 work/flights_2025.parquet --flights-2026 work/flights_2026.parquet \
      --tracks-2025 TRACKS_2025/ --tracks-2026-baseline TRACKS_2026/baseline/ \
      --tracks-2026-event TRACKS_2026/event/ --firs firs.geojson --out results/main
Outputs in --out:
  tables/S4_1_operator_relative_standoff.csv   R_k and comparable flights of the operators in the main panel
  tables/S4_2_rank_correlations.csv            rho and p for every reference x threshold x arc fraction
  tables/S4_3_reference_arcs.csv               OD pairs with a reference, flights with a profile
  private/operator_flights_*.parquet           flight-level m_f (do not share; see README)
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import shapely
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
from shapely.geometry import LineString

import common as C

EPISODES = ("2025", "2026")


# ----------------------------------------------------------------------------- geometry
def smoothed_outline(geom_xy):
    r = C.OUTLINE_SMOOTH_KM
    o = geom_xy.buffer(r, quad_segs=16).buffer(-r, quad_segs=16).simplify(C.OUTLINE_SIMPLIFY_KM)
    if o.geom_type == "MultiPolygon":
        o = max(o.geoms, key=lambda p: p.area)
    return o.exterior


def ring_distance(a, b, L):
    d = np.abs(a - b) % L
    return np.minimum(d, L - d)


# ----------------------------------------------------------------------------- trajectories
def load_tracks(dirs, ids):
    """Trajectory points of the given flights; segments split at gaps > GAP_S, short isolated
    segments removed; projected to km."""
    ids = set(ids)
    parts = []
    for d in dirs:
        for fp in sorted(Path(d).glob("*.parquet")):
            df = pd.read_parquet(fp, columns=["fr24_id", "timestamp", "lat", "lon"])
            parts.append(df[df.fr24_id.isin(ids)])
    df = pd.concat(parts, ignore_index=True)
    df["timestamp"] = pd.to_datetime(df.timestamp, utc=True)
    df = (df.drop_duplicates(["fr24_id", "timestamp"])
          .sort_values(["fr24_id", "timestamp"]).reset_index(drop=True))
    dt = df.groupby("fr24_id").timestamp.diff().dt.total_seconds()
    df["seg"] = (dt.isna() | (dt > C.GAP_S)).cumsum()
    df = df[df.groupby("seg").seg.transform("size") >= C.MIN_SEGMENT_POINTS].reset_index(drop=True)
    xy = C.project(df.lon.values, df.lat.values)
    df["x"], df["y"] = xy[:, 0], xy[:, 1]
    return df


def entry_exit(g, ring, max_gap=None):
    """Ring positions of the first and last intersection of the trajectory (joined across gaps)
    with the outline; None if fewer than two, or (max_gap) if the crossing contains a longer gap."""
    line = LineString(g[["x", "y"]].values)
    hit = line.intersection(ring)
    pts = [q for q in getattr(hit, "geoms", [hit]) if q.geom_type == "Point" and not q.is_empty]
    if len(pts) < 2:
        return None
    pos = np.array([line.project(q) for q in pts])
    if max_gap is not None:
        xy = g[["x", "y"]].values
        cum = np.r_[0, np.cumsum(np.hypot(*np.diff(xy, axis=0).T))]
        dt = np.diff(g.timestamp.values).astype("timedelta64[s]").astype(float)
        inside = (cum[1:] > pos.min()) & (cum[:-1] < pos.max())
        if (dt[inside] > max_gap).any():
            return None
    return ring.project(pts[int(pos.argmin())]), ring.project(pts[int(pos.argmax())])


def track_tree(g):
    out = [shapely.get_coordinates(shapely.segmentize(LineString(s[["x", "y"]].values), C.TRACK_DENSIFY_KM))
           for _, s in g.groupby("seg") if len(s) >= 2]
    return cKDTree(np.vstack(out)) if out else None


# ----------------------------------------------------------------------------- flight-level standoff
def flight_standoff(F, base_tracks, event_tracks, ring, max_gap=None):
    """m_f for every closure-exposed event flight with a reference arc. Returns (flights, refs)."""
    L = ring.length
    X = C.split_exposure(F)["exposed"]
    base_od = X[X.period == "baseline"].set_index("fr24_id").od
    ev = X[X.period == "event"].set_index("fr24_id")
    ab = {}
    for fid, g in base_tracks.groupby("fr24_id"):
        if fid in base_od.index:
            r = entry_exit(g, ring, max_gap)
            if r is not None:
                ab.setdefault(base_od[fid], []).append((fid, *r))
    te = dict(tuple(event_tracks.groupby("fr24_id")))
    rows = []
    n_near = None
    for od, lst in ab.items():
        R = pd.DataFrame(lst, columns=["fr24_id", "sA", "sB"])
        cost = [(ring_distance(r.sA, R.sA.values, L) + ring_distance(r.sB, R.sB.values, L)).sum()
                for r in R.itertuples()]
        m = R.iloc[int(np.argmin(cost))]
        arcs = {}
        for sign in (+1, -1):
            s = np.arange(0, ((m.sB - m.sA) * sign) % L + 1e-9, C.ARC_STEP_KM)
            arcs[sign] = shapely.get_coordinates(shapely.line_interpolate_point(ring, (m.sA + sign * s) % L))
        north = max(arcs, key=lambda k: arcs[k][:, 1].mean())
        name = {north: "N", -north: "S"}
        per = []
        for fid in ev.index[ev.od == od]:
            g = te.get(fid)
            tree = track_tree(g) if g is not None else None
            if tree is None:
                continue
            dist = {k: tree.query(arcs[k])[0] for k in arcs}
            side = min(dist, key=lambda k: dist[k].mean())
            per.append((fid, side, dist[side]))
        for side in arcs:
            ds = [d for _, sd, d in per if sd == side]
            if not ds:
                continue
            near = np.flatnonzero(np.min(ds, axis=0) <= C.ARC_END_KM)
            if not len(near):
                continue
            a, b = near[0], near[-1]
            for fid, sd, d in per:
                if sd != side:
                    continue
                v = np.sort(d[a:b + 1])
                r = dict(fr24_id=fid, od=od, operator=ev.at[fid, "operator"], direction=name[side],
                         segment_km=(b - a) * C.ARC_STEP_KM)
                for q in C.ARC_FRACTIONS:
                    r[f"m{int(q * 100)}_km"] = float(v[:max(1, int(round(q * len(v))))].mean())
                rows.append(r)
    refs = dict(exposed_od_pairs=int(X.od.nunique()), od_pairs_with_reference=len(ab),
                exposed_event_flights=int(len(ev)))
    return pd.DataFrame(rows), refs


# ----------------------------------------------------------------------------- operator scores
def relative_standoff(Fl, col):
    """R_k of Eq. (15) and the number of comparable flights per operator."""
    D = Fl.dropna(subset=["operator"]).assign(group=lambda x: x.od + "|" + x.direction)
    s = D.groupby(["group", "operator"])[col].agg(["sum", "size"])
    tot = s.groupby(level=0).sum()
    g = s.index.get_level_values(0)
    n_oth = tot["size"].reindex(g).values - s["size"].values
    ok = n_oth > 0
    mean_oth = (tot["sum"].reindex(g).values[ok] - s["sum"].values[ok]) / n_oth[ok]
    r = pd.DataFrame({"operator": s.index.get_level_values(1)[ok],
                      "dev": s["sum"].values[ok] - s["size"].values[ok] * mean_oth, "n": s["size"].values[ok]})
    t = r.groupby("operator")
    return pd.DataFrame({"R_km": t.dev.sum() / t.n.sum(), "comparable_flights": t.n.sum()})


def correlations(FL):
    """{(fraction, min_n): (rho, n_operators)} for one pair of flight tables {episode: Fl}."""
    out = {}
    for q in C.ARC_FRACTIONS:
        col = f"m{int(q * 100)}_km"
        j = relative_standoff(FL["2025"], col).join(relative_standoff(FL["2026"], col),
                                                    lsuffix="_2025", rsuffix="_2026", how="inner")
        for t in C.OPERATOR_MIN_FLIGHTS:
            k = j[(j.comparable_flights_2025 >= t) & (j.comparable_flights_2026 >= t)]
            out[(q, t)] = (spearmanr(k.R_km_2025, k.R_km_2026)[0] if len(k) >= 5 else np.nan, len(k))
    return out


def permutation_tests(FL, rng, B=C.B_PERMUTATION):
    obs = correlations(FL)
    idx = {e: list(pd.Series(np.arange(len(FL[e]))).groupby((FL[e].od + "|" + FL[e].direction).values)
                   .indices.values()) for e in EPISODES}
    null = {k: np.empty(B) for k in obs}
    for b in range(B):
        P = {}
        for e in EPISODES:
            a = FL[e].operator.values.copy()
            for ii in idx[e]:
                a[ii] = rng.permutation(a[ii])
            P[e] = FL[e].assign(operator=a)
        for k, (r, _) in correlations(P).items():
            null[k][b] = r
    rows = []
    for t in C.OPERATOR_MIN_FLIGHTS:
        fam = np.nanmax(np.abs(np.column_stack([null[(q, t)] for q in C.ARC_FRACTIONS])), axis=1)
        for q in C.ARC_FRACTIONS:
            r, n = obs[(q, t)]
            z = null[(q, t)]
            rows.append(dict(min_comparable_flights=t, arc_fraction=q, operators=n, rho_s=r,
                             p_two_sided=float(np.mean(np.abs(z) >= abs(r))),
                             p_family_two_sided=float(np.mean(fam >= abs(r)))))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- driver
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for ep in EPISODES:
        ap.add_argument(f"--flights-{ep}", required=True, help="output of step1_clearance.py")
        ap.add_argument(f"--tracks-{ep}", help="trajectory directory of both periods")
        ap.add_argument(f"--tracks-{ep}-baseline")
        ap.add_argument(f"--tracks-{ep}-event")
    ap.add_argument("--firs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--unit", choices=("nm", "km"), default="nm")
    a = vars(ap.parse_args())
    U = C.Unit(a["unit"])
    out = Path(a["out"]); (out / "tables").mkdir(parents=True, exist_ok=True); (out / "private").mkdir(exist_ok=True)
    firs = C.load_firs(a["firs"])

    data = {}
    for ep in EPISODES:
        _, ref_set = C.reference_closure(ep)
        F = C.analysis_sample(pd.read_parquet(a[f"flights_{ep}"]), ref_set)
        ring = smoothed_outline(C.closure_geometry(firs, ref_set))
        both = a[f"tracks_{ep}"]
        dirs = {p: [a[f"tracks_{ep}_{p}"] or both] for p in ("baseline", "event")}
        X = C.split_exposure(F)["exposed"]
        base_ids = X[(X.period == "baseline") & (X.clearance_km <= 0)].fr24_id
        tb = load_tracks(dirs["baseline"], base_ids)
        te = load_tracks(dirs["event"], X[X.period == "event"].fr24_id)
        data[ep] = (F, tb, te, ring)
        print(f"[step4/{ep}] closure-exposed OD pairs {X.od.nunique():,}; baseline crossings "
              f"{len(base_ids):,}; event flights {(X.period == 'event').sum():,}")

    rng = np.random.default_rng(C.SEED_PERMUTATION)
    tests, refrows, FLmain = [], [], None
    for ref_name, max_gap in (("all baseline crossings", None),
                              (f"crossings without gaps > {C.REFERENCE_MAX_GAP_S // 60} min", C.REFERENCE_MAX_GAP_S)):
        FL = {}
        for ep in EPISODES:
            F, tb, te, ring = data[ep]
            FL[ep], refs = flight_standoff(F, tb, te, ring, max_gap)
            n_cmp = relative_standoff(FL[ep], "m20_km").comparable_flights.sum()
            refrows.append(dict(reference=ref_name, episode=ep, **refs,
                                event_flights_with_profile=len(FL[ep]), comparable_flights=int(n_cmp)))
            print(f"[step4/{ep}] {ref_name}: {refs['od_pairs_with_reference']:,} OD pairs with a reference, "
                  f"{len(FL[ep]):,} event flights with a profile")
        tests.append(permutation_tests(FL, rng).assign(reference=ref_name))
        if max_gap is None:
            FLmain = FL

    # S4.1: main specification (all crossings, 20% of the arc)
    col = f"m{int(C.ARC_FRACTION_MAIN * 100)}_km"
    S = {ep: relative_standoff(FLmain[ep], col) for ep in EPISODES}
    T = S["2025"].join(S["2026"], lsuffix="_2025", rsuffix="_2026", how="outer")
    T["main_panel"] = (T.comparable_flights_2025 >= C.OPERATOR_MIN_MAIN) & (T.comparable_flights_2026 >= C.OPERATOR_MIN_MAIN)
    m = T[T.main_panel]
    u25, u26 = (m.R_km_2025.rank() - 0.5) / len(m), (m.R_km_2026.rank() - 0.5) / len(m)
    T["tercile_both"] = ""
    T.loc[m.index[(u25 > 2 / 3) & (u26 > 2 / 3)], "tercile_both"] = "upper"
    T.loc[m.index[(u25 < 1 / 3) & (u26 < 1 / 3)], "tercile_both"] = "lower"
    for c in ("comparable_flights_2025", "comparable_flights_2026"):
        T[c] = T[c].fillna(0).astype(int)
    T = T[T.main_panel].drop(columns="main_panel")
    T = T.reset_index().rename(columns={"index": "operator"}).sort_values("operator")
    for c in ("R_km_2025", "R_km_2026"):
        T[c] = U(T[c], 1)
    T = T.rename(columns={c: c.replace("_km", f"_{U.tag}") for c in T.columns})
    T.to_csv(out / "tables" / "S4_1_operator_relative_standoff.csv", index=False)

    R = pd.concat(tests)[["reference", "min_comparable_flights", "arc_fraction", "operators", "rho_s",
                          "p_two_sided", "p_family_two_sided"]]
    R.round(4).to_csv(out / "tables" / "S4_2_rank_correlations.csv", index=False)
    pd.DataFrame(refrows).to_csv(out / "tables" / "S4_3_reference_arcs.csv", index=False)
    for ep in EPISODES:
        FLmain[ep].to_parquet(out / "private" / f"operator_flights_{ep}.parquet", index=False)
    main_row = R[(R.reference == "all baseline crossings") & (R.min_comparable_flights == C.OPERATOR_MIN_MAIN)
                 & (R.arc_fraction == C.ARC_FRACTION_MAIN)].iloc[0]
    print(f"[step4] main: {int(main_row.operators)} operators, rho_S = {main_row.rho_s:.3f}, "
          f"two-sided p = {main_row.p_two_sided:.3f}, family p = {main_row.p_family_two_sided:.3f}")
    print(f"[step4] done -> {out}")


if __name__ == "__main__":
    main()
