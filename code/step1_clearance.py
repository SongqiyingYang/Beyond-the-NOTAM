"""Step 1: minimum horizontal clearance of every flight from the closure set (Eq. 2).

Input
  --tracks     directory of parquet files with one row per trajectory point on the 60-s grid:
               fr24_id, timestamp (UTC), lat, lon
  --flights    parquet or CSV with one row per flight:
               fr24_id, operator (ICAO airline designator), orig_icao, dest_icao,
               category ("Passenger" for passenger operations), flight (flight number; used
               only for the service-level exposure definition)
  --firs       FIR polygons with an IDENT column (see code/README.md)

Rules (Section II)
  * points outside the region of interest (25-75 E, 15-45 N) are dropped;
  * a day runs from 06:00 to 06:00 UTC; a flight belongs to the day of its first point in the
    region; points of days outside the baseline or evaluation window are dropped;
  * baseline flights are measured against the closure set of the reference evaluation day
    (the day with the most closed FIRs); event flights against the closure set of their day;
  * distances use an equirectangular projection centred at 30 N;
  * for each flight the output keeps the minimum clearance, the closed FIR nearest to that
    point, and the position at closest approach; a flight seen in several files keeps its
    minimum per period.

Variants (Section III-D)
  --trim-km R      ignore trajectory points within R km of the flight's own origin and
                   destination airports (airport coordinates from the airportsdata package,
                   optionally filled from --airports)
  --drop-fir FIR   remove a FIR from the closure set used for distances (e.g. OMAE in 2026,
                   OSTT in 2025); the analysis-sample screening is set in run_all.py

Output: parquet, one row per flight and period, with columns
  fr24_id, period, day, clearance_km, nearest_fir, near_lat, near_lon, n_pts,
  operator, orig_icao, dest_icao, category, flight, od

Usage: python step1_clearance.py --episode 2026 --tracks DIR --flights FILE --firs FILE
                                 --out flights_2026.parquet [--trim-km 50] [--drop-fir OMAE]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import shapely

import common as C


def airport_xy(codes, extra=None):
    """Projected airport positions (see common.airport_coordinates)."""
    ap = C.airport_coordinates(extra)
    ll = np.array([ap.get(str(c), (np.nan, np.nan)) for c in codes], float)
    return pd.DataFrame(C.project(ll[:, 1], ll[:, 0]), index=list(codes), columns=["x", "y"])


def closure_by_day(episode, firs, drop):
    ref, ref_set = C.reference_closure(episode)
    ref_set = [f for f in ref_set if f not in drop]

    def geoms(cs):
        cs = sorted(f for f in cs if f not in drop)
        return C.closure_geometry(firs, cs), {f: C.project_geom(firs.loc[f]) for f in cs}

    days = {d: geoms(C.closed_set(episode, d)) for d in C.EPISODES[episode]["eval_days"]}
    return ref, ref_set, geoms(ref_set), days


def flights_in_file(fp, episode, period, ref_geoms, day_geoms, names, meta_xy, trim_km):
    days = C.EPISODES[episode]["baseline_days" if period == "baseline" else "eval_days"]
    df = pd.read_parquet(fp, columns=["fr24_id", "timestamp", "lat", "lon"])
    df["timestamp"] = pd.to_datetime(df.timestamp, utc=True)
    df = df[df.lon.between(C.ROI["lon_min"], C.ROI["lon_max"])
            & df.lat.between(C.ROI["lat_min"], C.ROI["lat_max"])]
    if df.empty:
        return None
    df["day"] = (df.timestamp - pd.Timedelta(hours=C.DAY_SLICE_HOUR)).dt.date.astype(str)
    df = df[df.day.isin(days)].sort_values(["fr24_id", "timestamp"])
    if df.empty:
        return None
    xy = C.project(df.lon.values, df.lat.values)
    if trim_km > 0 and meta_xy is not None:
        keep = np.ones(len(df), bool)
        for end in ("orig", "dest"):
            a = meta_xy[end].reindex(df.fr24_id.values).to_numpy()
            d = np.hypot(xy[:, 0] - a[:, 0], xy[:, 1] - a[:, 1])
            keep &= ~(d < trim_km)                        # unknown airport (nan) -> not trimmed
        df, xy = df[keep], xy[keep]
        if df.empty:
            return None
    pts = shapely.points(xy[:, 0], xy[:, 1])
    first_day = df.groupby("fr24_id").day.first()
    dist = np.empty(len(df))
    per_fir = np.full((len(df), len(names)), np.inf)
    for d, idx in df.groupby("day").indices.items():
        g, fg = day_geoms.get(d, ref_geoms) if period == "event" else ref_geoms
        dist[idx] = shapely.distance(pts[idx], g)
        for j, nm in enumerate(names):
            if nm in fg:
                per_fir[idx, j] = shapely.distance(pts[idx], fg[nm])
    df = df.assign(dist=dist, near=np.array(names)[per_fir.argmin(1)])
    agg = df.groupby("fr24_id").agg(clearance_km=("dist", "min"), n_pts=("dist", "size"))
    at = df.loc[df.groupby("fr24_id").dist.idxmin()].set_index("fr24_id")
    agg["nearest_fir"], agg["near_lat"], agg["near_lon"] = at.near, at.lat, at.lon
    agg["day"], agg["period"] = first_day, period
    return agg.reset_index()


def build(episode, tracks, flights, firs, trim_km=0.0, drop=(), airports=None, periods=None):
    firs = C.load_firs(firs)
    meta = (pd.read_parquet(flights) if str(flights).endswith(".parquet") else pd.read_csv(flights))
    meta = meta.drop_duplicates("fr24_id").set_index("fr24_id")
    ref, ref_set, ref_geoms, day_geoms = closure_by_day(episode, firs, set(drop))
    names = sorted(ref_set)
    print(f"[step1/{episode}] reference day {ref}, closure set {names}"
          + (f", trimming {trim_km:g} km around airports" if trim_km else ""))
    meta_xy = None
    if trim_km > 0:
        meta_xy = {end: airport_xy(meta[f"{end}_icao"].astype(str).unique(), airports)
                   .reindex(meta[f"{end}_icao"].astype(str).values).set_axis(meta.index)
                   for end in ("orig", "dest")}
    periods = periods or {"baseline": Path(tracks), "event": Path(tracks)}
    rows = []
    for period, d in periods.items():
        for fp in sorted(Path(d).glob("*.parquet")):
            r = flights_in_file(fp, episode, period, ref_geoms, day_geoms, names, meta_xy, trim_km)
            if r is not None:
                rows.append(r)
    F = (pd.concat(rows, ignore_index=True).sort_values("clearance_km")
         .drop_duplicates(["fr24_id", "period"]).reset_index(drop=True))
    F = F.join(meta[[c for c in ("operator", "orig_icao", "dest_icao", "category", "flight")
                     if c in meta.columns]], on="fr24_id")
    F["od"] = F.orig_icao.astype(str) + ">" + F.dest_icao.astype(str)
    print(f"[step1/{episode}] {len(F):,} flight-periods")
    return F


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--episode", required=True, choices=sorted(C.EPISODES))
    ap.add_argument("--tracks", help="directory with the trajectory files of both periods")
    ap.add_argument("--tracks-baseline", help="separate directory for baseline trajectories")
    ap.add_argument("--tracks-event", help="separate directory for event-period trajectories")
    ap.add_argument("--flights", required=True)
    ap.add_argument("--firs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--trim-km", type=float, default=0.0)
    ap.add_argument("--drop-fir", action="append", default=[])
    ap.add_argument("--airports", help="optional CSV filling airportsdata: OurAirports airports.csv "
                    "or a table with columns icao, lat, lon")
    a = ap.parse_args()
    periods = None
    if a.tracks_baseline or a.tracks_event:
        periods = {"baseline": Path(a.tracks_baseline), "event": Path(a.tracks_event)}
    F = build(a.episode, a.tracks, a.flights, a.firs, a.trim_km, a.drop_fir, a.airports, periods)
    F.to_parquet(a.out, index=False)


if __name__ == "__main__":
    main()
