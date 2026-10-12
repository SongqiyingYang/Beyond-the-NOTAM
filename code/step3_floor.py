"""Step 3: revealed structural floor and residual standoff (Section III-B).

  anchor_cells        boundary anchor a_f of every flight (closest-approach position projected
                      onto the boundary of the closure set) and its 1.5-degree cell
  cell_floors         fifth percentile of positive baseline clearances per cell (Eq. 17) and
                      eligibility: >= 40 such flights, >= 20 of them within 50 km of the boundary
  fir_floors          FIR-level fifth percentile (>= 100 positive baseline clearances), used for
                      flights anchored to an ineligible cell
  residuals           residual standoff of closure-exposed event-period flights (Eq. 18)
  median_ci           OD-clustered bootstrap interval of the median residual (B = 1000)
  quantile_sensitivity, eligibility_sensitivity, validation_cells
"""
import numpy as np
import pandas as pd
import shapely

import common as C


def anchor_cells(F, boundary_geom_xy, deg=C.CELL_DEG):
    F = F.copy()
    ok = F.near_lat.notna().to_numpy()
    xy = C.project(F.near_lon[ok], F.near_lat[ok])
    line = shapely.shortest_line(boundary_geom_xy.boundary, shapely.points(xy[:, 0], xy[:, 1]))
    a = shapely.get_coordinates(line).reshape(-1, 2, 2)[:, 0, :]
    lon, lat = C.unproject(a)
    F.loc[ok, "anchor_lat"], F.loc[ok, "anchor_lon"] = lat, lon
    F["cell"] = cell_id(F.anchor_lat.to_numpy(), F.anchor_lon.to_numpy(), deg)
    return F


def cell_id(lat, lon, deg=C.CELL_DEG):
    return np.array([f"{np.floor(a / deg):.0f}_{np.floor(b / deg):.0f}" if a == a else None
                     for a, b in zip(lat, lon)], dtype=object)


def cell_floors(F, near_km=C.CELL_NEAR_KM, min_near=C.CELL_MIN_NEAR, q=C.FLOOR_Q):
    b = F[(F.period == "baseline") & (F.clearance_km > 0) & F.cell.notna()]
    g = b.groupby("cell").clearance_km
    tab = pd.DataFrame({"n": g.size(), "n_near": b[b.clearance_km < near_km].groupby("cell").size(),
                        "floor": g.min() if q == 0 else g.quantile(q)}).fillna({"n_near": 0})
    tab["eligible"] = (tab.n >= C.CELL_MIN_FLIGHTS) & (tab.n_near >= min_near)
    return tab


def fir_floors(F):
    b = F[(F.period == "baseline") & (F.clearance_km > 0)]
    g = b.groupby("nearest_fir").clearance_km.agg(["size", lambda s: np.percentile(s, 100 * C.FLOOR_Q)])
    g.columns = ["n", "floor"]
    return g[g.n >= C.FIR_MIN_FLIGHTS]


def exposed_event(F):
    E = C.split_exposure(F)["exposed"]
    return E[E.period == "event"].copy()


def residuals(F, cells, firs):
    """Closure-exposed event flights with their floor (cell first, then FIR) and residual."""
    E = exposed_event(F)
    fc = E.cell.map(cells.floor[cells.eligible])
    E["floor_source"] = np.where(fc.notna(), "cell", "fir")
    E["floor"] = fc.fillna(E.nearest_fir.map(firs.floor))
    n_without = int(E.floor.isna().sum())
    E = E[E.floor.notna()].copy()
    E["residual"] = E.clearance_km - E.floor
    return E, n_without


def median_ci(E, seed=C.SEED_RESIDUAL, B=C.B_RESIDUAL):
    rng = np.random.default_rng(seed)
    ods = E.od.values
    uo = np.array(sorted(set(ods)))
    idx = {o: np.where(ods == o)[0] for o in uo}
    v = E.residual.values
    med = [np.median(v[np.concatenate([idx[uo[i]] for i in rng.integers(0, len(uo), len(uo))])])
           for _ in range(B)]
    return tuple(np.percentile(med, [2.5, 97.5]))


def quantile_sensitivity(F, E, firs):
    """Median residual when the cell floor uses quantile q (eligibility and FIR fallback as in
    the main specification)."""
    out = {}
    elig = cell_floors(F).eligible
    b = F[(F.period == "baseline") & (F.clearance_km > 0) & F.cell.notna()]
    g = b.groupby("cell").clearance_km
    for q in (0, 0.01, 0.02, 0.05, 0.10):
        fl = (g.min() if q == 0 else g.quantile(q))[elig]
        fx = E.cell.map(fl).fillna(E.nearest_fir.map(firs.floor))
        out[q] = float(np.nanmedian(E.clearance_km - fx))
    return out


def eligibility_sensitivity(F, boundary_geom_xy, firs):
    rows = []
    for deg, near, mn in ((1.5, 50, 20), (1.5, 50, 30), (1.5, 50, 50), (1.5, 30, 20), (1.5, 30, 40),
                          (1.0, 50, 20), (2.5, 50, 20)):
        Fx = F if deg == C.CELL_DEG else anchor_cells(F.drop(columns=["cell"]), boundary_geom_xy, deg)
        E, _ = residuals(Fx, cell_floors(Fx, near, mn), firs)
        rows.append(dict(cell_size_deg=deg, near_boundary_km=near, min_near_boundary_flights=mn,
                         main=(deg, near, mn) == (1.5, 50, 20),
                         share_cell_floor=(E.floor_source == "cell").mean(),
                         median_residual_km=E.residual.median()))
    return pd.DataFrame(rows)


def validation_cells(F, cells):
    """Eligible cells with >= 20 closure-exposed event flights: baseline floor and event P10."""
    E = exposed_event(F)
    E = E[E.cell.isin(cells.index[cells.eligible])]
    g = E.groupby("cell").clearance_km
    t = pd.DataFrame({"floor": cells.floor, "n_event": g.size(), "p10": g.quantile(0.10)}).dropna()
    return t[t.n_event >= C.VALIDATION_MIN_EVENT]
