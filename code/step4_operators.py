"""Step 4: operator-associated heterogeneity and cross-event persistence (Section III-C).

  operator_effects   delta_a from the market-operator fixed-effects model (Eq. 19), with
                     standard errors clustered by OD market and the joint F-test of equal
                     operator effects
  rank_correlations  Spearman correlation of delta_a across episodes (Eq. 20) on all operators
                     observed in both, and on operators with >= 100 or >= 200 flights in each
                     period of each episode
  exposed_panel      operators with >= 20 closure-exposed flights in each period of each
                     episode; median raw event-period clearance and median residual standoff
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr

import common as C


def operator_effects(F):
    """F: analysis sample of one episode. Returns (table of delta_a, summary dict)."""
    b = F[F.period == "baseline"].groupby("od").clearance_km.agg(["median", "size"])
    markets = b[(b["median"] < C.MARKET_MAX_MEDIAN_KM) & (b["size"] >= C.MARKET_MIN_BASELINE)].index
    R = F[F.od.isin(markets)].copy()
    cnt = R.pivot_table(index="operator", columns="period", values="fr24_id", aggfunc="count").fillna(0)
    keep = cnt[(cnt.get("baseline", 0) >= C.OPERATOR_MIN_FLIGHTS)
               & (cnt.get("event", 0) >= C.OPERATOR_MIN_FLIGHTS)].index
    R = R[R.operator.isin(keep)].copy()
    R["cell"] = R.od.astype(str) + "|" + R.operator.astype(str)
    both = R.groupby("cell").period.nunique() == 2
    R = R[R.cell.isin(both[both].index)].copy()

    # within-cell transformation removes alpha_{o x a}; delta_a are the coefficients of
    # E_f x 1{a(f) = a} (Frisch-Waugh-Lovell)
    ops = sorted(R.operator.unique())
    E = (R.period == "event").to_numpy(float)
    X = np.column_stack([E * (R.operator.to_numpy() == a) for a in ops])
    cell = R.cell.to_numpy()

    def demean(v):
        return v - pd.Series(v).groupby(cell).transform("mean").to_numpy()

    y = demean(R.clearance_km.to_numpy())
    Xd = np.column_stack([demean(X[:, j]) for j in range(X.shape[1])])
    mod = sm.OLS(y, Xd).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(R.od)[0]})
    n = R.groupby(["operator", "period"]).size().unstack().fillna(0)
    res = pd.DataFrame({"operator": ops, "delta_km": mod.params, "se_km": mod.bse,
                        "ci_lo_km": mod.params - 1.96 * mod.bse, "ci_hi_km": mod.params + 1.96 * mod.bse,
                        "p_value": mod.pvalues,
                        "baseline_flights": n.reindex(ops)["baseline"].to_numpy().astype(int),
                        "event_flights": n.reindex(ops)["event"].to_numpy().astype(int)})
    Rm = np.zeros((len(ops) - 1, len(ops)))
    for i in range(len(ops) - 1):
        Rm[i, i], Rm[i, i + 1] = 1, -1
    ft = mod.f_test(Rm)
    summary = dict(operators=len(ops), flights=len(R), market_operator_cells=R.cell.nunique(),
                   markets=R.od.nunique(), F=float(np.squeeze(ft.fvalue)),
                   p_value=float(np.squeeze(ft.pvalue)),
                   span_delta_km=float(res.delta_km.max() - res.delta_km.min()))
    return res.sort_values("delta_km", ascending=False).reset_index(drop=True), summary


def _rho(x, y):
    r, p = spearmanr(x, y)
    return float(r), float(p), int(len(x))


def rank_correlations(eff, panel=None):
    """eff: {episode: operator_effects table}; delta_a enter rounded to 0.1 km, the precision
    at which they are reported. panel: output of exposed_panel."""
    a = eff["2025"].set_index("operator"); b = eff["2026"].set_index("operator")
    j = a.join(b, lsuffix="_25", rsuffix="_26", how="inner")
    rows = []
    r, p, n = _rho(j.delta_km_25.round(1), j.delta_km_26.round(1))
    rows.append(dict(measure="delta_a", panel="all operators observed in both episodes",
                     rho_s=r, operators=n, p_value=p))
    for th in (100, 200):
        m = j[(j[["baseline_flights_25", "event_flights_25", "baseline_flights_26",
                  "event_flights_26"]] >= th).all(1)]
        r, p, n = _rho(m.delta_km_25.round(1), m.delta_km_26.round(1))
        rows.append(dict(measure="delta_a", panel=f">= {th} flights in each period of each episode",
                         rho_s=r, operators=n, p_value=p))
    if panel is not None:
        pj = panel["2025"].join(panel["2026"], lsuffix="_25", rsuffix="_26", how="inner")
        for col, name in (("raw_median_km", "raw event-period clearance"),
                          ("residual_median_km", "residual standoff")):
            r, p, n = _rho(pj[f"{col}_25"], pj[f"{col}_26"])
            rows.append(dict(measure=name, panel=f"{n}-operator closure-exposed panel",
                             rho_s=r, operators=n, p_value=p))
        k = j.loc[j.index.intersection(pj.index)]
        r, p, n = _rho(k.delta_km_25.round(1), k.delta_km_26.round(1))
        rows.append(dict(measure="delta_a", panel=f"{n}-operator closure-exposed panel",
                         rho_s=r, operators=n, p_value=p))
    return pd.DataFrame(rows)


def exposed_panel(F, E):
    """F: analysis sample with exposure split applied; E: residual table from step 3."""
    X = C.split_exposure(F)["exposed"]
    cnt = X.pivot_table(index="operator", columns="period", values="fr24_id", aggfunc="count").fillna(0)
    keep = cnt[(cnt.get("baseline", 0) >= C.PANEL_MIN_EXPOSED)
               & (cnt.get("event", 0) >= C.PANEL_MIN_EXPOSED)].index
    ev = X[(X.period == "event") & X.operator.isin(keep)]
    out = pd.DataFrame({"raw_median_km": ev.groupby("operator").clearance_km.median(),
                        "residual_median_km": E[E.operator.isin(keep)].groupby("operator").residual.median()})
    return out.dropna()
