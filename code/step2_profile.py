"""Step 2: revealed avoidance profile (Section III-A).

  profile_with_ci      h_hat(d) of Eq. (8) with OD-pair-clustered bootstrap intervals (B = 2000)
  benchmark_fits       two-parameter hard buffer A 1{d <= xi_bar} (Eq. 11) and exponential
                       A exp(-d/sigma) (Eq. 12) fitted on 13 thresholds, 20-200 km, by weighted
                       least squares on the original scale with inverse-variance weights 1/se^2;
                       A is estimated separately for each benchmark
  half_decay           d_1/2 of Eq. (13), anchored at 50 km
  exposure_definitions the seven exposure definitions (Eq. 14 and continuing services),
                       each with an OD-pair-clustered bootstrap interval for d_1/2 (B = 1500)

The bootstrap resamples OD pairs with replacement. Random streams: the profile bootstrap
draws first for the unexposed and then for the closure-exposed OD pairs from one generator
seeded with SEED_PROFILE; the exposure definitions draw, in the order 2026 then 2025 and the
order of DEFINITIONS, from one generator seeded with SEED_DHALF.
"""
import numpy as np
import pandas as pd

import common as C


# ----------------------------------------------------------------------------- profile
def count_matrices(F, unit="od", grid=C.DGRID):
    """Per unit and period: number of flights and number with clearance below each d."""
    us = sorted(F[unit].unique())
    ui = {u: i for i, u in enumerate(us)}
    cb = np.zeros((len(us), len(grid))); ce = np.zeros_like(cb)
    nb = np.zeros(len(us)); ne = np.zeros(len(us))
    for per, cc, nn in (("baseline", cb, nb), ("event", ce, ne)):
        s = F[F.period == per]
        idx = s[unit].map(ui).values
        clr = s.clearance_km.values
        np.add.at(nn, idx, 1.0)
        for k, d in enumerate(grid):
            np.add.at(cc, (idx, k), (clr < d).astype(float))
    od_of_unit = F.groupby(unit).od.first().reindex(us).values
    return np.array(us), cb, ce, nb, ne, od_of_unit


def h_hat(cb, ce, nb, ne, rows=None):
    """Eq. (8): 1 - S_E(d) / S_B(d), with S_E standardised to the baseline OD-pair mix (Eq. 7)."""
    if rows is not None:
        cb, ce, nb, ne = cb[rows], ce[rows], nb[rows], ne[rows]
    ok = (nb > 0) & (ne > 0)
    if ok.sum() == 0:
        return np.full(cb.shape[1], np.nan)
    Sb = cb[ok].sum(0) / nb[ok].sum()
    w = nb[ok]
    Se = (w[:, None] * (ce[ok] / ne[ok][:, None])).sum(0) / w.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(Sb > 0, 1.0 - Se / Sb, np.nan)


def profile_with_ci(F, rng, B=C.B_PROFILE):
    """h_hat on DGRID with percentile intervals and bootstrap standard errors."""
    us, cb, ce, nb, ne, _ = count_matrices(F)
    M = len(us)
    g = h_hat(cb, ce, nb, ne)
    boot = np.array([h_hat(cb, ce, nb, ne, rng.integers(0, M, M)) for _ in range(B)])
    lo, hi = np.nanpercentile(boot, [2.5, 97.5], axis=0)
    ev = F[F.period == "event"].clearance_km.values
    bs = F[F.period == "baseline"].clearance_km.values
    P = pd.DataFrame(dict(d_km=C.DGRID, h_hat=g, ci_lo=lo, ci_hi=hi, se=np.nanstd(boot, axis=0),
                          event_flights_within_d=[(ev < d).sum() for d in C.DGRID],
                          baseline_flights_within_d=[(bs < d).sum() for d in C.DGRID]))
    return P, boot


# ----------------------------------------------------------------------------- benchmarks
def _wls_amplitude(g, shapes, w):
    """Closed-form weighted least-squares amplitude for each candidate shape (rows of `shapes`)."""
    A = (w * g * shapes).sum(1) / np.maximum((w * shapes ** 2).sum(1), 1e-12)
    return A, g - A[:, None] * shapes


def fit_shape(d, g, se, shape):
    """Two-parameter benchmark fitted on the original scale (Eqs. 11-12).

    exponential: A exp(-d/sigma), sigma on SIGMA_GRID;
    hard:        A 1{d <= xi_bar}, xi_bar on a grid of XI_BAR_STEP across the fitting window.
    The hard-buffer loss is constant between two adjacent thresholds, so the edge is identified
    only as an interval [xi_lo, xi_hi), where xi_hi is the next threshold beyond the optimum.
    Returns a dict with the parameter(s), the amplitude and the weighted RMSE over the nodes."""
    w = 1.0 / np.maximum(se, 1e-9) ** 2
    if shape == "exponential":
        grid = C.SIGMA_GRID
        shapes = np.exp(-d[None, :] / grid[:, None])
    else:
        grid = np.arange(d.min(), d.max() + C.XI_BAR_STEP, C.XI_BAR_STEP)
        shapes = (d[None, :] <= grid[:, None]).astype(float)
    A, res = _wls_amplitude(g, shapes, w)
    loss = (w * res ** 2).sum(1)
    j = int(np.argmin(loss))
    out = dict(A=float(A[j]), wrmse=float(np.sqrt(loss[j] / len(d))))
    if shape == "exponential":
        out["sigma_km"] = float(grid[j])
    else:
        tie = grid[np.isclose(loss, loss[j])]
        nxt = d[d > tie.max()]
        out["xi_bar_lo_km"] = float(tie.min())
        out["xi_bar_hi_km"] = float(nxt.min() if len(nxt) else tie.max())
    return out


def benchmark_fits(P):
    win = P[P.d_km.isin(C.FIT_NODES)].sort_values("d_km")
    d, g, se = win.d_km.values, win.h_hat.values, win.se.values
    e, h = fit_shape(d, g, se, "exponential"), fit_shape(d, g, se, "hard")
    return dict(fit_nodes=len(win),
                A_exponential=e["A"], sigma_km=e["sigma_km"], wrmse_exponential=e["wrmse"],
                A_hard_buffer=h["A"], xi_bar_lo_km=h["xi_bar_lo_km"], xi_bar_hi_km=h["xi_bar_hi_km"],
                wrmse_hard_buffer=h["wrmse"], ratio_hard_to_exponential=h["wrmse"] / e["wrmse"])


# ----------------------------------------------------------------------------- spatial scale
def half_decay(d, g, anchor=C.D_HALF_ANCHOR):
    """Eq. (13): first distance beyond the anchor at which h_hat falls to half its anchor value
    (linear interpolation between nodes)."""
    tgt = 0.5 * np.interp(anchor, d, g)
    below = np.where((d > anchor) & (g <= tgt))[0]
    if len(below) == 0:
        return np.nan
    j = below[0]
    if j == 0:
        return float(d[0])
    x0, x1, y0, y1 = d[j - 1], d[j], g[j - 1], g[j]
    return float(x0 + (tgt - y0) * (x1 - x0) / (y1 - y0)) if y1 != y0 else float(x1)


# ----------------------------------------------------------------------------- exposure definitions
DEFINITIONS = [("OD pair, r_o > 0 (main)", "od", 0.0), ("OD pair, r_o >= 0.25", "od", 0.25),
               ("OD pair, r_o >= 0.50", "od", 0.50), ("OD pair, r_o >= 0.75", "od", 0.75),
               ("OD pair, r_o = 1", "od", 1.0), ("continuing service, r > 0", "svc", 0.0),
               ("continuing service, r = 1", "svc", 1.0)]


def _select(F, unit, tau):
    if unit == "od":
        r = F.od.map(C.penetration_rate(F))
    else:
        b = F[F.period == "baseline"]
        r = F.svc.map(b.assign(z=b.clearance_km <= 0).groupby("svc").z.mean())
    if tau == 0:
        return F[r > 0]
    return F[r >= (0.999 if tau == 1.0 else tau)]           # r = 1 allows for float rounding


def evaluate_definition(sub, unit, rng, B=C.B_DHALF):
    us, cb, ce, nb, ne, od_of_unit = count_matrices(sub, unit)
    g = h_hat(cb, ce, nb, ne)
    groups = pd.Series(range(len(us))).groupby(pd.Series(od_of_unit)).apply(list)
    boot = np.empty((B, len(C.DGRID)))
    for b in range(B):
        pick = rng.integers(0, len(groups), len(groups))
        boot[b] = h_hat(cb, ce, nb, ne, np.concatenate([groups.iloc[i] for i in pick]))
    dh = np.array([half_decay(C.DGRID, b) for b in boot]); dh = dh[np.isfinite(dh)]
    m = np.isfinite(g) & (C.DGRID >= 20) & (C.DGRID <= 200) & (g > 0)
    return dict(units=len(us), event_flights=int((sub.period == "event").sum()),
                monotone=bool(np.all(np.diff(g[m]) <= 1e-9)), d_half_km=half_decay(C.DGRID, g),
                d_half_lo_km=float(np.percentile(dh, 2.5)), d_half_hi_km=float(np.percentile(dh, 97.5)))


def exposure_definitions(samples, rng):
    """samples: {episode: analysis sample with OD pairs in both periods and, for the
    service-level rows, the flight-number column 'flight'}."""
    rows = []
    for ep in ("2026", "2025"):
        F = samples[ep]
        for name, unit, tau in DEFINITIONS:
            if unit == "svc":
                if "flight" not in F.columns or F.flight.isna().all():
                    continue
                S = F[F.flight.notna() & (F.flight.astype(str).str.len() > 2)].copy()
                S["svc"] = S.operator.astype(str) + "|" + S.flight.astype(str)
                alive = S.groupby("svc").period.nunique() == 2
                S = S[S.svc.isin(alive[alive].index)]
                sub = _select(S, "svc", tau)
            else:
                sub = _select(F, "od", tau)
            rows.append(dict(episode=ep, definition=name, **evaluate_definition(sub, unit, rng)))
    return pd.DataFrame(rows)
