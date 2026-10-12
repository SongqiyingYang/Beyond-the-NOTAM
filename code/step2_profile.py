"""Step 2: revealed avoidance profile (Section III-A).

  profile_with_ci      h_hat(d) of Eq. (9) with market-clustered bootstrap intervals (B = 2000)
  benchmark_fits       hard buffer (Eq. 12) and exponential (Eq. 13) fitted on 13 thresholds,
                       20-200 km, weighted least squares on the log scale
  local_slopes         beta(d) of Eq. (14) over 20-60, 60-150 and 150-300 km
  half_decay           d_1/2 of Eq. (15), anchored at 50 km
  exposure_definitions the seven exposure definitions (Eq. 16 and continuing services),
                       each with a market-clustered bootstrap interval for d_1/2 (B = 1500)

The bootstrap resamples OD markets with replacement. Random streams: the profile bootstrap
draws first for the unexposed and then for the closure-exposed markets from one generator
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
    market = F.groupby(unit).od.first().reindex(us).values
    return np.array(us), cb, ce, nb, ne, market


def h_hat(cb, ce, nb, ne, rows=None):
    """Eq. (9): 1 - S_E(d) / S_B(d), with S_E standardised to the baseline market mix (Eq. 8)."""
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
def fit_shape(d, g, se, shape):
    """Grid search of the one-parameter shape minimising the weighted squared log error;
    returns (parameter, unweighted log-RMSE over the same nodes)."""
    grid = C.SIGMA_GRID if shape == "exponential" else C.XI_BAR_GRID
    curves = (np.array([np.exp(-d / s) for s in grid]) if shape == "exponential"
              else np.array([(d <= q).astype(float) for q in grid]))
    lg = np.log(np.maximum(g, C.LOG_FLOOR))
    lc = np.log(np.maximum(curves, C.LOG_FLOOR))
    w = (g / np.maximum(se, 1e-9)) ** 2                  # delta method: 1 / var(log h_hat)
    j = int(np.argmin((((lc - lg) ** 2) * w).sum(1)))
    return float(grid[j]), float(np.sqrt(((lc[j] - lg) ** 2).mean()))


def benchmark_fits(P):
    win = P[P.d_km.isin(C.FIT_NODES)].sort_values("d_km")
    d, g, se = win.d_km.values, win.h_hat.values, win.se.values
    sigma, r_exp = fit_shape(d, g, se, "exponential")
    xi_bar, r_hard = fit_shape(d, g, se, "hard")
    return dict(fit_nodes=len(win), sigma_km=sigma, log_rmse_exponential=r_exp,
                xi_bar_km=xi_bar, log_rmse_hard_buffer=r_hard, ratio_hard_to_exponential=r_hard / r_exp)


# ----------------------------------------------------------------------------- slopes and scale
def local_slopes(d, g, se):
    """Eq. (14): weighted least-squares slope of log h_hat on log d within each window."""
    out = []
    for lo, hi in C.SLOPE_WINDOWS:
        m = (d >= lo) & (d <= hi) & (g > 0)
        if m.sum() < 3:
            out.append(np.nan); continue
        x, y = np.log(d[m]), np.log(g[m])
        w = (g[m] / np.maximum(se[m], 1e-9)) ** 2
        xm, ym = np.average(x, weights=w), np.average(y, weights=w)
        out.append(float((w * (x - xm) * (y - ym)).sum() / (w * (x - xm) ** 2).sum()))
    return np.array(out)


def slopes_with_ci(P, boot):
    d, se = P.d_km.values, P.se.values
    point = local_slopes(d, P.h_hat.values, se)
    bs = np.array([local_slopes(d, b, se) for b in boot])
    lo, hi = np.nanpercentile(bs, [2.5, 97.5], axis=0)
    return point, lo, hi


def half_decay(d, g, anchor=C.D_HALF_ANCHOR):
    """Eq. (15): first distance beyond the anchor at which h_hat falls to half its anchor value
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
DEFINITIONS = [("market, r_o > 0 (main)", "od", 0.0), ("market, r_o >= 0.25", "od", 0.25),
               ("market, r_o >= 0.50", "od", 0.50), ("market, r_o >= 0.75", "od", 0.75),
               ("market, r_o = 1", "od", 1.0), ("continuing service, r > 0", "svc", 0.0),
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
    us, cb, ce, nb, ne, market = count_matrices(sub, unit)
    g = h_hat(cb, ce, nb, ne)
    groups = pd.Series(range(len(us))).groupby(pd.Series(market)).apply(list)
    boot = np.empty((B, len(C.DGRID)))
    for b in range(B):
        pick = rng.integers(0, len(groups), len(groups))
        boot[b] = h_hat(cb, ce, nb, ne, np.concatenate([groups.iloc[i] for i in pick]))
    se = np.nanstd(boot, axis=0)
    dh = np.array([half_decay(C.DGRID, b) for b in boot]); dh = dh[np.isfinite(dh)]
    m = np.isfinite(g) & (C.DGRID >= 20) & (C.DGRID <= 200) & (g > 0)
    P = pd.DataFrame(dict(d_km=C.DGRID, h_hat=g, se=se))
    fits = benchmark_fits(P)
    return dict(units=len(us), event_flights=int((sub.period == "event").sum()),
                monotone=bool(np.all(np.diff(g[m]) <= 1e-9)), d_half_km=half_decay(C.DGRID, g),
                d_half_lo_km=float(np.percentile(dh, 2.5)), d_half_hi_km=float(np.percentile(dh, 97.5)),
                ratio_hard_to_exponential=fits["ratio_hard_to_exponential"])


def exposure_definitions(samples, rng):
    """samples: {episode: analysis sample with markets in both periods and, for the
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
