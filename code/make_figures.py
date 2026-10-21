"""Figures 5 to 8 of the paper from the outputs of run_all.py (Figs. 5-7) and
step4_operators.py (Fig. 8), written to the same --run directory.

Usage: python make_figures.py --run results/main --out figures/
Distances are drawn in the unit of the run's tables (NM by default).
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
import numpy as np                                  # noqa: E402
import pandas as pd                                 # noqa: E402
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D                 # noqa: E402

import common as C                                  # noqa: E402

COLOR = {"2026": "#d62728", "2025": "#1f77b4"}
MARKER = {"2026": "s", "2025": "o"}
CLASS_COLOR = {"high": "#6a3d9a", "low": "#1b9e9e", "other": "#c8c8c8"}
EPISODES = ("2025", "2026")


def style():
    matplotlib.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.labelsize": 6.8,
                                "legend.fontsize": 6.0, "axes.unicode_minus": False,
                                "pdf.fonttype": 42, "ps.fonttype": 42})


def finish(ax, fig, out, stem):
    ax.tick_params(labelsize=6.3)
    ax.grid(alpha=0.2, which="major")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{stem}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {stem}.pdf")


def unit_of(columns, stem):
    tag = [c for c in columns if c.startswith(stem + "_")][0].split("_")[-1]
    return tag, ("NM" if tag == "nm" else "km"), (C.NM_PER_KM if tag == "nm" else 1.0)


def fig5(run, out):
    P = pd.read_csv(run / "tables" / "S2_1_avoidance_profile.csv", dtype={"episode": str})
    Fi = pd.read_csv(run / "tables" / "S2_2_benchmark_fits.csv", dtype={"episode": str}).set_index("episode")
    tag, lab, f = unit_of(P.columns, "d")
    d = f"d_{tag}"
    xmax = 300 * f
    fig, ax = plt.subplots(figsize=(3.45, 2.55))
    for ev in EPISODES:
        col, mk = COLOR[ev], MARKER[ev]
        pe_ = P[(P.episode == ev) & (P.od_pairs == "closure-exposed") & (P[d] <= xmax + 0.05)]
        pu = P[(P.episode == ev) & (P.od_pairs == "unexposed") & (P[d] <= xmax + 0.05)]
        ax.plot(pu[d], pu.h_hat, "--", color=col, lw=0.8, marker=mk, ms=2.4, mfc="white", mec=col, mew=0.6,
                label=f"{ev} unexposed")
        ax.fill_between(pe_[d], pe_.ci_lo, pe_.ci_hi, color=col, alpha=0.18, lw=0)
        ax.plot(pe_[d], pe_.h_hat, "-", color=col, lw=1.3, marker=mk, ms=2.8, label=f"{ev} closure-exposed")
        dh = Fi.loc[ev, f"d_half_{tag}"]
        ax.axvline(dh, color=col, lw=0.7, ls=":")
        ax.text(dh + 3 * f, 0.30 if ev == "2026" else 0.40, fr"$d_{{1/2}}^{{{ev}}}$ = {dh:.1f} {lab}",
                color=col, fontsize=6, va="center",
                bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none", alpha=0.85))
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlim(0, xmax); ax.set_ylim(-0.1, 0.8)
    ax.set_xlabel(f"distance $d$ ({lab})"); ax.set_ylabel(r"$\widehat h(d)$")
    ax.legend(loc="upper right", frameon=False, handlelength=2.3, labelspacing=0.3)
    finish(ax, fig, out, "fig5_profile")

    # (b) two-parameter benchmarks on the original scale (Eqs. 11-12)
    fig, ax = plt.subplots(figsize=(3.45, 2.55))
    ax.axvspan(20 * f, 200 * f, color="#f2f2f2", zorder=0)
    dd = np.linspace(0, xmax, 601)
    handles, labels, edges = [], [], []
    for ev in EPISODES:
        col, mk = COLOR[ev], MARKER[ev]
        q = P[(P.episode == ev) & (P.od_pairs == "closure-exposed") & (P[d] <= xmax + 0.05)]
        ax.vlines(q[d], q.ci_lo, q.ci_hi, color=col, lw=0.6, alpha=0.5, zorder=3)
        ax.plot(q[d], q.h_hat, mk, color=col, ms=3.0, zorder=4)
        r = Fi.loc[ev]
        sig, xlo, xhi = r[f"sigma_{tag}"], r[f"xi_bar_lo_{tag}"], r[f"xi_bar_hi_{tag}"]
        edges.append((xlo, xhi))
        ax.plot(dd, r.A_exponential * np.exp(-dd / sig), "--", color=col, lw=0.9, alpha=0.9, zorder=5)
        xm = (xlo + xhi) / 2                      # edge drawn at the middle of its interval
        ax.plot([0, xm, xm, xmax], [r.A_hard_buffer, r.A_hard_buffer, 0, 0], "-", color=col, lw=0.9,
                alpha=0.6, zorder=5)
        handles += [(Line2D([], [], ls="", marker=mk, color=col, ms=3.0), Line2D([], [], ls="--", color=col, lw=0.9)),
                    Line2D([], [], ls="-", color=col, lw=0.9, alpha=0.6)]
        labels += [fr"{ev} observed, exponential ($\sigma$={sig:.1f} {lab})", f"{ev} hard buffer"]
    ax.axvspan(min(e[0] for e in edges), max(e[1] for e in edges), color="#d9d9d9", alpha=0.6, lw=0, zorder=1)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlim(0, xmax); ax.set_ylim(-0.05, 0.85)
    ax.set_xlabel(f"distance $d$ ({lab})"); ax.set_ylabel(r"$\widehat h(d)$")
    ax.legend(handles, labels, loc="upper right", fontsize=5.6, frameon=False, handlelength=1.8, borderpad=0.2,
              labelspacing=0.25, handler_map={tuple: HandlerTuple(ndivide=None, pad=0)})
    finish(ax, fig, out, "fig5_benchmark")
    return tag, lab, f


def fig6_7(run, out, lab, f):
    fig, ax = plt.subplots(figsize=(3.33, 2.55))
    for ev in EPISODES:
        v = pd.read_csv(run / "private" / f"validation_cells_{ev}.csv")
        col, mk = COLOR[ev], MARKER[ev]
        ax.scatter(v.floor * f, v.p10 * f, s=4 + v.n_event / 11, marker=mk, facecolor=col, edgecolor="white",
                   lw=0.4, alpha=0.8, zorder=4, label=f"{ev} ({len(v)} cells)")
    lo, hi = (0.4, 80) if lab == "NM" else (0.7, 150)
    ax.plot([lo, hi], [lo, hi], "--", color="k", lw=0.7, zorder=2)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    ax.set_xlabel(rf"baseline structural floor $\widehat\phi(a)$, P5 ({lab})")
    ax.set_ylabel(f"event-period clearance, P10 ({lab})")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=6.2, frameon=False,
              markerscale=0.8, borderaxespad=0, handletextpad=0.2)
    finish(ax, fig, out, "fig6_floor_attainability")

    R = pd.read_parquet(run / "private" / "residuals.parquet")
    fig, ax = plt.subplots(figsize=(3.45, 2.6))
    for ev in EPISODES:
        e = np.sort(R[R.episode == ev].residual.values) * f
        y = np.arange(1, len(e) + 1) / len(e)
        col = COLOR[ev]
        ax.step(e, y, where="post", color=col, lw=1.2, label=f"{ev} (n={len(e):,})")
        m = np.median(e)
        ax.plot([m, m], [0, 0.5], ":", color=col, lw=0.8)
        ax.plot(m, 0.5, "o", color=col, ms=3)
        ax.text(m + (4 if ev == "2025" else 14) * f, 0.44 if ev == "2025" else 0.54,
                f"median = {m:.1f} {lab} ($q$ = 0.05)", color=col, fontsize=6, va="center")
    ax.axvline(0, color="k", lw=0.7)
    ax.set_xlim(*((-25, 110) if lab == "NM" else (-50, 200))); ax.set_ylim(0, 1)
    ax.set_xlabel(rf"residual standoff $\widehat e_f$ ({lab})")
    ax.set_ylabel(r"$P(\widehat e_f \leq e)$")
    ax.legend(loc="lower right", fontsize=6.3, frameon=False)
    finish(ax, fig, out, "fig6_residual_ecdf")


def fig8(run, out):
    """Operator-relative standoff R_k (Eq. 15), 2025 against 2026, main specification."""
    T = pd.read_csv(run / "tables" / "S4_1_operator_relative_standoff.csv").set_index("operator")
    tag = [c for c in T.columns if c.startswith("R_")][0].split("_")[1]
    lab = "NM" if tag == "nm" else "km"
    j = T                                             # S4_1 holds the main panel only
    x, y = j[f"R_{tag}_2025"], j[f"R_{tag}_2026"]
    R = pd.read_csv(run / "tables" / "S4_2_rank_correlations.csv")
    rc = R[(R.reference == "all baseline crossings") & (R.min_comparable_flights == C.OPERATOR_MIN_MAIN)
           & (R.arc_fraction == C.ARC_FRACTION_MAIN)].iloc[0]
    cls = j.tercile_both.fillna("").map({"upper": "high", "lower": "low"}).fillna("other")
    nt = j.comparable_flights_2025 + j.comparable_flights_2026
    size = 12 + 40 * np.sqrt(nt / nt.max())
    fig, ax = plt.subplots(figsize=(3.6, 3.4))
    ax.axhline(0, color="#999", lw=0.5, zorder=0); ax.axvline(0, color="#999", lw=0.5, zorder=0)
    for c in ("other", "low", "high"):
        m = cls == c
        ax.scatter(x[m], y[m], s=size[m], color=CLASS_COLOR[c], lw=0.4, edgecolor="white",
                   zorder=3 if c != "other" else 2)
    for op in cls.index[cls != "other"]:
        ax.annotate(op, (x[op], y[op]), xytext=(4, 2), textcoords="offset points", fontsize=6.5,
                    color=CLASS_COLOR[cls[op]], fontweight="bold")
    ax.text(0.02, 0.98, f"Spearman ρ = {rc.rho_s:.2f}\ntwo-sided permutation p = {rc.p_two_sided:.2f}"
            f"\nn = {int(rc.operators)} operators", transform=ax.transAxes, ha="left", va="top", fontsize=7)
    pct = int(C.ARC_FRACTION_MAIN * 100)
    tl = f"Distance relative to other operators\nnearest {pct}% of the detour arc, {{}} ({lab})"
    ax.set_xlabel(tl.format(2025), fontsize=7.5); ax.set_ylabel(tl.format(2026), fontsize=7.5)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=CLASS_COLOR[k], ms=5, label=t) for k, t in
                       (("high", f"farther than peers in both ({(cls == 'high').sum()})"),
                        ("low", f"closer than peers in both ({(cls == 'low').sum()})"),
                        ("other", "other operators"))],
              fontsize=6, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.83))
    finish(ax, fig, out, "fig8_operator_relative_standoff")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="output directory of run_all.py")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    run, out = Path(a.run), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    style()
    _, lab, f = fig5(run, out)
    fig6_7(run, out, lab, f)
    fig8(run, out)


if __name__ == "__main__":
    main()
