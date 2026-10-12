"""Figures 5 to 8 of the paper from the outputs of run_all.py.

Usage: python make_figures.py --run results/main --out figures/
Distances are drawn in the unit of the run's tables (NM by default).
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe                 # noqa: E402
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
# label offsets (points) used in the paper's Fig. 8; other coloured operators get (4, 0)
LABEL_OFF = {"RBG": (-4, 0, "right", "center"), "NMA": (-4, 0, "right", "center"),
             "VIR": (4, 0, "left", "center"), "SWR": (4, 0, "left", "center"),
             "MSR": (0, 3.5, "center", "bottom"), "KLM": (3, -4, "left", "center"),
             "AFR": (-3, -3.5, "center", "top"), "THY": (1, 3, "left", "bottom"),
             "WMT": (-4, 1, "right", "center"), "WZZ": (-3, -2.5, "right", "top"),
             "SXS": (-4, 1, "right", "center"), "AEE": (4, -2, "left", "top")}
LABEL_TEXT = {"AEE": "AEE, FIE"}
SKIP_LABEL = {"FIE"}                                # shares the AEE label in the paper


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
        pe_ = P[(P.episode == ev) & (P.markets == "closure-exposed") & (P[d] <= xmax + 0.05)]
        pu = P[(P.episode == ev) & (P.markets == "unexposed") & (P[d] <= xmax + 0.05)]
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

    fig, ax = plt.subplots(figsize=(3.45, 2.55))
    ax.axvspan(20 * f, 200 * f, color="#efefef", zorder=0)
    dd = np.linspace(0, xmax, 301)
    handles, labels = [], []
    fmt = ".1f" if tag == "nm" else "g"
    for ev in EPISODES:
        col, mk = COLOR[ev], MARKER[ev]
        q = P[(P.episode == ev) & (P.markets == "closure-exposed") & (P[d] <= xmax + 0.05) & (P.h_hat > 0)]
        ax.vlines(q[d], np.maximum(q.ci_lo, 1e-4), q.ci_hi, color=col, lw=0.6, alpha=0.5, zorder=3)
        ax.plot(q[d], q.h_hat, mk, color=col, ms=3.0, zorder=4)
        sig = Fi.loc[ev, f"sigma_{tag}"]
        ax.plot(dd, np.exp(-dd / sig), "--", color=col, lw=0.9, alpha=0.85)
        handles.append((Line2D([], [], ls="--", color=col, lw=0.9), Line2D([], [], ls="", marker=mk, color=col, ms=3)))
        labels.append(fr"{ev}, exponential $\sigma$={sig:{fmt}} {lab}")
    qb = Fi[f"xi_bar_{tag}"].iloc[0]
    ax.plot([0, qb, qb, xmax], [1, 1, 1e-6, 1e-6], color="#9a9a9a", lw=1.0, alpha=0.7)
    handles.append(Line2D([], [], color="#9a9a9a", lw=1.0, alpha=0.7))
    labels.append(fr"hard buffer $\bar\xi$={qb:{fmt}} {lab}")
    ax.set_yscale("log"); ax.set_xlim(0, xmax); ax.set_ylim(2e-3, 1.5)
    ax.set_xlabel(f"distance $d$ ({lab})"); ax.set_ylabel(r"$\widehat h(d)$ (log scale)")
    ax.legend(handles, labels, loc="lower left", frameon=False, handlelength=1.8, borderpad=0.2,
              labelspacing=0.2, handler_map={tuple: HandlerTuple(ndivide=None, pad=0)})
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
    T = {ev: pd.read_csv(run / "tables" / f"S4_1_operator_effects_{ev}.csv").set_index("operator")
         for ev in EPISODES}
    tag, lab, _ = unit_of(T["2025"].columns, "delta")
    j = T["2025"].join(T["2026"], lsuffix="_25", rsuffix="_26", how="inner")
    x, y = j[f"delta_{tag}_25"], j[f"delta_{tag}_26"]
    rc = pd.read_csv(run / "tables" / "S4_3_rank_correlations.csv").iloc[0]   # computed in km
    rho, p, n = rc.rho_s, rc.p_value, len(j)
    u25, u26 = (x.rank() - 0.5) / n, (y.rank() - 0.5) / n
    cnt = j[["baseline_flights_25", "event_flights_25", "baseline_flights_26", "event_flights_26"]]
    elig = cnt.min(1) >= 50
    cls = pd.Series("other", index=j.index)
    cls[(u25 > 2 / 3) & (u26 > 2 / 3) & elig] = "high"
    cls[(u25 < 1 / 3) & (u26 < 1 / 3) & elig] = "low"
    size = 9 + 16 * np.sqrt(cnt.sum(1) / cnt.sum(1).max())
    fig, ax = plt.subplots(figsize=(3.45, 3.2))
    ax.axhline(0, color="#999", lw=0.4, zorder=0); ax.axvline(0, color="#999", lw=0.4, zorder=0)
    for c in ("other", "low", "high"):
        m = cls == c
        ax.scatter(x[m], y[m], s=size[m], facecolor=CLASS_COLOR[c], edgecolor="white", lw=0.4,
                   alpha=0.85 if c == "other" else 0.95, zorder=3 if c == "other" else 4)
    halo = [pe.withStroke(linewidth=1.6, foreground="white")]
    for op in cls.index[cls != "other"]:
        if op in SKIP_LABEL:
            continue
        dx, dy, ha, va = LABEL_OFF.get(op, (4, 0, "left", "center"))
        ax.annotate(LABEL_TEXT.get(op, op), (x[op], y[op]), xytext=(dx, dy), textcoords="offset points",
                    ha=ha, va=va, fontsize=5.6, color=CLASS_COLOR[cls[op]], fontweight="bold",
                    path_effects=halo, zorder=6)
    lim = (-27, 49) if tag == "nm" else (-50, 90)
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
    ax.set_xlabel(rf"$\widehat\delta_a^{{\,2025}}$ ({lab})"); ax.set_ylabel(rf"$\widehat\delta_a^{{\,2026}}$ ({lab})")
    mant, ex = f"{p:.1e}".split("e")
    ax.text(0.03, 0.97, fr"$\rho_S$ = {rho:.3f}" + "\n" + fr"$p$ = {mant}$\times10^{{{int(ex)}}}$" + "\n"
            + fr"$n$ = {n}", transform=ax.transAxes, ha="left", va="top", fontsize=6.5)
    handles = [Line2D([], [], ls="", marker="o", ms=4.5, mfc=CLASS_COLOR[k], mec="none", label=t) for k, t in
               (("high", f"upper tercile in both ({(cls == 'high').sum()})"),
                ("low", f"lower tercile in both ({(cls == 'low').sum()})"), ("other", "other operators"))]
    ax.legend(handles=handles, loc="lower right", fontsize=5.8, frameon=False, handletextpad=0.2)
    finish(ax, fig, out, "fig7_operator_scatter")


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
