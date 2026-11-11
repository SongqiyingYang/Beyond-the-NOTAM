# Beyond the NOTAM: Supplementary Repository

Supplementary material for the paper *Beyond the NOTAM: Revealed Airspace-Avoidance Behaviour During Major Airspace Closures*.

The paper uses before-and-after commercial flight trajectories from two conflict-related airspace-closure episodes in the Middle East (June 2025 and March 2026). It measures the routing response beyond the published closure boundary with three tools: a revealed avoidance profile, a baseline-revealed structural floor, and operator-relative standoff along the detour arcs around the closure. This repository contains two parts:

- **`tables/`**: the results behind the main text and the sensitivity analyses referred to in the paper;
- **`code/`**: the code that produces these tables and Figs. 5–8 from flight trajectories.

All tables are aggregate results. No flight-level data are included (see [Data availability](#data-availability)).

Section, equation, figure, and table numbers refer to the paper.

## Contents

| File | Content | Paper |
|---|---|---|
| **S1. Closure coding** | | |
| `tables/S1_1_closure_calendar.csv` | FIR-level closure calendar of both episodes, with coding, type of evidence, and sources | Sec. II-D, Table I |
| `tables/S1_2_transit_coding_B.csv` | Transit traffic through the two FIRs coded as operationally closed (OMAE 2026, OSTT 2025) | Sec. II-D |
| **S2. Revealed avoidance profile** | | |
| `tables/S2_1_avoidance_profile.csv` | $\widehat h(d)$ with 95% intervals, closure-exposed and unexposed OD pairs, 5.4–486.0 NM | Eq. (8), Fig. 5(a) |
| `tables/S2_2_benchmark_fits.csv` | Hard-buffer and exponential benchmark fits, half-decay distance | Eqs. (11)–(13), Fig. 5(b) |
| `tables/S2_3_exposure_definitions.csv` | Profile summaries under the seven exposure definitions | Eq. (14) |
| **S3. Revealed structural floor** | | |
| `tables/S3_1_boundary_cells.csv` | Every boundary-anchor cell: baseline counts, floor, eligibility, event-period lower tail | Eq. (4), Figs. 4 and 6 |
| `tables/S3_2_fir_fallback_floors.csv` | FIR-level fallback floors | Sec. III-B |
| `tables/S3_3_decomposition.csv` | Median clearance, floor, and residual standoff with 95% intervals | Eq. (5), Fig. 7 |
| `tables/S3_4_floor_quantile.csv` | Median residual standoff for $q \in \{0, 0.01, 0.02, 0.05, 0.10\}$ | Sec. IV-B |
| `tables/S3_5_eligibility_cell_size.csv` | Median residual standoff under alternative cell-eligibility rules and cell sizes | Sec. III-D |
| `tables/S3_6_residual_quantiles.csv` | Quantiles of the residual-standoff distribution | Fig. 7 |
| **S4. Operator-relative standoff** | | |
| `tables/S4_1_operator_relative_standoff.csv` | $R_k$ and number of comparable flights of the 24 operators in both episodes | Eq. (15), Fig. 8 |
| `tables/S4_2_rank_correlations.csv` | Cross-event rank correlations and permutation tests under the operator sensitivity analyses | Secs. III-D, IV-C |
| `tables/S4_3_reference_arcs.csv` | OD pairs with a reference arc and flights with a boundary-distance profile | Sec. III-C |
| **S5. Sensitivity analyses** | | |
| `tables/S5_1_robustness_variants.csv` | Headline quantities under airport trimming and the alternative OMAE codings | Secs. II-D, III-D |
| **Code** | | |
| `code/` | Reproduction code; see [`code/README.md`](code/README.md) | Secs. II–III |

## Conventions

- **Episodes.** The 2025 episode compares the baseline window 6–11 June 2025 with the evaluation window 13–17 June 2025; its closure set is LLLL, OIIX, ORBB, OSTT. The 2026 episode compares 21–27 February 2026 with 2–6 March 2026; its closure set is LLLL, OBBB, OIIX, OKAC, OMAE, ORBB, OSTT, OTDF. Trajectories are on a 60-s grid.
- **OD pairs.** A directional origin–destination (OD) pair is *closure-exposed* when at least one of its baseline flights entered the airspace that subsequently closed ($r_o > 0$, Eqs. (2)–(3)). Otherwise it is *unexposed*. Only OD pairs observed in both periods are used.
- **Distances.** All distances are horizontal distances from the boundary of the FIR-level closure set $S_0$, reported in nautical miles (NM; 1 NM = 1.852 km). Column names carry the unit (`_nm`).
  - The analysis was carried out in kilometers. Its thresholds were set as round values in km and appear here at their converted values: for example, 50 km is 27.0 NM, and the fitting window of 20–200 km is 10.8–108.0 NM.
  - Statistics without a distance unit do not depend on this conversion: profile values, amplitudes, error ratios, $\rho_S$, and $p$.
- **Intervals.** Intervals are 95% percentile intervals from nonparametric bootstraps clustered by OD pair:

  | Quantity | Replications |
  |---|---|
  | Profile | $B = 2000$ |
  | Half-decay distance | $B = 1500$ |
  | Median residual standoff | $B = 1000$ |

- **Operators** are identified by ICAO three-letter airline designators.

## S1. Closure coding

`S1_1_closure_calendar.csv` lists every FIR-level restriction recorded for the two episodes, with its start and end (UTC), its coding, the type of evidence, and the sources. The column `code` gives the coding in machine-readable form:

| `code` | Meaning | In $S_0$ |
|---|---|---|
| `A` | Closed | Yes |
| `B` | Operationally closed: transit confined to designated routes | Yes |
| `night_only` | Closed at night only | No |
| `short_event` | Closed for a few hours | No |
| `sub_fir` | Restriction published below FIR level | No |

A FIR enters $S_0$ on a given day when it is closed for at least half of the 24-h day slice starting at 06:00 UTC.

The type of evidence is one of four:

- **NOTAM:** the NOTAM text is available.
- **Official notice:** a statement by the civil aviation authority.
- **Industry report:** contemporaneous closure tables and reporting, mainly the Flightradar24 blog.
- **Traffic data:** continuity supported only by the absence of transit traffic.

Two FIRs are coded as operationally closed. In both, the transit traffic in `S1_2_transit_coding_B.csv` almost disappeared. Transit traffic here means passenger flights with neither endpoint in the country and at least one trajectory point inside the FIR.

| FIR | Baseline | Evaluation | Change |
|---|---|---|---|
| OMAE (2026) | 553.0 flights per day | 35.0 flights per day | −94% |
| OSTT (2025) | 82.5 flights per day | 0.4 flights per day | −100% |

Two restrictions in the 2026 evaluation window lie outside $S_0$ because they were published below FIR level:

- the partial closure of an area of the Jeddah FIR (OEJD) bordering Iraq and the Persian Gulf, published without lateral limits;
- the closure of the southern sector of the Baku FIR (UBBA) on 5–6 March.

## S2. Revealed avoidance profile

**Profile (`S2_1`).** $\widehat h(d)$ of Eq. (8) is tabulated at every distance node from 5.4 to 486.0 NM, for closure-exposed and unexposed OD pairs in both episodes. The columns `event_flights_within_d` and `baseline_flights_within_d` give the number of flights within distance $d$ in each period. On closure-exposed pairs the profile decreases monotonically in both episodes. On unexposed pairs it does not:

- 2025: $\widehat h(54.0) = 0.014 < \widehat h(108.0) = 0.032$;
- 2026: $\widehat h(27.0) = 0.083 < \widehat h(54.0) = 0.114$.

**Benchmarks (`S2_2`).** The hard buffer $A\,\mathbf 1\{d \le \bar\xi\}$ (Eq. 11) and the exponential $A\exp(-d/\sigma)$ (Eq. 12) are fitted to the closure-exposed profile on the 13 distance thresholds from 10.8 to 108.0 NM.

- **Criterion:** weighted least squares on the original scale, with inverse-variance weights $1/\mathrm{se}^2$, where se is the standard deviation of the profile bootstrap (`S2_1`).
- **Amplitude:** $A$ is estimated separately for each benchmark.
- **Grids:** $\sigma$ from 2.7 to 540.0 NM (5–1000 km) in steps of 0.54 NM (1 km); $\bar\xi$ across the fitting window in the same steps. The hard-buffer loss does not change between two adjacent thresholds, so its edge is identified only as an interval: `xi_bar_lo` is the smallest optimal value and `xi_bar_hi` the next threshold.
- **Fit quality:** the weighted root-mean-square error $\sqrt{\sum_i w_i r_i^2 / 13}$, with the same weights $w_i$ as the fit.

| Episode | $\sigma$ (NM) | $A_{\exp}$ | Hard-buffer edge (NM) | $A_{\text{hard}}$ | Weighted RMSE ratio, hard / exponential | $d_{1/2}$ (NM) |
|---|---|---|---|---|---|---|
| 2025 | 24.8 | 0.86 | 32.4–40.5 | 0.35 | 4.4 | 49.0 [46.9, 53.9] |
| 2026 | 28.6 | 0.80 | 32.4–40.5 | 0.39 | 7.9 | 40.5 [38.3, 54.9] |

**Exposure definitions (`S2_3`).** The table covers seven definitions:

- the main rule $r_o > 0$;
- the stricter thresholds $r_o \ge \tau$ of Eq. (14);
- two definitions at the level of continuing services, where a service is an operator and flight number observed in both periods.

The closure-exposed profile is monotone under all seven. $d_{1/2}$ ranges from 47.0 to 49.1 NM in 2025 and from 38.0 to 43.8 NM in 2026.

## S3. Revealed structural floor

**Cells (`S3_1`).** Each flight is anchored to the point of the closure boundary nearest to its closest approach. Anchors are grouped on a 1.5° × 1.5° grid; `lat_min` to `lon_max` give the bounds of each cell.

- **Eligibility.** A cell is eligible when it has at least 40 positive-clearance baseline flights, of which at least 20 came within 27.0 NM of the boundary. Its floor is the fifth percentile of those clearances. There are 20 eligible cells of 59 in 2025 and 15 of 51 in 2026.
- **Validation cells.** Eligible cells with at least 20 closure-exposed event-period flights are used for the floor-attainability comparison of Fig. 6 (15 cells in 2025, 10 in 2026). For these cells, `event_p10_nm` gives the tenth percentile of event-period clearance.
- **FIR fallback (`S3_2`).** Flights anchored to an ineligible cell use the FIR-level fifth percentile, provided the FIR has at least 100 positive-clearance baseline flights.

**Decomposition (`S3_3`).** Cell floors cover 96.2% (2025) and 90.9% (2026) of the closure-exposed event flights. Four 2026 flights have neither a cell nor a FIR floor and are excluded.

| Episode | Median clearance | Median floor | Median residual standoff [95% CI] | Share of negative residuals |
|---|---|---|---|---|
| 2025 | 20.2 NM | 3.2 NM | 9.8 NM [9.3, 9.9] | 10.3% |
| 2026 | 21.9 NM | 7.6 NM | 9.9 NM [6.8, 13.1] | 17.3% |

**Sensitivity (`S3_4`, `S3_5`).** The median residual standoff is positive under every floor quantile, eligibility rule, and cell size examined, but its value depends on these choices.

| Choice varied | 2025 | 2026 |
|---|---|---|
| Floor quantile $q$ | 6.3–16.3 NM | 1.6–14.8 NM |
| Eligibility rule and cell size | 9.6–9.9 NM | 8.4–12.7 NM |

The eligibility rules require at least 20, 30, 40, or 50 baseline flights within 16.2 or 27.0 NM of the boundary; the cell sizes are 1.0°, 1.5°, and 2.5°.

**Distribution (`S3_6`).** Quantiles of the residual standoff from the 5th to the 95th percentile, corresponding to Fig. 7.

## S4. Operator-relative standoff

**Method (Sec. III-C).** For each closure-exposed OD pair, the baseline flights that entered the closure define an entry point A and an exit point B on a smoothed outline of the closure set (morphological closing with radius 108.0 NM). A and B are the first and last intersections of the trajectory with the outline, with the trajectory joined across data gaps; the reference of the OD pair is the crossing whose (A, B) lies closest to those of the other crossings. A and B split the outline into two detour arcs, sampled every 2.7 NM.

- Each event-period flight is assigned to the arc with the smaller mean distance to its trajectory (its detour direction).
- $d_f(s)$ is the distance from arc position $s$ to the trajectory.
- Flights of the same OD pair and detour direction are compared on a common arc segment: from the first to the last $s$ at which the lowest $d_f(s)$ among them is at most 108.0 NM.
- The flight-level standoff $m_f$ is the mean of the smallest 20% of $d_f(s)$ on that segment.
- $R_k$ (Eq. 15) averages, over operator $k$'s flights, $m_f$ minus the mean $m$ of the other operators in the same OD-pair × detour-direction group. A flight is comparable when its group contains another operator.

Trajectory segments of fewer than three points between data gaps longer than 5 min are treated as position outliers and dropped.

**Operators (`S4_1`).** $R_k$ and the number of comparable flights of the 24 operators with at least 10 comparable flights in each episode. `tercile_both` marks the operators in the upper (farther than peers) or lower (closer than peers) tercile in both episodes, as in Fig. 8.

**Rank correlations (`S4_2`).** Spearman $\rho_S$ of $R_k$ between the episodes for every combination of:

- the reference-arc rule: all baseline crossings (main), or only crossings without a data gap longer than 30 min between A and B (OD pairs without such a crossing are dropped);
- the minimum number of comparable flights per operator in each episode: 10 (main), 15, 20;
- the share of the arc used for $m_f$: 10%, 20% (main), 30%, 50%.

$p$-values come from 2,000 permutations of operator labels within each OD-pair × detour-direction group, which preserve every operator's number of flights per group. `p_two_sided` is $P(|\rho^*| \ge |\rho|)$. `p_family_two_sided` adjusts for the four arc shares by comparing $|\rho|$ with the maximum $|\rho^*|$ over the shares in the same permutation.

| Reference arcs | Min. flights | Operators | $\rho_S$ (10%–50% of the arc) | Two-sided $p$ | Family-adjusted $p$ |
|---|---|---|---|---|---|
| All crossings (main) | 10 | 24 | 0.27–0.35 (20%: 0.32) | 0.13–0.25 (20%: 0.17) | 0.21–0.39 |
| All crossings | 15 or 20 | 17 | 0.00–0.11 | 0.70–1.00 | 0.90–1.00 |
| Crossings without gaps > 30 min | 10 | 18 | 0.49–0.53 | 0.04–0.06 | 0.09–0.12 |
| Crossings without gaps > 30 min | 15 or 20 | 16 | 0.38–0.44 | 0.11–0.18 | 0.20–0.29 |

The correlations are positive throughout; none is significant at the 5% level after the adjustment for arc shares.

**Reference arcs (`S4_3`).**

| Reference arcs | Episode | OD pairs with a reference | Event flights with a profile | Comparable flights |
|---|---|---|---|---|
| All crossings | 2025 | 1,062 of 1,066 | 7,779 of 8,032 | 4,116 |
| All crossings | 2026 | 425 of 432 | 2,666 of 2,785 | 1,624 |
| Without gaps > 30 min | 2025 | 861 of 1,066 | 6,653 of 8,032 | 3,560 |
| Without gaps > 30 min | 2026 | 351 of 432 | 2,332 of 2,785 | 1,470 |

## S5. Sensitivity analyses

`S5_1_robustness_variants.csv` repeats the headline quantities of the profile and structural-floor analyses under the variants below. Each variant recomputes the clearances and re-runs every analysis step.

| Variant | Description |
|---|---|
| `main` | Main specification |
| `trim_50`, `trim_100` | Clearance recomputed after excluding trajectory points within 27.0 or 54.0 NM (50 or 100 km) of each flight's origin and destination airports |
| `omae_boundary` (2026) | OMAE removed from $S_0$; flights to and from UAE airports still excluded |
| `omae_open` (2026) | OMAE removed from $S_0$; flights to and from UAE airports admitted |

For each variant and affected episode the table reports the exposed sample, the monotonicity of the profile, $d_{1/2}$ with its interval and its range over the seven exposure definitions, and the median clearance, floor, and residual standoff with its interval, together with the residual at $q = 0$ and $q = 0.10$.

| Variant | Episode | $d_{1/2}$ (NM) | Median residual standoff (NM) |
|---|---|---|---|
| main | 2025 | 49.0 [46.9, 53.9] | 9.8 [9.3, 9.9] |
| trim_50 | 2025 | 49.4 [47.0, 52.9] | 9.3 [9.1, 9.9] |
| trim_100 | 2025 | 51.9 [47.9, 58.5] | 7.5 [7.2, 9.6] |
| main | 2026 | 40.5 [38.3, 54.9] | 9.9 [6.8, 13.1] |
| trim_50 | 2026 | 43.2 [38.7, 52.9] | 10.0 [6.9, 13.7] |
| trim_100 | 2026 | 45.5 [38.7, 54.1] | 4.8 [1.5, 6.8] |
| omae_boundary | 2026 | 67.2 [62.5, 71.9] | 30.8 [20.2, 40.6] |
| omae_open | 2026 | 57.2 [50.3, 65.6] | 14.5 [13.0, 20.3] |

The closure-exposed profile is monotone and the median residual standoff is positive in every variant. Under the OMAE codings, flights that passed around the United Arab Emirates are measured against more distant parts of the closure, which raises $d_{1/2}$ and the residual standoff.

## Code

[`code/`](code/) contains the code that computes every table in S2–S5 and Figs. 5–8 from flight trajectories, flight metadata, and FIR polygons. It reads the closure calendar from `tables/S1_1_closure_calendar.csv`. [`code/README.md`](code/README.md) describes the required inputs, the commands for the main specification and the sensitivity variants, and the random-number settings.

## Data availability

The trajectories were obtained from Flightradar24 under a license that does not permit redistribution of flight-level data. This repository therefore contains aggregate results and code only.

The FIR boundaries used in the paper come from public sources: the eurocontrol-atlas FIR/UIR data derived from the EUROCONTROL EAD, and the VATSpy Data Project for OTDF. They are not redistributed here.

## Citation

Citation details will be added after the review process.
