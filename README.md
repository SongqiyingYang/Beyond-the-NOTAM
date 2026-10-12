# Beyond the NOTAM: Supplementary Repository

Supplementary material for the paper *Beyond the NOTAM: Revealed Airspace-Avoidance Behaviour During Major Airspace Closures* (under double-blind review; author information withheld).

The paper uses before-and-after commercial flight trajectories from two conflict-related airspace-closure episodes in the Middle East (June 2025 and March 2026). It measures the routing response beyond the published closure boundary through three tools: a revealed avoidance profile, a baseline-revealed structural floor, and market–operator fixed effects. This repository contains two parts:

- **`tables/`**: the complete results behind the main text. These cover the closure coding, the full profiles and benchmark fits, the structural-floor construction, the operator estimates, and the sensitivity analyses referred to in the paper.
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
| `tables/S2_1_avoidance_profile.csv` | $\widehat h(d)$ with 95% intervals, closure-exposed and unexposed markets, 5.4–486.0 NM | Eq. (9), Fig. 5(a) |
| `tables/S2_2_benchmark_fits.csv` | Hard-buffer and exponential benchmark fits, half-decay distance | Eqs. (12)–(13), (15), Fig. 5(b) |
| `tables/S2_3_local_slopes.csv` | Local log–log slope $\beta(d)$ over three distance windows | Eq. (14) |
| `tables/S2_4_exposure_definitions.csv` | Profile summaries under the seven exposure definitions | Eq. (16) |
| **S3. Revealed structural floor** | | |
| `tables/S3_1_boundary_cells.csv` | Every boundary-anchor cell: baseline counts, floor, eligibility, event-period lower tail | Eq. (17), Figs. 4 and 6 |
| `tables/S3_2_fir_fallback_floors.csv` | FIR-level fallback floors | Sec. III-B |
| `tables/S3_3_decomposition.csv` | Median clearance, floor, and residual standoff with 95% intervals | Eq. (18), Fig. 7 |
| `tables/S3_4_floor_quantile.csv` | Median residual standoff for $q \in \{0, 0.01, 0.02, 0.05, 0.10\}$ | Sec. IV-B |
| `tables/S3_5_eligibility_cell_size.csv` | Median residual standoff under alternative cell-eligibility rules and cell sizes | Sec. III-D |
| `tables/S3_6_residual_quantiles.csv` | Quantiles of the residual-standoff distribution | Fig. 7 |
| **S4. Operator-associated heterogeneity** | | |
| `tables/S4_1_operator_effects_2025.csv`, `..._2026.csv` | $\widehat\delta_a$ for every operator, with standard errors and 95% intervals | Eq. (19), Fig. 8 |
| `tables/S4_2_homogeneity_tests.csv` | Estimation samples and joint $F$-tests of equal operator effects | Sec. IV-C |
| `tables/S4_3_rank_correlations.csv` | Cross-event Spearman correlations on the full panel and on restricted panels | Eq. (20), Sec. III-C |
| **S5. Sensitivity analyses** | | |
| `tables/S5_1_robustness_variants.csv` | Headline quantities under the airport-proximity treatments and the alternative closure codings | Sec. III-D |
| **Code** | | |
| `code/` | Reproduction code; see [`code/README.md`](code/README.md) | Secs. II–III |

## Conventions

- **Episodes.** The 2025 episode compares the baseline window 6–11 June 2025 with the evaluation window 13–17 June 2025; its closure set is LLLL, OIIX, ORBB, OSTT. The 2026 episode compares 21–27 February 2026 with 2–6 March 2026; its closure set is LLLL, OBBB, OIIX, OKAC, OMAE, ORBB, OSTT, OTDF. Trajectories are on a 60-s grid.
- **Markets.** A directional origin–destination (OD) market is *closure-exposed* when at least one of its baseline flights entered the airspace that subsequently closed ($r_o > 0$, Eq. (4)). Otherwise it is *unexposed*. Only markets observed in both periods are used.
- **Distances.** All distances are horizontal distances from the boundary of the FIR-level closure set $S_0$, reported in nautical miles (NM; 1 NM = 1.852 km). Column names carry the unit (`_nm`).
  - The analysis was carried out in kilometers. Its thresholds were set as round values in km and appear here at their converted values: for example, 50 km is 27.0 NM, and the fitting window of 20–200 km is 10.8–108.0 NM.
  - Statistics without a distance unit do not depend on this conversion: profile values, error ratios, $F$, $\rho_S$, and $p$.
- **Intervals.** Intervals are 95% percentile intervals from nonparametric bootstraps clustered by OD market:

  | Quantity | Replications |
  |---|---|
  | Profile | $B = 2000$ |
  | Half-decay distance | $B = 1500$ |
  | Median residual standoff | $B = 1000$ |

  Operator effects use OD-clustered standard errors.
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

Such restrictions can move event-period traffic outward without moving the structural floor. Under a literal reading of the Jeddah restriction, 33 of the 2,785 closure-exposed 2026 event flights (1.2%) have their closest approach on the affected stretch of the boundary.

## S2. Revealed avoidance profile

**Profile (`S2_1`).** $\widehat h(d)$ of Eq. (9) is tabulated at every distance node from 5.4 to 486.0 NM, for closure-exposed and unexposed markets in both episodes. The columns `event_flights_within_d` and `baseline_flights_within_d` give the number of flights within distance $d$ in each period. On closure-exposed markets the profile decreases monotonically in both episodes. On unexposed markets it does not:

- 2025: $\widehat h(54.0) = 0.014 < \widehat h(108.0) = 0.032$;
- 2026: $\widehat h(27.0) = 0.083 < \widehat h(54.0) = 0.114$.

**Benchmarks (`S2_2`).** Both reference shapes are fitted by weighted least squares on the logarithmic scale, on the 13 distance thresholds between 10.8 and 108.0 NM.

- **Weights:** $(\widehat h/\mathrm{se})^2$ (delta method).
- **Zero predictions:** the hard buffer predicts zero beyond $\bar\xi$; such values are floored at $10^{-6}$ before taking logarithms.
- **Parameter grids:**

  | Parameter | Grid |
  |---|---|
  | $\bar\xi$ | 10.8–324.0 NM in steps of 10.8 NM (20–600 km in steps of 20 km) |
  | $\sigma$ | 5.4, 8.1, 10.8, 13.5, 16.2, 21.6, 27.0, 35.1, 43.2, 54.0, 70.2, 86.4, 108.0, 135.0, 172.8, 216.0, 270.0, 351.0, 432.0 NM |

- **Fit quality:** the unweighted root-mean-square error of $\log \widehat h$ over the same 13 thresholds.

| Episode | $\sigma$ (NM) | $\bar\xi$ (NM) | log-RMSE, exponential | log-RMSE, hard buffer | Ratio | $d_{1/2}$ (NM) |
|---|---|---|---|---|---|---|
| 2025 | 21.6 | 108.0 | 0.406 | 2.269 | 5.6 | 49.0 [46.9, 53.9] |
| 2026 | 27.0 | 108.0 | 0.174 | 2.128 | 12.3 | 40.5 [38.3, 54.9] |

**Local slopes (`S2_3`).** On closure-exposed markets the slope steepens with distance in both episodes:

| Episode | 10.8–32.4 NM | 32.4–81.0 NM | 81.0–162.0 NM |
|---|---|---|---|
| 2025 | −0.99 | −1.46 | −3.18 |
| 2026 | −0.69 | −1.48 | −2.40 |

On unexposed markets it is flatter, and in 2025 it does not steepen at all (−1.10, −1.04, −1.11).

**Exposure definitions (`S2_4`).** The table covers seven definitions:

- the main rule $r_o > 0$;
- the stricter thresholds $r_o \ge \tau$ of Eq. (16);
- two definitions at the level of continuing services, where a service is an operator and flight number observed in both periods.

The closure-exposed profile is monotone under all seven.

| Episode | Range of $d_{1/2}$ | Hard-buffer error relative to exponential |
|---|---|---|
| 2025 | 47.0–49.1 NM | 5.5–10.7 times |
| 2026 | 38.0–43.8 NM | 4.9–15.0 times |

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

## S4. Operator-associated heterogeneity

**Estimation sample (`S4_2`).** Model (19) is estimated on market–operator cells observed in both periods. The OD markets are those with a baseline median clearance below 216.0 NM (400 km) and at least three baseline flights; both closure-exposed and unexposed markets are included. Operators are included if they have at least 25 flights in each period on these markets.

| Episode | Operators | Flights | Market–operator cells | Markets | $F$ | $p$ | Span of $\widehat\delta_a$ |
|---|---|---|---|---|---|---|---|
| 2025 | 110 | 59,536 | 4,503 | 2,771 | 6.72 | $4.5 \times 10^{-80}$ | 71.2 NM |
| 2026 | 81 | 32,520 | 2,399 | 1,522 | 5.08 | $1.0 \times 10^{-38}$ | 52.0 NM |

**Operator effects (`S4_1`).** $\widehat\delta_a$ is the operator-associated before-and-after change in clearance relative to the operator's own baseline routing on the same markets. Seventy-six operators appear in both episodes.

**Rank correlations (`S4_3`).**

| Measure | Panel | $\rho_S$ | $p$ |
|---|---|---|---|
| $\widehat\delta_a$ | 76 operators observed in both episodes | 0.415 | $1.9 \times 10^{-4}$ |
| $\widehat\delta_a$ | 24 operators with ≥ 100 flights in each period of each episode | 0.519 | 0.009 |
| $\widehat\delta_a$ | 12 operators with ≥ 200 flights in each period of each episode | 0.685 | 0.014 |
| Raw event-period clearance | 23 operators with ≥ 20 closure-exposed flights in each period of each episode | 0.511 | 0.013 |
| Residual standoff $\widehat e_f$ | same 23 operators | 0.472 | 0.023 |
| $\widehat\delta_a$ | same 23 operators | 0.175 | 0.42 |

The 23-operator panel is restricted to closure-exposed flights, whereas $\widehat\delta_a$ is estimated on all markets below 216.0 NM. The two orderings therefore need not coincide on the same operators.

## S5. Sensitivity analyses

`S5_1_robustness_variants.csv` repeats the headline quantities of each analysis step under the variants below. Each variant re-runs every analysis step; the trimming and closure-coding variants first recompute the clearances.

| Variant | Description |
|---|---|
| `main` | Main specification |
| `trim_50`, `trim_100` | Clearance recomputed after excluding trajectory points within 27.0 or 54.0 NM (50 or 100 km) of each flight's origin and destination airports; all markets are kept |
| `markets_dropped_50`, `markets_dropped_100` | Markets served from an airport within 27.0 or 54.0 NM of the closure set dropped |
| `omae_boundary` (2026) | OMAE removed from $S_0$; flights to and from UAE airports still excluded |
| `omae_open` (2026) | OMAE removed from $S_0$; flights to and from UAE airports admitted |
| `ostt_boundary`, `ostt_open` (2025) | The same two codings for OSTT |
| `*_trim_50` | The coding variants combined with the 27.0 NM trajectory trimming |

The market-dropping variants are reported as a comparison only, because they remove a large share of the exposed traffic: Dubai alone carries 30% of the 2025 exposed traffic.

The table reports, for each variant and each affected episode:

- the exposed sample and the monotonicity of the profile;
- $d_{1/2}$ with its interval, and its range over the seven exposure definitions;
- the ratio of the hard-buffer error to the exponential error;
- the median clearance, the floor, and the residual standoff with its interval, together with the residual at $q = 0$ and $q = 0.10$;
- the cross-event $\rho_S$ of $\widehat\delta_a$.

Across all variants:

- **Profile.** The closure-exposed profile is monotone, and the hard buffer fits at least 3.8 times worse than the exponential decay.
- **Half-decay distance in 2025.** $d_{1/2}$ lies between 48.4 and 53.2 NM.
- **Half-decay distance in 2026.** Under the airport treatments and the OSTT codings, $d_{1/2}$ lies between 40.4 and 45.5 NM. Under the OMAE codings it rises to between 49.7 and 67.2 NM. When OMAE is removed from the boundary, flights that passed around the United Arab Emirates are measured against more distant parts of the closure.
- **Residual standoff.** The median residual standoff is positive throughout, from 4.8 to 30.8 NM. The lowest value, 2026 with 54.0 NM trimming, is driven by one boundary cell near Muscat whose baseline sample lies exactly at the eligibility threshold. With an eligibility threshold of 30 near-boundary flights, the value is 11.0 NM.
- **Operator ordering.** The cross-event $\rho_S$ of $\widehat\delta_a$ ranges from 0.337 to 0.460, with $p < 0.002$ in every variant.

## Code

[`code/`](code/) contains the code that computes every table in S2–S5 and Figs. 5–8 from flight trajectories, flight metadata, and FIR polygons. It reads the closure calendar from `tables/S1_1_closure_calendar.csv`.

Run on the data used in the paper, it reproduces all published numbers. The only exception is the local slopes and their intervals, whose bootstrap differs slightly from the archived run; the differences are documented in [`code/README.md`](code/README.md). That README also describes the required inputs, the commands for the main specification and for each sensitivity variant, and the random-number settings.

## Data availability

The trajectories were obtained from Flightradar24 under a license that does not permit redistribution of flight-level data. This repository therefore contains aggregate results and code only.

The FIR boundaries used in the paper come from public sources: the eurocontrol-atlas FIR/UIR data derived from the EUROCONTROL EAD, and the VATSpy Data Project for OTDF. They are not redistributed here.

## Citation

Citation details will be added after the review process.
