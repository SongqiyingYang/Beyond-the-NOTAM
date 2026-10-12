# Reproduction code

This folder contains the code that produces every result table in `../tables/` (S2 to S5) and Figs. 5–8 of the paper from flight trajectories. It implements the methods of Sections II and III. No data are included: the trajectories, flight metadata, and FIR polygons must be supplied by the user (see [Inputs](#inputs)).

Run on the trajectory data used in the paper, the code reproduces the published results. Point estimates, sample sizes, $F$-tests, and rank correlations are identical. Confidence intervals are identical as well, with one exception: the intervals of the local slopes, explained in [Random numbers](#random-numbers).

## Files

| File | Paper | Content |
|---|---|---|
| `common.py` | Sec. II | Episodes, region of interest, thresholds, projection, closure calendar, analysis sample, closure exposure |
| `step1_clearance.py` | Eq. (2) | Minimum clearance of every flight from the closure set; airport trimming and alternative closure sets |
| `step2_profile.py` | Sec. III-A | Revealed avoidance profile, benchmark fits, local slopes, half-decay distance, exposure definitions |
| `step3_floor.py` | Sec. III-B | Boundary anchors, cell and FIR floors, residual standoff, sensitivity analyses |
| `step4_operators.py` | Sec. III-C | Market–operator fixed effects, $F$-test, cross-event rank correlations, closure-exposed operator panel |
| `run_all.py` | — | Runs steps 2–4 for both episodes and writes tables S2–S4 |
| `collect_variants.py` | Sec. III-D | Combines several runs into Table S5.1 |
| `make_figures.py` | Figs. 5–8 | Draws the figures from the output of `run_all.py` |

All computations use kilometers, as in the paper's analysis. Tables and figures are reported in nautical miles by default (`--unit nm`; 1 NM = 1.852 km). `--unit km` reports kilometers instead.

## Inputs

1. **Trajectories.** Parquet files with one row per trajectory point on a 60-s grid. Columns: `fr24_id`, `timestamp` (UTC), `lat`, `lon`. The paper uses Flightradar24 ADS-B/MLAT data, which cannot be redistributed.

2. **Flight metadata.** Parquet or CSV with one row per flight. Columns:

   | Column | Content |
   |---|---|
   | `fr24_id` | Flight identifier, matching the trajectories |
   | `operator` | ICAO airline designator of the operating carrier |
   | `orig_icao`, `dest_icao` | ICAO codes of the origin and destination airports |
   | `category` | `Passenger` for passenger operations |
   | `flight` | Flight number; used only for the two service-level exposure definitions |

3. **FIR polygons.** Any file readable by geopandas, with an `IDENT` column holding the ICAO FIR designator. The paper uses the FIR/UIR boundaries of the eurocontrol-atlas project (derived from the EUROCONTROL EAD), with OTDF taken from the VATSpy Data Project. Both are public.

4. **Closure calendar.** Read from `../tables/S1_1_closure_calendar.csv`. The column `code` holds the coding: `A` (closed) and `B` (operationally closed) enter the closure set.

5. **Airport coordinates.** Needed only for the airport-proximity variants. They are taken from the `airportsdata` package and filled, through `--airports`, from OurAirports' public `airports.csv` (or any CSV with columns `icao`, `lat`, `lon`). The published variants use both sources; with `airportsdata` alone, a few airports have no position and their markets drop out of the market-dropping variants.

## Running the analysis

```bash
pip install -r requirements.txt

# Step 1, once per episode (about 10 min for the 2026 episode on a laptop)
python step1_clearance.py --episode 2025 --tracks TRACKS_2025/ --flights flights_2025.parquet \
       --firs firs.geojson --out work/flights_2025.parquet
python step1_clearance.py --episode 2026 --tracks-baseline TRACKS_2026/baseline/ \
       --tracks-event TRACKS_2026/event/ --flights flights_2026.parquet \
       --firs firs.geojson --out work/flights_2026.parquet

# Steps 2-4 and the figures (about 1 min)
python run_all.py --flights-2025 work/flights_2025.parquet --flights-2026 work/flights_2026.parquet \
       --firs firs.geojson --out results/main
python make_figures.py --run results/main --out figures/
```

`run_all.py` writes three kinds of output:

- `results/main/tables/`: the tables S2.1–S4.3 in the format of `../tables/`;
- `results/main/summary.json`: the headline quantities used by `collect_variants.py`;
- `results/main/private/`: flight-level residuals and cell-level validation data, used by the figures.

The files in `private/` are flight-level results derived from your input trajectories. Treat them like the input data: keep them out of any public repository if your data licence forbids redistribution of flight-level data. The `.gitignore` of this repository already excludes them.

## Sensitivity variants (Table S5.1)

| Variant | Step 1 options | `run_all.py` options |
|---|---|---|
| `trim_50`, `trim_100` | `--trim-km 50` (or `100`) `--airports airports.csv`, both episodes | — |
| `markets_dropped_50`, `markets_dropped_100` | as main | `--drop-markets-km 50` (or `100`) `--airports airports.csv` |
| `omae_boundary` | 2026: `--drop-fir OMAE` | `--geometry-drop 2026:OMAE` |
| `omae_open` | 2026: `--drop-fir OMAE` | `--geometry-drop 2026:OMAE --endpoint-drop 2026:OMAE` |
| `ostt_boundary` | 2025: `--drop-fir OSTT` | `--geometry-drop 2025:OSTT` |
| `ostt_open` | 2025: `--drop-fir OSTT` | `--geometry-drop 2025:OSTT --endpoint-drop 2025:OSTT` |
| `*_trim_50` | the coding variant combined with `--trim-km 50` | as the coding variant |

The two kinds of closure-coding option act at different stages:

- `--geometry-drop` removes the FIR from the closure boundary, the boundary from which clearance, anchors, and floors are measured.
- `--endpoint-drop` also readmits flights to and from that FIR's airports into the analysis sample.

The episode prefix matters because OSTT belongs to the closure set in both episodes but is coded as operationally closed only in 2025. Then combine the runs:

```bash
python collect_variants.py --out ../tables/S5_1_robustness_variants.csv \
       --variant main=results/main --variant trim_50=results/trim_50 \
       --variant omae_boundary=results/omae_boundary@2026 ...
```

`NAME=DIR@EPISODE` keeps only the episode that the closure-coding variant changes.

## Random numbers

All bootstraps resample OD markets with replacement and use NumPy's default generator.

| Quantity | Replications | Seed and drawing order |
|---|---|---|
| Profile | 2,000 | 20260825, separately for each episode; unexposed markets are drawn first, then closure-exposed markets |
| Half-decay distance, for all seven exposure definitions | 1,500 | 20260826, one stream: 2026 first, then 2025, in the order of `DEFINITIONS` |
| Median residual standoff | 1,000 | 20260825, separately for each episode |

The intervals of the local slopes use the profile bootstrap. In the archived run they came from a separate bootstrap of 1,200 replications, so the archived slopes differ from those of this code by at most 0.01 on closure-exposed markets and 0.10 on unexposed markets, and their interval bounds by at most 0.15. The paper's main text reports no slope values.

Rank correlations of $\widehat\delta_a$ are computed on the effects rounded to 0.1 km, the precision at which they are reported.
