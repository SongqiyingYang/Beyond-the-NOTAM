# Reproduction code

This folder contains the code that produces every result table in `../tables/` (S2 to S5) and Figs. 5–8 of the paper from flight trajectories. It implements the methods of Sections II and III. No data are included: the trajectories, flight metadata, and FIR polygons must be supplied by the user (see [Inputs](#inputs)).

## Files

| File | Paper | Content |
|---|---|---|
| `common.py` | Sec. II | Episodes, region of interest, thresholds, projection, closure calendar, analysis sample, closure exposure |
| `step1_clearance.py` | Eq. (1) | Minimum clearance of every flight from the closure set; airport trimming and the alternative OMAE coding |
| `step2_profile.py` | Sec. III-A | Revealed avoidance profile, two-parameter benchmark fits, half-decay distance, exposure definitions |
| `step3_floor.py` | Sec. III-B | Boundary anchors, cell and FIR floors, residual standoff, sensitivity analyses |
| `step4_operators.py` | Sec. III-C | Detour arcs, boundary-distance profiles, operator-relative standoff, cross-event rank correlations with permutation tests |
| `run_all.py` | — | Runs steps 2–3 for both episodes and writes tables S2–S3 |
| `collect_variants.py` | Sec. III-D | Combines several runs into Table S5.1 |
| `make_figures.py` | Figs. 5–8 | Draws the figures from the output of `run_all.py` and `step4_operators.py` |

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

5. **Airport coordinates.** Needed only for the airport-trimming variants. They are taken from the `airportsdata` package and can be filled, through `--airports`, from OurAirports' public `airports.csv` (or any CSV with columns `icao`, `lat`, `lon`).

## Running the analysis

```bash
pip install -r requirements.txt

# Step 1, once per episode (about 10 min for the 2026 episode on a laptop)
python step1_clearance.py --episode 2025 --tracks TRACKS_2025/ --flights flights_2025.parquet \
       --firs firs.geojson --out work/flights_2025.parquet
python step1_clearance.py --episode 2026 --tracks-baseline TRACKS_2026/baseline/ \
       --tracks-event TRACKS_2026/event/ --flights flights_2026.parquet \
       --firs firs.geojson --out work/flights_2026.parquet

# Steps 2-3 (about 2 min)
python run_all.py --flights-2025 work/flights_2025.parquet --flights-2026 work/flights_2026.parquet \
       --firs firs.geojson --out results/main

# Step 4 reads the trajectories again (about 8 min, mostly the permutation tests)
python step4_operators.py --flights-2025 work/flights_2025.parquet --flights-2026 work/flights_2026.parquet \
       --tracks-2025 TRACKS_2025/ --tracks-2026-baseline TRACKS_2026/baseline/ \
       --tracks-2026-event TRACKS_2026/event/ --firs firs.geojson --out results/main

python make_figures.py --run results/main --out figures/
```

`run_all.py` writes three kinds of output:

- `results/main/tables/`: the tables S2.1–S3.6 in the format of `../tables/` (S4.1–S4.3 are added by `step4_operators.py`);
- `results/main/summary.json`: the headline quantities used by `collect_variants.py`;
- `results/main/private/`: flight-level residuals, cell-level validation data, and flight-level standoff $m_f$, used by the figures.

The files in `private/` are flight-level results derived from your input trajectories. Treat them like the input data: keep them out of any public repository if your data licence forbids redistribution of flight-level data.

## Sensitivity variants (Table S5.1)

| Variant | Step 1 options | `run_all.py` options |
|---|---|---|
| `trim_50`, `trim_100` | `--trim-km 50` (or `100`), both episodes | — |
| `omae_boundary` | 2026: `--drop-fir OMAE` | `--geometry-drop 2026:OMAE` |
| `omae_open` | 2026: `--drop-fir OMAE` | `--geometry-drop 2026:OMAE --endpoint-drop 2026:OMAE` |

The two closure-coding options act at different stages:

- `--geometry-drop` removes OMAE from the closure boundary, the boundary from which clearance, anchors, and floors are measured.
- `--endpoint-drop` also readmits flights to and from UAE airports into the analysis sample.

Then combine the runs:

```bash
python collect_variants.py --out ../tables/S5_1_robustness_variants.csv \
       --variant main=results/main --variant trim_50=results/trim_50 --variant trim_100=results/trim_100 \
       --variant omae_boundary=results/omae_boundary@2026 --variant omae_open=results/omae_open@2026
```

`NAME=DIR@EPISODE` keeps only the episode that the closure-coding variant changes.

## Random numbers

All bootstraps resample OD pairs with replacement and use NumPy's default generator.

| Quantity | Replications | Seed and drawing order |
|---|---|---|
| Profile | 2,000 | 20260825, separately for each episode; unexposed OD pairs are drawn first, then closure-exposed pairs |
| Half-decay distance, for all seven exposure definitions | 1,500 | 20260826, one stream: 2026 first, then 2025, in the order of `DEFINITIONS` |
| Median residual standoff | 1,000 | 20260825, separately for each episode |
| Operator permutation tests | 2,000 per reference-arc rule | 20261021, one stream: all baseline crossings first, then crossings without long gaps; within each permutation, 2025 before 2026 and groups in sorted order |

In each permutation, the same relabelled flights serve all arc shares and operator thresholds, so that the family-adjusted $p$-value compares like with like.
