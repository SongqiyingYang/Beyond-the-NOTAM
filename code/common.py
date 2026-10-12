"""Shared settings and helpers for the reproduction code.

All computations are carried out in kilometers, as in the paper's analysis; tables and
figures are reported in nautical miles (1 NM = 1.852 km) unless --unit km is given.
Every threshold below is the value used in the paper (in km).
"""
from pathlib import Path

import numpy as np
import pandas as pd
import shapely

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
CALENDAR = REPO / "tables" / "S1_1_closure_calendar.csv"

# ----------------------------------------------------------------------------- episodes
EPISODES = {
    "2025": dict(baseline_days=[f"2025-06-{d:02d}" for d in range(6, 12)],
                 eval_days=[f"2025-06-{d:02d}" for d in range(13, 18)]),
    "2026": dict(baseline_days=[f"2026-02-{d:02d}" for d in range(21, 28)],
                 eval_days=[f"2026-03-{d:02d}" for d in range(2, 7)]),
}
ROI = dict(lon_min=25.0, lon_max=75.0, lat_min=15.0, lat_max=45.0)
DAY_SLICE_HOUR = 6            # a "day" runs from 06:00 to 06:00 UTC
LAT0 = 30.0                   # latitude of the equirectangular projection (centre of the ROI)
KM_DEG = 111.195              # km per degree of latitude

# ----------------------------------------------------------------------------- Section III-A
DGRID = np.array([10, 15, 20, 25, 30, 40, 50, 60, 75, 90, 100, 120, 150, 175,
                  200, 225, 250, 275, 300, 350, 400, 450, 500, 600, 750, 900], float)
FIT_NODES = np.array([20, 25, 30, 40, 50, 60, 75, 90, 100, 120, 150, 175, 200], float)
SIGMA_GRID = np.array([10, 15, 20, 25, 30, 40, 50, 65, 80, 100, 130, 160, 200,
                       250, 320, 400, 500, 650, 800], float)
XI_BAR_GRID = np.arange(20, 620, 20, dtype=float)
LOG_FLOOR = 1e-6              # predicted values are floored before taking logarithms
SLOPE_WINDOWS = [(20, 60), (60, 150), (150, 300)]
D_HALF_ANCHOR = 50.0
EXPOSURE_TAUS = (0.25, 0.50, 0.75, 1.0)
B_PROFILE, SEED_PROFILE = 2000, 20260825
B_DHALF, SEED_DHALF = 1500, 20260826

# ----------------------------------------------------------------------------- Section III-B
CELL_DEG = 1.5
CELL_MIN_FLIGHTS = 40
CELL_NEAR_KM = 50.0
CELL_MIN_NEAR = 20
FIR_MIN_FLIGHTS = 100
FLOOR_Q = 0.05
VALIDATION_MIN_EVENT = 20
B_RESIDUAL, SEED_RESIDUAL = 1000, 20260825

# ----------------------------------------------------------------------------- Section III-C
MARKET_MAX_MEDIAN_KM = 400.0
MARKET_MIN_BASELINE = 3
OPERATOR_MIN_FLIGHTS = 25
PANEL_MIN_EXPOSED = 20

# ICAO location-indicator prefix -> FIR, used to place origin and destination airports
PREFIX_FIR = {
    "OE": "OEJD", "OM": "OMAE", "OB": "OBBB", "OT": "OTDF", "OK": "OKAC",
    "OI": "OIIX", "OR": "ORBB", "OO": "OOMM", "OJ": "OJAC", "OS": "OSTT",
    "OY": "OYSC", "LL": "LLLL", "LC": "LCCC", "LG": "LGGG", "LT": "LTAA",
    "HE": "HECC", "HS": "HSSS", "HH": "HHAA", "HL": "HLLL", "OP": "OPKR",
    "OA": "OAKX", "UB": "UBBA", "UG": "UGGG", "UT": "UTAK",
}

NM_PER_KM = 1.0 / 1.852


# ----------------------------------------------------------------------------- units
class Unit:
    """Conversion of reported distances; the analysis itself always runs in km."""

    def __init__(self, name="nm"):
        assert name in ("km", "nm"), name
        self.name = name
        self.f = 1.0 if name == "km" else NM_PER_KM
        self.tag = name                       # column suffix
        self.label = "km" if name == "km" else "NM"

    def __call__(self, x, nd=1):
        """Convert a distance (scalar, array or Series) from km and round."""
        if np.ndim(x):
            return np.round(np.asarray(x, float) * self.f, nd) + 0.0      # + 0.0 turns -0.0 into 0.0
        return round(float(x) * self.f, nd) + 0.0

    def thr(self, km):
        """Label of a threshold set in km: '50 km' or '27.0 NM'."""
        return f"{km:g} km" if self.name == "km" else f"{km * self.f:.1f} NM"


# ----------------------------------------------------------------------------- geometry
def project(lon, lat):
    """Equirectangular projection to km (x = lon cos(lat0), y = lat)."""
    return np.c_[np.asarray(lon, float) * np.cos(np.radians(LAT0)) * KM_DEG,
                 np.asarray(lat, float) * KM_DEG]


def project_geom(g):
    return shapely.transform(g, lambda c: project(c[:, 0], c[:, 1]))


def unproject(xy):
    xy = np.asarray(xy, float)
    return xy[:, 0] / (np.cos(np.radians(LAT0)) * KM_DEG), xy[:, 1] / KM_DEG


def load_firs(path):
    """FIR polygons (GeoJSON or any format read by geopandas) with an IDENT column."""
    import geopandas as gpd
    g = gpd.read_file(path)
    if "IDENT" not in g.columns:
        raise ValueError(f"{path}: the FIR file needs an IDENT column with the ICAO FIR designator")
    return g.set_index("IDENT").geometry


def airport_coordinates(extra=None):
    """ICAO code -> (lat, lon) from the airportsdata package, filled from an optional CSV:
    either OurAirports' airports.csv (icao_code, latitude_deg, longitude_deg) or a table
    with columns icao, lat, lon."""
    import airportsdata
    ap = {k: (v["lat"], v["lon"]) for k, v in airportsdata.load("ICAO").items()}
    if extra:
        x = pd.read_csv(extra, low_memory=False)
        cols = (("icao_code", "latitude_deg", "longitude_deg") if "icao_code" in x.columns
                else ("icao", "lat", "lon"))
        for code, lat, lon in x[list(cols)].dropna().itertuples(index=False):
            ap.setdefault(str(code), (float(lat), float(lon)))
    return ap


def closure_geometry(firs, closure_set):
    """Union of the closed FIRs in projected km coordinates."""
    return project_geom(firs.loc[sorted(closure_set)].union_all())


def airport_fir(code):
    if not isinstance(code, str) or len(code) < 2:
        return None
    return PREFIX_FIR.get(code[:2])


# ----------------------------------------------------------------------------- closure calendar
def load_calendar(path=CALENDAR):
    cal = pd.read_csv(path, dtype=str).fillna("")
    cal = cal[cal.code.isin(["A", "B"])].copy()
    cal["closed_from_utc"] = pd.to_datetime(cal.closed_from_utc)
    cal["closed_to_utc"] = pd.to_datetime(cal.closed_to_utc)
    return cal


def closed_set(episode, day, calendar=None):
    """FIRs closed (coding A) or operationally closed (coding B) for at least half of the
    24-h day slice that starts at 06:00 UTC on `day`."""
    cal = load_calendar() if calendar is None else calendar
    cal = cal[cal.episode == episode]
    d0 = pd.Timestamp(day) + pd.Timedelta(hours=DAY_SLICE_HOUR)
    d1 = d0 + pd.Timedelta(hours=24)
    out = set()
    for fir, sub in cal.groupby("fir"):
        iv = sorted((max(a, d0), min(b, d1)) for a, b in zip(sub.closed_from_utc, sub.closed_to_utc)
                    if a < d1 and b > d0)
        if not iv:
            continue
        merged = [list(iv[0])]
        for a, b in iv[1:]:
            if a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        if sum((b - a).total_seconds() for a, b in merged) / 86400 >= 0.5:
            out.add(fir)
    return out


def reference_closure(episode, calendar=None):
    """Evaluation day with the most closed FIRs and its closure set. Baseline clearances are
    measured against this set (the airspace that subsequently became unavailable)."""
    days = EPISODES[episode]["eval_days"]
    sets = {d: closed_set(episode, d, calendar) for d in days}
    ref = max(days, key=lambda d: len(sets[d]))
    return ref, sorted(sets[ref])


# ----------------------------------------------------------------------------- samples
def analysis_sample(F, endpoint_exclusion):
    """Passenger flights with an identified operator whose origin and destination airports
    are not inside a FIR of the closure set (Section II-E)."""
    F = F[(F.category == "Passenger") & F.operator.notna()].copy()
    of = F.orig_icao.astype(str).map(airport_fir)
    df = F.dest_icao.astype(str).map(airport_fir)
    return F[~(of.isin(endpoint_exclusion) | df.isin(endpoint_exclusion))].copy()


def split_exposure(F):
    """Markets observed in both periods, split into closure-exposed (r_o > 0: at least one
    baseline flight entered the subsequently closed airspace) and unexposed (r_o = 0)."""
    touched = (F[F.period == "baseline"].assign(z=lambda x: x.clearance_km <= 0.0)
               .groupby("od").z.any())
    both = F.groupby("od").period.nunique() == 2
    keep = both[both].index
    exp = [o for o in keep if touched.get(o, False)]
    une = [o for o in keep if not touched.get(o, False)]
    return {"exposed": F[F.od.isin(exp)].copy(), "unexposed": F[F.od.isin(une)].copy(),
            "both": F[F.od.isin(keep)].copy()}


def penetration_rate(F):
    """r_o = share of a market's baseline flights with zero clearance (Eq. 3)."""
    b = F[F.period == "baseline"]
    return b.assign(z=b.clearance_km <= 0).groupby("od").z.mean()
