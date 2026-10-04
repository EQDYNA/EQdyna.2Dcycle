"""saf_model_to_geo.py: the inverse of the published SAF geometry chain."""
import numpy as np

import saf_model_to_geo as g


def test_matches_published_matlab_values():
    assert g.check() == 0


def test_model_utm_round_trip():
    rng = np.random.default_rng(0)
    x, y = rng.uniform(-250e3, 150e3, 200), rng.uniform(-30e3, 60e3, 200)
    x2, y2 = g.utm_to_model(*g.model_to_utm(x, y))
    np.testing.assert_allclose(x2, x, atol=1e-6)
    np.testing.assert_allclose(y2, y, atol=1e-6)


def test_lonlat_stays_in_zone_11_west_of_120():
    # the SAF's NW end: zone-11 easting ~200 km, longitude west of -120 deg
    lon, lat = g.utm_to_lonlat(201787.687, 3960338.182)
    assert -121.0 < float(lon) < -120.0 and 35.0 < float(lat) < 36.5
    e, n, _ = g.deg2utm(lat, lon, g.ZONE)
    assert abs(float(e) - 201787.687) < 1e-3 and abs(float(n) - 3960338.182) < 1e-3
