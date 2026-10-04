"""stf_read.py against the binary layout documented in src/stf_output.f90.

The file here is assembled byte by byte from that documented layout, so a
reader that drifts from the writer's contract fails these tests.
"""
import struct

import numpy as np
import pytest

from stf_read import STFCase, STFFile, MU_DEFAULT, WIDTH_DEFAULT

NVAR = 7
NAMES = ["slip_rate", "slip_rate_t", "slip", "slip_t", "shear_stress", "normal_stress", "friction"]
UNITS = ["m/s", "m/s", "m", "m", "Pa", "Pa", "1"]
# node table: fault, local, x, y, tx, ty, length
NODES = [(1, 1, 0.0, 0.0, 1.0, 0.0, 100.0), (1, 2, 100.0, 0.0, 1.0, 0.0, 200.0),
         (2, 1, 0.0, 500.0, 1.0, 0.0, 150.0), (2, 2, 100.0, 500.0, 1.0, 0.0, 50.0)]


def header(icstart=1, origin=1):
    b = b"EQDYNSTF" + struct.pack("<5i", 1, icstart, 2, len(NODES), NVAR)
    b += struct.pack("<d i d d i", 0.01, 1, 1e-3, 0.0, origin)
    b += b"".join(n.ljust(16).encode() for n in NAMES) + b"".join(u.ljust(16).encode() for u in UNITS)
    for nd in NODES:
        b += struct.pack("<2i5d", *nd)
    return b


def event(eqid, time_yr, interval, nodes, data, rupt, vpeak, trailer=None):
    """nodes are 1-based; data has shape (nvar, nout, nt)."""
    nout, nt = len(nodes), (data.shape[2] if nodes else 0)
    nbytes = 72 + 20 * nout + 4 * NVAR * nout * nt + 4      # stf_output.f90's formula
    b = b"EVNT" + struct.pack("<q", nbytes)
    b += struct.pack("<i d d i d d d d d i i", eqid, time_yr, interval, nodes[0] if nodes else 1,
                     0.0, 0.0, 0.01, 0.01, 0.01 * nt, nt, nout)
    if nout:
        b += np.asarray(nodes, "<i4").tobytes() + np.asarray(rupt, "<f8").tobytes()
        b += np.asarray(vpeak, "<f8").tobytes() + np.asarray(data, "<f4").tobytes()
    b += struct.pack("<i", eqid if trailer is None else trailer)
    assert len(b) == 12 + nbytes
    return b


def sample_data(nout, nt, seed):
    rng = np.random.default_rng(seed)
    return rng.random((NVAR, nout, nt)).astype("<f4")


@pytest.fixture
def two_events(tmp_path):
    d1 = sample_data(2, 3, 1)
    p = tmp_path / "stf.bin1"
    p.write_bytes(header()
                  + event(1, 100.0, 100.0, [1, 3], d1, [0.5, 1.5], [1.2, 0.8])
                  + event(2, 130.0, 30.0, [], np.zeros((NVAR, 0, 0), "<f4"), [], []))
    return p, d1


def test_header_and_index(two_events):
    p, _ = two_events
    f = STFFile(str(p))
    assert (f.version, f.icstart, f.ntotft, f.nnode, f.nvar) == (1, 1, 2, 4, NVAR)
    assert f.varnames == NAMES and f.varunits == UNITS
    assert sorted(f.index) == [1, 2]
    assert f.index[1]["time_yr"] == 100.0 and f.index[2]["interval_yr"] == 30.0


def test_event_values_and_node_mapping(two_events):
    p, d1 = two_events
    ev = STFCase(str(p)).event(1)
    assert list(ev["node"]) == [0, 2]                       # 1-based on disk, 0-based here
    assert list(ev["fault"]) == [1, 2] and list(ev["length"]) == [100.0, 150.0]
    np.testing.assert_array_equal(ev["slip_rate"], d1[0])   # var 1, (nout, nt)
    np.testing.assert_array_equal(ev["friction"], d1[6])    # var 7
    np.testing.assert_allclose(ev["time"], [0.01, 0.02, 0.03])
    assert list(ev["rupture_time"]) == [0.5, 1.5]


def test_empty_event(two_events):
    p, _ = two_events
    ev = STFCase(str(p)).event(2)
    assert ev["nout"] == 0 and ev["slip_rate"].shape == (0, 0)


def test_moment_rate(two_events):
    p, d1 = two_events
    t, mdot = STFCase(str(p)).moment_rate(1)
    want = MU_DEFAULT * WIDTH_DEFAULT * (d1[0] * np.array([100.0, 150.0])[:, None]).sum(0)
    np.testing.assert_allclose(mdot, want, rtol=1e-6)


def test_partial_last_record_is_ignored(two_events, tmp_path):
    p, _ = two_events
    full = event(3, 150.0, 20.0, [2], sample_data(1, 4, 3), [0.1], [2.0])
    p.write_bytes(p.read_bytes() + full[:-10])               # a run still writing event 3
    assert sorted(STFFile(str(p)).index) == [1, 2]


def test_trailer_mismatch_is_detected(tmp_path):
    p = tmp_path / "stf.bin1"
    p.write_bytes(header() + event(1, 1.0, 1.0, [1], sample_data(1, 2, 4), [0.1], [1.0], trailer=99))
    with pytest.raises(ValueError, match="trailer"):
        STFCase(str(p)).event(1)


def test_segments_merge_later_wins(tmp_path):
    (tmp_path / "stf.bin1").write_bytes(header(1) + event(1, 10.0, 10.0, [1], sample_data(1, 2, 5), [0.1], [1.0]))
    (tmp_path / "stf.bin2").write_bytes(header(2) + event(2, 20.0, 10.0, [2], sample_data(1, 2, 6), [0.2], [1.0]))
    c = STFCase(str(tmp_path))
    assert c.eqids == [1, 2] and c.name == tmp_path.name
