#!/usr/bin/env python3
"""Read EQdyna.2Dcycle source-time-function output (stf.bin<icstart>).

The binary layout is defined in src/stf_output.f90; this module implements it.
One file per run segment holds a header with the static fault-node table and
one record per earthquake: for every fault node that slipped, the full dynamic
time series of

    slip_rate      |slip rate|                       m/s
    slip_rate_t    slip rate along the fault tangent  m/s  (signed)
    slip           |slip|                            m
    slip_t         slip along the fault tangent       m    (signed)
    shear_stress                                      Pa
    normal_stress  (negative = compression)           Pa
    friction       friction coefficient               -

plus the event id (= catalog.csv eqId), its time in the earthquake sequence,
the preceding interseismic interval, nucleation node, and per-node rupture time
and peak slip rate.

Library use:
    from stf_read import STFCase
    case = STFCase("work/my_case")            # finds stf.bin* in case dir / aRawSimuData
    case.summary()                             # one row per event
    ev = case.event(38)                        # dict of numpy arrays
    t, mdot = case.moment_rate(38)             # moment-rate function, N m / s

Command line:
    stf_read.py <case_dir|stf.bin*> --list
    stf_read.py <case_dir|stf.bin*> --event 38 --netcdf ev38.nc
    stf_read.py <case_dir|stf.bin*> --event 38 --npz ev38.npz

Moment rate uses the same convention as plotRuptureDynamics: mu = rou*vs^2 and
a 22 km seismogenic width, integrated over the tributary length of each node.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import struct
import sys

import numpy as np

MAGIC = b"EQDYNSTF"
MU_DEFAULT = 2670.0 * 3464.0 ** 2      # rou * vs^2, as plotRuptureDynamics
WIDTH_DEFAULT = 22.0e3                  # seismogenic width, m


class STFFile:
    """One stf.bin<icstart> file: header, node table and an event index."""

    def __init__(self, path: str):
        self.path = path
        with open(path, "rb") as f:
            self._read_header(f)
            self._index(f)

    # -- header -------------------------------------------------------------
    def _read_header(self, f):
        if f.read(8) != MAGIC:
            raise ValueError(f"{self.path}: not an EQdyna STF file")
        (self.version, self.icstart, self.ntotft, self.nnode,
         self.nvar) = struct.unpack("<5i", f.read(20))
        if self.version != 1:
            raise ValueError(f"{self.path}: unsupported STF version {self.version}")
        self.dt_solver, = struct.unpack("<d", f.read(8))
        self.stf_every, = struct.unpack("<i", f.read(4))
        self.vmin, self.t_seq_start_yr = struct.unpack("<2d", f.read(16))
        self.origin_known = bool(struct.unpack("<i", f.read(4))[0])
        raw = f.read(32 * self.nvar)
        names = [raw[16 * i:16 * i + 16].decode().strip() for i in range(2 * self.nvar)]
        self.varnames, self.varunits = names[:self.nvar], names[self.nvar:]
        rec = np.dtype([("fault", "<i4"), ("local", "<i4"), ("x", "<f8"), ("y", "<f8"),
                        ("tx", "<f8"), ("ty", "<f8"), ("length", "<f8")])
        self.nodes = np.frombuffer(f.read(rec.itemsize * self.nnode), dtype=rec)
        self._data_start = f.tell()

    # -- event index --------------------------------------------------------
    _EVHDR = struct.Struct("<i d d i d d d d d i i")

    def _index(self, f):
        self.index = {}
        f.seek(self._data_start)
        while True:
            tag = f.read(4)
            if not tag:
                break
            if tag != b"EVNT":
                raise ValueError(f"{self.path}: corrupt record at byte {f.tell() - 4}")
            nbytes, = struct.unpack("<q", f.read(8))
            start = f.tell()
            if start + nbytes > os.path.getsize(self.path):
                break                           # partially written last record
            h = self._EVHDR.unpack(f.read(self._EVHDR.size))
            self.index[h[0]] = dict(
                eqid=h[0], time_yr=h[1], interval_yr=h[2], nuc_node=h[3],
                nuc_x=h[4], nuc_y=h[5], t0=h[6], dt=h[7], t_end=h[8],
                nt=h[9], nout=h[10], offset=start)
            f.seek(start + nbytes)

    def event(self, eqid: int) -> dict:
        m = self.index[eqid]
        nt, nout, nv = m["nt"], m["nout"], self.nvar
        with open(self.path, "rb") as f:
            f.seek(m["offset"] + self._EVHDR.size)
            node = np.frombuffer(f.read(4 * nout), "<i4") - 1       # 0-based
            rtime = np.frombuffer(f.read(8 * nout), "<f8")
            vpeak = np.frombuffer(f.read(8 * nout), "<f8")
            data = np.frombuffer(f.read(4 * nt * nout * nv), "<f4").reshape(nv, nout, nt)
            trailer, = struct.unpack("<i", f.read(4))
        if trailer != eqid:
            raise ValueError(f"{self.path}: event {eqid} trailer mismatch ({trailer})")
        nd = self.nodes[node] if nout else self.nodes[:0]
        ev = dict(m)
        ev.update(
            time=m["t0"] + m["dt"] * np.arange(nt),
            node=node, fault=nd["fault"].copy(), x=nd["x"].copy(), y=nd["y"].copy(),
            tx=nd["tx"].copy(), ty=nd["ty"].copy(), length=nd["length"].copy(),
            rupture_time=rtime.copy(), peak_slip_rate=vpeak.copy(),
            origin_known=self.origin_known)
        for i, name in enumerate(self.varnames):
            ev[name] = data[i]                  # shape (nout, nt)
        return ev


class STFCase:
    """All STF segments of a case (stf.bin1, stf.bin1256, ...), merged by eqid."""

    def __init__(self, path: str, name: str | None = None):
        # case name for labelling figures: the case directory's name
        base = path if os.path.isdir(path) else os.path.dirname(os.path.abspath(path))
        if os.path.basename(base) == "aRawSimuData":
            base = os.path.dirname(base)
        self.name = name or os.path.basename(os.path.abspath(base))
        if os.path.isdir(path):
            files = []
            for d in (path, os.path.join(path, "aRawSimuData")):
                files += glob.glob(os.path.join(d, "stf.bin*"))
        else:
            files = [path]
        if not files:
            raise FileNotFoundError(f"no stf.bin* under {path}")
        key = lambda p: int(re.search(r"stf\.bin(\d+)$", p).group(1)) if re.search(r"stf\.bin(\d+)$", p) else 0
        self.files = [STFFile(p) for p in sorted(set(files), key=key)]
        self.owner = {}
        for sf in self.files:                   # later segments win on duplicates
            for eqid in sf.index:
                self.owner[eqid] = sf
        self.nodes = self.files[0].nodes
        self.varnames = self.files[0].varnames
        self.varunits = self.files[0].varunits

    @property
    def eqids(self):
        return sorted(self.owner)

    def meta(self, eqid):
        return self.owner[eqid].index[eqid]

    def event(self, eqid: int) -> dict:
        return self.owner[eqid].event(eqid)

    def moment_rate(self, eqid, mu=MU_DEFAULT, width=WIDTH_DEFAULT):
        """Moment-rate function M0dot(t) = mu * W * sum_i |slip rate|_i * length_i."""
        ev = self.event(eqid)
        mdot = mu * width * (ev["slip_rate"] * ev["length"][:, None]).sum(axis=0)
        return ev["time"], mdot

    def summary(self, mu=MU_DEFAULT, width=WIDTH_DEFAULT):
        rows = []
        for eqid in self.eqids:
            ev = self.event(eqid)
            m0 = mu * width * float((ev["slip"][:, -1] * ev["length"]).sum()) if ev["nout"] else 0.0
            # same magnitude formula as plotRuptureDynamics (M0 in dyne-cm)
            mw = 2.0 / 3.0 * np.log10(m0 * 1.0e7) - 10.7 if m0 > 0 else float("nan")
            rows.append((eqid, ev["time_yr"], ev["interval_yr"], ev["nout"], ev["nt"],
                         ev["t_end"], m0, mw))
        return rows


def to_netcdf(ev: dict, varnames, varunits, out: str):
    import xarray as xr
    coords = dict(time=("time", ev["time"]), node=("node", ev["node"] + 1),
                  x=("node", ev["x"]), y=("node", ev["y"]), fault=("node", ev["fault"]),
                  tx=("node", ev["tx"]), ty=("node", ev["ty"]), length=("node", ev["length"]))
    data = {n: (("node", "time"), ev[n], {"units": u}) for n, u in zip(varnames, varunits)}
    data["rupture_time"] = ("node", ev["rupture_time"], {"units": "s", "note": "1000 = never ruptured"})
    data["peak_slip_rate"] = ("node", ev["peak_slip_rate"], {"units": "m/s"})
    attrs = {k: ev[k] for k in ("eqid", "time_yr", "interval_yr", "nuc_node", "nuc_x", "nuc_y",
                                "t0", "dt", "t_end")}
    attrs["origin_known"] = int(ev["origin_known"])
    attrs["source"] = "EQdyna.2Dcycle stf.bin (src/stf_output.f90)"
    xr.Dataset(data, coords=coords, attrs=attrs).to_netcdf(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path", help="case directory or an stf.bin* file")
    ap.add_argument("--list", action="store_true", help="one line per event")
    ap.add_argument("--event", type=int, nargs="+", help="event id(s) to export")
    ap.add_argument("--netcdf", help="write the selected event to NetCDF (one event)")
    ap.add_argument("--npz", help="write the selected event to .npz (one event)")
    a = ap.parse_args()

    case = STFCase(a.path)
    if a.list or not a.event:
        print(f"{len(case.eqids)} events in {len(case.files)} file(s); "
              f"variables: {', '.join(case.varnames)}")
        print(f"{'eqid':>6} {'time_yr':>11} {'interval':>9} {'nodes':>6} {'nt':>6} "
              f"{'t_end_s':>8} {'M0_Nm':>10} {'Mw':>5}")
        for r in case.summary():
            print(f"{r[0]:6d} {r[1]:11.1f} {r[2]:9.1f} {r[3]:6d} {r[4]:6d} "
                  f"{r[5]:8.2f} {r[6]:10.3e} {r[7]:5.2f}")
    if a.event:
        if (a.netcdf or a.npz) and len(a.event) != 1:
            sys.exit("--netcdf/--npz take exactly one --event")
        ev = case.event(a.event[0])
        if a.netcdf:
            to_netcdf(ev, case.varnames, case.varunits, a.netcdf)
            print(f"wrote {a.netcdf}")
        if a.npz:
            np.savez_compressed(a.npz, **{k: v for k, v in ev.items() if not isinstance(v, dict)})
            print(f"wrote {a.npz}")


if __name__ == "__main__":
    main()
