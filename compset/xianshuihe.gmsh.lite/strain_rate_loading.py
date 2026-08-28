#!/usr/bin/env python3
"""Max shear strain rate and its direction along the Xianshuihe faults.

Supplies the two loading columns of `nsmpGeoPhys.txt` that cannot be read
off a published slip-rate profile:

    ftLoadMaxShear   far-field maximum shear strain rate, s^-1
    ftLoadAngle      angle between the loading direction and the fault
                     tangent, degrees

Source is the Global Strain Rate Model v2.1 (Kreemer et al. 2014), which
publishes the strain-rate TENSOR per grid cell, so both magnitude and
orientation follow directly. Wang & Shen (2020) publish velocities only;
their strain figures are interpolated with a variable smoothing length
(40 to >250 km), which would make the smoothing our modelling choice.

From the tensor (exx, eyy, exy):

    gamma_max = sqrt( ((exx - eyy)/2)^2 + exy^2 )
    theta_p   = 0.5 * atan2(2*exy, exx - eyy)        principal axis
    max-shear planes lie at 45 deg to the principal axes

GSRM is in a geographic frame (x east, y north) and in units of 1e-9/yr;
both are converted here. The fault tangents are in the local rotated frame
(122.9 deg), so the tensor is rotated into that frame before the angle
relative to the fault tangent is taken.

All nine digitised polylines are sampled, not just the meshed chain.

Usage:
    python3 strain_rate_loading.py [--out-prefix xianshuihe_strain]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

from plot_fault_trace import (assign_fault_ids, chain, principal_frame,
                              read_kml, rotate, to_local_km)


def load_mesh_faults(case_dir: Path):
    """[(name, (n,2) node coords)] per fault, from a meshed case.

    Coordinates are already in the local rotated frame -- meshgen.py consumed
    the ftN.gmt.txt files written in that frame.
    """
    def find(name):
        for d in (case_dir, case_dir / "fem_mesh_output"):
            if (d / name).exists():
                return d / name
        raise SystemExit(f"error: {name} not found under {case_dir}")

    vert = np.loadtxt(find("vert.txt"))[:, :2]
    nsmp = np.loadtxt(find("nsmp.txt"), dtype=int)
    if nsmp.min() == 1:
        nsmp = nsmp - 1
    rows = [l.split() for l in open(find("meshGeneralInfo.txt")) if l.strip()]
    nfn = [int(x) for x in rows[1]]
    mx = max(nfn)
    return [(f"ft{i+1}", vert[nsmp[i * mx: i * mx + n, 0]]) for i, n in enumerate(nfn)]

SEC_PER_YR = 365.25 * 24 * 3600.0
GSRM_UNIT = 1.0e-9           # GSRM values are 1e-9/yr
R_EARTH_KM = 6371.0
MM = 1.0 / 25.4


def load_gsrm(path: Path):
    d = np.loadtxt(path)
    return dict(lat=d[:, 0], lon=d[:, 1], exx=d[:, 2], eyy=d[:, 3], exy=d[:, 4],
                e1=d[:, 8], e2=d[:, 9], azi_e1=d[:, 10])


def max_shear(exx, eyy, exy):
    """(magnitude, principal azimuth in rad) in the same frame as the input."""
    g = np.hypot(0.5 * (exx - eyy), exy)
    theta_p = 0.5 * np.arctan2(2.0 * exy, exx - eyy)
    return g, theta_p


def to_si(v_per_yr_1e9):
    return v_per_yr_1e9 * GSRM_UNIT / SEC_PER_YR


def main() -> None:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gsrm", type=Path,
                    default=here / "strain_rate_input" / "GSRM_v2.1_xianshuihe_region.txt")
    ap.add_argument("--out-prefix", default="xianshuihe_strain")
    ap.add_argument("--case", type=Path, default=None,
                    help="Meshed case dir. With it, sampling is at the MESH fault "
                         "nodes (what nsmpGeoPhys.txt needs); without it, at the "
                         "KML control points.")
    args = ap.parse_args()

    g = load_gsrm(args.gsrm)
    gamma_geo, theta_p_geo = max_shear(g["exx"], g["eyy"], g["exy"])
    gamma_si = to_si(gamma_geo)
    print(f"GSRM cells: {len(g['lat'])}")
    print(f"max shear strain rate (s^-1): min={gamma_si.min():.3e} "
          f"med={np.median(gamma_si):.3e} max={gamma_si.max():.3e}")
    print(f"  for scale, the SAF compsets use ftLoadMaxShear = 1.427e-14 s^-1")

    # fault geometry, local rotated frame
    segs = read_kml(here / "user_fault_geometry_input" / "xianshuihe_fault_trace.kml")
    lon0 = np.concatenate([s[:, 0] for s in segs]).mean()
    lat0 = np.concatenate([s[:, 1] for s in segs]).mean()
    local = [to_local_km(s, lon0, lat0) for s in segs]
    centre, theta = principal_frame(np.vstack(local))
    rotated = [rotate(p, centre, theta) for p in local]

    # GSRM cells in the same local rotated frame
    gx = np.radians(g["lon"] - lon0) * R_EARTH_KM * np.cos(np.radians(lat0))
    gy = np.radians(g["lat"] - lat0) * R_EARTH_KM
    gxy = rotate(np.column_stack([gx, gy]), centre, theta)
    # rotating the frame by -theta rotates tensor directions by -theta
    theta_p_loc = theta_p_geo - theta

    # Linear interpolation over the GSRM cells, with nearest-neighbour only to
    # fill outside the convex hull. Nearest-cell alone is wrong here: the grid
    # is 0.1 deg (~9-11 km) while mesh nodes are 400 m apart, so ~25
    # consecutive nodes would otherwise share one identical value.
    # Interpolate the TENSOR COMPONENTS, never the principal angle. Orientation
    # is pi-periodic, so interpolating theta across a wrap (+89 deg to -89 deg)
    # produces meaningless values -- it put +/-90 deg spikes through the
    # direction and compression panels. The components are single-valued
    # fields, and gamma and theta are derived from the interpolated tensor.
    # Rotation into the local frame is a fixed linear map, so it commutes with
    # interpolation and is applied after.
    comp = np.column_stack([to_si(g["exx"]), to_si(g["eyy"]), to_si(g["exy"])])
    interp_c = LinearNDInterpolator(gxy, comp)
    near_c = NearestNDInterpolator(gxy, comp)

    if args.case:
        faults = load_mesh_faults(args.case)
        source = f"mesh fault nodes ({args.case})"
    else:
        faults = []
        for fid, group in assign_fault_ids(rotated):
            pieces, _ = chain(rotated, group)
            faults.append((fid, np.vstack(pieces)))
        source = "KML control points"

    rows = []
    n_filled = 0
    print(f"\nsampling at {source}")
    print(f"{'fault':>6} {'nodes':>6} {'gamma_max (s^-1)':>18} {'load angle (deg)':>18}")
    for fid, pts in faults:
        c = interp_c(pts)
        bad = ~np.isfinite(c[:, 0])
        if bad.any():
            n_filled += int(bad.sum())
            c[bad] = near_c(pts[bad])
        exx_i, eyy_i, exy_i = c[:, 0], c[:, 1], c[:, 2]
        gam = np.hypot(0.5 * (exx_i - eyy_i), exy_i)
        thp = 0.5 * np.arctan2(2.0 * exy_i, exx_i - eyy_i) - theta
        tang = np.gradient(pts, axis=0)
        tang /= np.hypot(*tang.T)[:, None]
        tang_ang = np.arctan2(tang[:, 1], tang[:, 0])
        shear_dir = thp + np.pi / 4.0
        ang = np.degrees((shear_dir - tang_ang + np.pi / 2) % np.pi - np.pi / 2)
        print(f"{fid:>6} {len(pts):>6} {gam.mean():>18.3e} {ang.mean():>18.1f}")
        rows.append((fid, pts, gam, ang))
    if n_filled:
        print(f"  {n_filled} node(s) outside the GSRM hull, filled by nearest cell")

    csv = here / f"{args.out_prefix}_loading.csv"
    with open(csv, "w") as f:
        f.write("fault,node,x_km,y_km,gamma_max_per_s,load_angle_deg\n")
        for fid, pts, gam, ang in rows:
            for i, (p, gm, a) in enumerate(zip(pts, gam, ang)):
                f.write(f"{fid},{i},{p[0]:.4f},{p[1]:.4f},{gm:.6e},{a:.3f}\n")
    print(f"\nWrote {csv}")

    # ---- Figure-2-style panels ------------------------------------------
    # Mirrors Figure 2 of Liu et al. (2022): (a) on-fault max shear strain
    # rate, (b) its direction, (c) local fault strike, (d) angle of
    # compression = strike - shear direction, (e) fault geometry. That paper
    # took the strain field from Smith-Konter & Sandwell (2009) for the SAF;
    # GSRM v2.1 plays the same role here.
    #
    # Angles are azimuths in the local rotated frame, whose x axis runs along
    # the system -- the same reference the paper's "relative to the northern
    # end of the fault system" provides.
    plt.rcParams.update({"font.size": 12.5, "axes.titlesize": 13.5,
                         "axes.labelsize": 13, "xtick.labelsize": 11.5,
                         "ytick.labelsize": 11.5, "legend.fontsize": 11,
                         "savefig.bbox": "tight"})
    fig, axes = plt.subplots(5, 1, figsize=(190 * MM, 250 * MM), sharex=True)
    cols = plt.cm.tab10(np.linspace(0, 1, 10))
    letters = "abcde"

    for k, (fid, pts, gam, ang) in enumerate(rows):
        c = cols[k % 10]
        tang = np.gradient(pts, axis=0)
        tang /= np.hypot(*tang.T)[:, None]
        strike = np.degrees(np.arctan2(tang[:, 1], tang[:, 0]))
        shear_dir = strike + ang          # by construction of ang
        axes[0].plot(pts[:, 0], gam, "-", lw=2.0, color=c, label=fid)
        axes[1].plot(pts[:, 0], shear_dir, "-", lw=2.0, color=c)
        axes[2].plot(pts[:, 0], strike, "-", lw=2.0, color=c)
        axes[3].plot(pts[:, 0], ang, "-", lw=2.0, color=c)
        axes[4].plot(pts[:, 0], pts[:, 1], "-", lw=2.2, color=c)

    axes[0].axhline(1.427e-14, color="0.35", ls=":", lw=1.8)
    axes[0].set_yscale("log")
    axes[0].set_ylabel("$\\dot\\gamma_{max}$ (s$^{-1}$)")
    axes[0].set_title("On-fault maximum shear strain rate (GSRM v2.1); "
                      "dotted = SAF compset 1.427e-14")
    axes[0].legend(ncol=7, frameon=False, fontsize=10.5, loc="lower center")
    axes[1].set_ylabel("shear dir. ($^\\circ$)")
    axes[1].set_title("Direction of maximum shear strain rate")
    axes[2].set_ylabel("fault strike ($^\\circ$)")
    axes[2].set_title("Local fault strike")
    axes[3].axhline(0, color="0.5", lw=1.0)
    axes[3].set_ylabel("compression ($^\\circ$)")
    axes[3].set_title("Angle of compression = strike $-$ shear direction "
                      "(this is ftLoadAngle)")
    axes[4].set_ylabel("y (km)")
    axes[4].set_title("Fault geometry, local rotated frame")
    axes[4].set_aspect("equal")
    axes[4].set_xlabel("Along-strike distance (km)")
    for k, ax in enumerate(axes):
        ax.grid(alpha=0.25, lw=0.7)
        ax.text(0.010, 0.94, f"({letters[k]})", transform=ax.transAxes,
                fontsize=15, fontweight="bold", va="top")
    out = here / f"{args.out_prefix}_figure2.png"
    fig.savefig(out, dpi=250)
    print(f"Wrote {out}")

    # ---- map of the strain field ----------------------------------------
    fig2 = plt.figure(figsize=(190 * MM, 150 * MM))
    ax_m = fig2.add_subplot(111)
    sc = ax_m.scatter(g["lon"], g["lat"], c=gamma_si, s=9, cmap="magma",
                      norm=matplotlib.colors.LogNorm())
    cb = fig2.colorbar(sc, ax=ax_m, pad=0.02, fraction=0.035)
    cb.set_label("max shear strain rate (s$^{-1}$)")
    step = 13
    ang_geo = theta_p_geo + np.pi / 4.0
    ax_m.quiver(g["lon"][::step], g["lat"][::step],
                np.cos(ang_geo[::step]), np.sin(ang_geo[::step]),
                headwidth=0, headlength=0, headaxislength=0, pivot="mid",
                scale=32, width=0.0026, color="0.15")
    for s_ in segs:
        ax_m.plot(s_[:, 0], s_[:, 1], "-", lw=2.0, color="tab:cyan")
    ax_m.set_aspect(1.0 / np.cos(np.radians(lat0)))
    ax_m.set_xlabel("Longitude ($^\\circ$E)")
    ax_m.set_ylabel("Latitude ($^\\circ$N)")
    ax_m.set_title("GSRM v2.1 max shear strain rate; ticks give its direction")
    out2 = here / f"{args.out_prefix}_map.png"
    fig2.savefig(out2, dpi=250)
    print(f"Wrote {out2}")


if __name__ == "__main__":
    main()
