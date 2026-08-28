#!/usr/bin/env python3
"""Write the GSRM-derived loading into a meshed case's nsmpGeoPhys.txt.

meshgen.py writes nsmpGeoPhys.txt during meshing, but the strain sampling
needs the mesh to already exist, so the loading columns are patched in
afterwards rather than re-entering meshgen. Row order is the nsmp order, the
same order strain_rate_loading.py --case writes its CSV in, so the mapping is
by construction and not by coordinate matching.

Columns rewritten (1-based, per the C_mesh=3 convention):

    6  ftLoadMaxShear   per-node max shear strain rate, s^-1
    7  ftLoadAngle      per-node angle between loading and fault tangent, deg
    9  ftVis            per-node viscosity, Pa s

ftVis follows the SAF construction. In interstress.f90 the SAF sets
ant = ant0*str/rd, so the product (loading rate x viscosity) is CONSTANT
across nodes: the asymptotic shear stress is uniform and only the angle
varies, while the local rate sets the approach timescale. Here that means

    ftVis_i = TARGET_SHEAR_STRESS_PA / gamma_i

Leaving ftVis uniform while ftLoadMaxShear varies would NOT reproduce that
construction -- the asymptotic stress would then vary by the full range of
the strain field.

TARGET_SHEAR_STRESS_PA is the tunable. The SAF uses 1.427e-14 x 8.4e21 =
120 MPa; gulang needed less (1.427e-14 x 5.0e21 = 71 MPa) because 8.4e21
drove the normal stress tensile at large strike angles. Eastern Tibet has
the same 100 MPa ambient normal stress as gulang, so gulang's value is the
safer starting point.

Usage:
    python3 apply_strain_loading.py --case ../../work/xsh_mesh7
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

# 1.427e-14 s^-1 x 5.0e21 Pa s, i.e. gulang's asymptotic shear stress.
TARGET_SHEAR_STRESS_PA = 1.427e-14 * 5.0e21

COL_MAXSHEAR, COL_ANGLE, COL_VIS = 5, 6, 8   # 0-based


def main() -> None:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", type=Path, required=True, help="Meshed case directory")
    ap.add_argument("--csv", type=Path, default=here / "xianshuihe_strain_loading.csv")
    ap.add_argument("--target-stress", type=float, default=TARGET_SHEAR_STRESS_PA,
                    help="Asymptotic shear stress in Pa that gamma*ftVis is held at")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    def find(name: str) -> Path:
        for d in (args.case, args.case / "fem_mesh_output"):
            if (d / name).exists():
                return d / name
        raise SystemExit(f"error: {name} not found under {args.case}")

    rows = list(csv.DictReader(open(args.csv)))
    if not rows:
        raise SystemExit(f"error: {args.csv} is empty")

    info = [l.split() for l in open(find("meshGeneralInfo.txt")) if l.strip()]
    nfn = [int(x) for x in info[1]]
    maxftnode = max(nfn)

    gp_path = find("nsmpGeoPhys.txt")
    gp = np.loadtxt(gp_path)
    if gp.shape[0] != maxftnode * len(nfn):
        raise SystemExit(f"error: nsmpGeoPhys.txt has {gp.shape[0]} rows, expected "
                         f"{maxftnode*len(nfn)} (maxftnode {maxftnode} x {len(nfn)} faults)")

    # group CSV rows by fault, preserving order
    per_fault: dict[str, list[dict]] = {}
    for r in rows:
        per_fault.setdefault(r["fault"], []).append(r)
    if len(per_fault) != len(nfn):
        raise SystemExit(f"error: CSV has {len(per_fault)} faults, mesh has {len(nfn)}")

    print(f"case   : {args.case}")
    print(f"target : gamma * ftVis = {args.target_stress:.3e} Pa "
          f"({args.target_stress/1e6:.0f} MPa asymptotic shear stress)")
    print(f"{'fault':>6} {'nodes':>6} {'gamma range (1/s)':>26} {'ftVis range (Pa s)':>26}")
    for i, (fid, frows) in enumerate(sorted(per_fault.items(),
                                            key=lambda kv: int(kv[0][2:]))):
        if len(frows) != nfn[i]:
            raise SystemExit(f"error: {fid} has {len(frows)} CSV rows, mesh has {nfn[i]}")
        gam = np.array([float(r["gamma_max_per_s"]) for r in frows])
        ang = np.array([float(r["load_angle_deg"]) for r in frows])
        vis = args.target_stress / gam
        sl = slice(i * maxftnode, i * maxftnode + nfn[i])
        gp[sl, COL_MAXSHEAR] = gam
        gp[sl, COL_ANGLE] = ang
        gp[sl, COL_VIS] = vis
        print(f"{fid:>6} {nfn[i]:>6} {gam.min():>12.3e}-{gam.max():<13.3e} "
              f"{vis.min():>12.3e}-{vis.max():<13.3e}")

    if args.dry_run:
        print("\ndry run, nothing written")
        return
    np.savetxt(gp_path, gp, fmt="%.6e")
    print(f"\nWrote {gp_path}")
    other = args.case / "nsmpGeoPhys.txt"
    if gp_path != other and other.exists():
        np.savetxt(other, gp, fmt="%.6e")
        print(f"Wrote {other}")
if __name__ == "__main__":
    main()
