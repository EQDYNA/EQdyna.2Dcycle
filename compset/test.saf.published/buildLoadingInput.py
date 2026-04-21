#! /usr/bin/env python3
"""
Build 4-column fault geometry input files for test.saf (C_mesh=3) by mapping
the 2.0.2 Rate_direction.txt loading values to the controlling points.

Uses the actual 2.0.2 fault mesh node positions (vert.txt + nsmp.txt) to do
nearest-neighbor lookup for each controlling point.

Output: compset/test.saf/user_fault_geometry_input/<name>.gmt.txt  (4 cols)
        x(km)  y(km)  loading_rate(s^-1)  loading_angle(deg)
"""
import numpy as np
from scipy.spatial import KDTree
import os

ROOT = os.path.join(os.path.dirname(__file__), '..', '..')

# ---- input paths ----
mesh_dir  = os.path.join(ROOT, 'work', 'test.saf.2.0.2')
pub_dir   = os.path.join(ROOT, 'compset', 'test.saf.published')
out_dir   = os.path.join(ROOT, 'compset', 'test.saf', 'user_fault_geometry_input')

# ---- read 2.0.2 mesh ----
vert = np.loadtxt(os.path.join(mesh_dir, 'vert.txt'))   # (numnp, 2) in km
nsmp = np.loadtxt(os.path.join(mesh_dir, 'nsmp.txt'), dtype=int)  # (maxftnode*ntotft, 2)

# ---- read Rate_direction.txt ----
rd = np.loadtxt(os.path.join(pub_dir, 'Rate_direction.txt'))  # (maxftnode*ntotft, 2)
ntotft   = 3
maxftnode = rd.shape[0] // ntotft   # = 1769
nfnode   = []
for k in range(ntotft):
    blk = rd[k*maxftnode:(k+1)*maxftnode, 0]
    nfnode.append(int(np.sum(blk != 0)))
print(f'maxftnode={maxftnode}, nfnode={nfnode}')

# ---- fault blocks: slave node positions + rate/angle per fault ----
# nsmp rows are 1-indexed (Fortran); convert to 0-indexed for vert lookup
# nsmp[:,0] = slave node (1-indexed in Fortran convention from 2.0.2 nsmp.txt)
# Check if nsmp is 0-indexed or 1-indexed
print('nsmp first few rows:', nsmp[:3])
# The 2.0.2 nsmp.txt from meshgen1 is 1-indexed → subtract 1 for Python
# (if first entry is 0, it's already 0-indexed)
offset = 0 if nsmp[0, 0] == 0 else -1

fault_names = ['sjfn', 'sjfs', 'ssaf']   # order matches Rate_direction.txt faults 1,2,3
gmt_files = {
    'sjfn': ['sjfn.gmt.txt'],
    'sjfs': ['sjfs.gmt.txt'],
    'ssaf': ['ssaf1.gmt.txt', 'ssaf2.gmt.txt'],
}

for k, fname in enumerate(fault_names):
    # Fault node positions from nsmp slave nodes
    n = nfnode[k]
    slave_idx = nsmp[k*maxftnode : k*maxftnode + n, 0] + offset  # 0-indexed into vert
    ft_xy = vert[slave_idx, :] / 1000.0   # m → km

    # Rate/angle for this fault's nodes
    blk_rate  = rd[k*maxftnode : k*maxftnode + n, 0]
    blk_angle = rd[k*maxftnode : k*maxftnode + n, 1]

    # KD-tree over fault nodes
    tree = KDTree(ft_xy)

    # Process each controlling point file for this fault
    for gmt_file in gmt_files[fname]:
        ctrl = np.loadtxt(os.path.join(out_dir, gmt_file))
        ctrl_xy = ctrl[:, :2]   # always use first 2 cols (x,y in km)
        dists, idxs = tree.query(ctrl_xy)
        rate_at_ctrl  = blk_rate[idxs]
        angle_at_ctrl = blk_angle[idxs]

        out = np.column_stack([ctrl_xy, rate_at_ctrl, angle_at_ctrl])
        out_path = os.path.join(out_dir, gmt_file)
        np.savetxt(out_path, out, fmt='%.10e')
        print(f'{gmt_file}: {len(ctrl)} pts, max dist={dists.max():.2f} km, '
              f'rate=[{rate_at_ctrl.min():.3e},{rate_at_ctrl.max():.3e}], '
              f'angle=[{angle_at_ctrl.min():.2f},{angle_at_ctrl.max():.2f}] deg')
