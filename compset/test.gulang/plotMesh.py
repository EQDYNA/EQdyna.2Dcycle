#!/usr/bin/env python3
"""
plotMesh.py — Standalone mesh visualization for EQdyna.2Dcycle.

Reads output files from fem_mesh_output/ and fault geometry from
user_fault_geometry_input/, then generates vector PDF figures.
Each panel is a separate figure file.

Usage:
    cd compset/test.gulang
    python3 plotMesh.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import matplotlib.colors as mcolors

# ============================================================
# Configuration
# ============================================================
output_dir = 'fem_mesh_output'
plot_dir = os.path.join(output_dir, 'meshPlots')
os.makedirs(plot_dir, exist_ok=True)
fault_input_dir = 'user_fault_geometry_input'
ftNames = ['ft1', 'ft2', 'ft3', 'ft4', 'ft5']

# ============================================================
# Data loading functions
# ============================================================

def load_mesh(output_dir):
    """Load node coordinates and element connectivity."""
    points = np.loadtxt(os.path.join(output_dir, 'vert.txt'))
    cells  = np.loadtxt(os.path.join(output_dir, 'fac.txt'), dtype=int)
    return points, cells


def load_mesh_info(output_dir):
    """Load mesh metadata from meshGeneralInfo.txt."""
    with open(os.path.join(output_dir, 'meshGeneralInfo.txt')) as f:
        lines = f.readlines()
    toks0 = lines[0].split()
    totalNodes = int(toks0[0])
    totalElems = int(toks0[1])
    faultNodeCounts = [int(x) for x in lines[1].split()]
    toks2 = lines[2].split()
    bounds = {
        'xmin': float(toks2[0]),
        'xmax': float(toks2[1]),
        'ymin': float(toks2[2]),
        'ymax': float(toks2[3]),
    }
    return totalNodes, totalElems, faultNodeCounts, bounds


def load_fault_traces(fault_input_dir, ftNames):
    """Load fault trace coordinates from .gmt.txt files."""
    traces = {}
    for ftName in ftNames:
        fname = os.path.join(fault_input_dir, ftName + '.gmt.txt')
        data = np.loadtxt(fname)
        traces[ftName] = (data[:, 0], data[:, 1])
    return traces


def load_split_nodes(output_dir, ftNames, faultNodeCounts):
    """Load master/slave node-id pairs from nsmp.txt."""
    nsmp = np.loadtxt(os.path.join(output_dir, 'nsmp.txt'), dtype=int)
    nFt = len(ftNames)
    maxN = nsmp.shape[0] // nFt

    master_ids = {}
    slave_ids  = {}
    for iFt, ftName in enumerate(ftNames):
        n = faultNodeCounts[iFt]
        block = nsmp[iFt * maxN : iFt * maxN + n, :]
        master_ids[ftName] = block[:, 0]
        slave_ids[ftName]  = block[:, 1]
    return master_ids, slave_ids


def load_tangent_data(output_dir, ftNames, faultNodeCounts):
    """Load tangent vectors from nsmpTanlen.txt."""
    data = np.loadtxt(os.path.join(output_dir, 'nsmpTanlen.txt'))
    nFt = len(ftNames)
    maxN = data.shape[0] // nFt

    tangents = {}
    for iFt, ftName in enumerate(ftNames):
        n = faultNodeCounts[iFt]
        block = data[iFt * maxN : iFt * maxN + n, :]
        tangents[ftName] = (block[:, 0], block[:, 1])
    return tangents

# ============================================================
# Geometry helpers
# ============================================================

def compute_aspect_ratio(points, cells):
    """Compute aspect ratio (max_edge / min_edge) for each quad element."""
    verts = points[cells]
    edges = np.empty((len(cells), 4))
    for i in range(4):
        j = (i + 1) % 4
        diff = verts[:, j, :] - verts[:, i, :]
        edges[:, i] = np.sqrt(diff[:, 0]**2 + diff[:, 1]**2)
    max_e = edges.max(axis=1)
    min_e = edges.min(axis=1)
    min_e = np.maximum(min_e, 1e-30)
    return max_e / min_e


def get_fault_bbox(fault_traces, buffer_km):
    """Bounding box around all fault traces with buffer (km)."""
    all_x = np.concatenate([t[0] for t in fault_traces.values()])
    all_y = np.concatenate([t[1] for t in fault_traces.values()])
    return (all_x.min() - buffer_km, all_x.max() + buffer_km,
            all_y.min() - buffer_km, all_y.max() + buffer_km)


def get_single_fault_bbox(fault_traces, ftName, buffer_km):
    """Bounding box around one fault trace with buffer."""
    x, y = fault_traces[ftName]
    return (x.min() - buffer_km, x.max() + buffer_km,
            y.min() - buffer_km, y.max() + buffer_km)

# ============================================================
# Plotting helpers
# ============================================================

FAULT_COLORS = ['#e6194b', '#3cb44b', '#4363d8', '#f58231', '#911eb4']


def _make_poly_collection(points, cells, **kwargs):
    """Return a PolyCollection of quads for fast rendering."""
    verts = points[cells]
    defaults = dict(edgecolors='black', linewidths=0.15,
                    facecolors='none')
    defaults.update(kwargs)
    return PolyCollection(verts, **defaults)


def _draw_fault_traces(ax, fault_traces, ftNames, lw=2.5, legend=False):
    """Overlay bold fault-trace lines on *ax*."""
    handles = []
    for i, ftName in enumerate(ftNames):
        x, y = fault_traces[ftName]
        line, = ax.plot(x, y, color=FAULT_COLORS[i % len(FAULT_COLORS)],
                        linewidth=lw, solid_capstyle='round',
                        label=ftName, zorder=5)
        handles.append(line)
    if legend:
        ax.legend(handles=handles, loc='best', fontsize=8)


def _set_bbox(ax, bbox):
    """Set axis limits from (xmin, xmax, ymin, ymax)."""
    ax.set_xlim(bbox[0], bbox[1])
    ax.set_ylim(bbox[2], bbox[3])


def _save_pdf(fig, filename):
    """Save figure as PDF and close."""
    outpath = os.path.join(plot_dir, filename)
    fig.savefig(outpath, format='pdf', bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved {outpath}')


# ============================================================
# Figure 1a: Full domain mesh
# ============================================================

def plot_full_domain(points, cells, fault_traces, ftNames, bounds):
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.add_collection(_make_poly_collection(points, cells))
    _draw_fault_traces(ax, fault_traces, ftNames, lw=2.0)
    ax.set_xlim(bounds['xmin'], bounds['xmax'])
    ax.set_ylim(bounds['ymin'], bounds['ymax'])
    ax.set_aspect('equal')
    ax.set_title('Full Domain Mesh')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    fig.tight_layout()
    _save_pdf(fig, 'fig1a_full_domain.pdf')


# ============================================================
# Figure 1b: Zoomed fault region
# ============================================================

def plot_fault_region_zoom(points, cells, fault_traces, ftNames):
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.add_collection(_make_poly_collection(points, cells))
    _draw_fault_traces(ax, fault_traces, ftNames, lw=2.5)
    bbox = get_fault_bbox(fault_traces, buffer_km=5.0)
    _set_bbox(ax, bbox)
    ax.set_aspect('equal')
    ax.set_title('Fault Region (5 km buffer)')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    fig.tight_layout()
    _save_pdf(fig, 'fig1b_fault_region_zoom.pdf')


# ============================================================
# Figure 2: Per-fault zoom — one figure per fault
# ============================================================

def plot_per_fault_zoom(points, cells, fault_traces, ftNames):
    for idx, ftName in enumerate(ftNames):
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.add_collection(_make_poly_collection(points, cells))
        # highlight this fault thick, others thin grey
        for j, fn in enumerate(ftNames):
            x, y = fault_traces[fn]
            if fn == ftName:
                ax.plot(x, y, color=FAULT_COLORS[j % len(FAULT_COLORS)],
                        linewidth=3.0, solid_capstyle='round', zorder=5,
                        label=fn)
            else:
                ax.plot(x, y, color='grey', linewidth=0.5, alpha=0.4,
                        zorder=4)
        bbox = get_single_fault_bbox(fault_traces, ftName, buffer_km=2.0)
        _set_bbox(ax, bbox)
        ax.set_aspect('equal')
        ax.set_title(f'{ftName} — mesh detail', fontsize=12, fontweight='bold')
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        fig.tight_layout()
        _save_pdf(fig, f'fig2_{ftName}_zoom.pdf')


# ============================================================
# Figure 3a: Element quality — full domain
# ============================================================

def plot_quality_full(points, cells, fault_traces, ftNames, bounds, ar):
    vmin, vmax = 1.0, 5.0
    cmap = plt.cm.RdYlGn_r

    fig, ax = plt.subplots(figsize=(10, 8))
    verts = points[cells]
    facecolors = cmap((np.clip(ar, vmin, vmax) - vmin) / (vmax - vmin))
    pc = PolyCollection(verts, edgecolors='face', linewidths=0.1,
                        facecolors=facecolors)
    ax.add_collection(pc)
    _draw_fault_traces(ax, fault_traces, ftNames, lw=1.5)
    ax.set_xlim(bounds['xmin'], bounds['xmax'])
    ax.set_ylim(bounds['ymin'], bounds['ymax'])
    ax.set_aspect('equal')
    ax.set_title('Element Quality — Full Domain')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')

    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.7)
    cbar.set_label('Aspect Ratio (max edge / min edge)')

    fig.tight_layout()
    _save_pdf(fig, 'fig3a_quality_full.pdf')


# ============================================================
# Figure 3b: Element quality — fault region zoom
# ============================================================

def plot_quality_zoom(points, cells, fault_traces, ftNames, ar):
    vmin, vmax = 1.0, 5.0
    cmap = plt.cm.RdYlGn_r

    fig, ax = plt.subplots(figsize=(14, 6))
    verts = points[cells]
    facecolors = cmap((np.clip(ar, vmin, vmax) - vmin) / (vmax - vmin))
    pc = PolyCollection(verts, edgecolors='face', linewidths=0.1,
                        facecolors=facecolors)
    ax.add_collection(pc)
    _draw_fault_traces(ax, fault_traces, ftNames, lw=1.5)
    bbox = get_fault_bbox(fault_traces, buffer_km=5.0)
    _set_bbox(ax, bbox)
    ax.set_aspect('equal')
    ax.set_title('Element Quality — Fault Region')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')

    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.7)
    cbar.set_label('Aspect Ratio (max edge / min edge)')

    fig.tight_layout()
    _save_pdf(fig, 'fig3b_quality_fault_zoom.pdf')


# ============================================================
# Figure 4: Fault trace overlay on mesh with legend
# ============================================================

def plot_fault_overlay(points, cells, fault_traces, ftNames):
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.add_collection(_make_poly_collection(points, cells,
                                            edgecolors='grey',
                                            linewidths=0.1))
    _draw_fault_traces(ax, fault_traces, ftNames, lw=3.0, legend=True)
    bbox = get_fault_bbox(fault_traces, buffer_km=5.0)
    _set_bbox(ax, bbox)
    ax.set_aspect('equal')
    ax.set_title('Fault Traces on Mesh')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    fig.tight_layout()
    _save_pdf(fig, 'fig4_fault_overlay.pdf')


# ============================================================
# Figure 5: Split-node verification — one figure per fault
# ============================================================

def plot_split_nodes(points, cells, fault_traces, ftNames,
                     master_ids, slave_ids, tangents):
    offset_km = 0.3

    for idx, ftName in enumerate(ftNames):
        fig, ax = plt.subplots(figsize=(12, 6))

        # mesh background
        ax.add_collection(_make_poly_collection(points, cells,
                                                edgecolors='lightgrey',
                                                linewidths=0.1))
        # fault trace
        fx, fy = fault_traces[ftName]
        ax.plot(fx, fy, color=FAULT_COLORS[idx % len(FAULT_COLORS)],
                linewidth=2.0, zorder=4)

        # master nodes
        mids = master_ids[ftName]
        mx = points[mids, 0]
        my = points[mids, 1]
        ax.scatter(mx, my, s=18, c='blue', marker='o', zorder=6,
                   label='master', edgecolors='navy', linewidths=0.3)

        # slave nodes — offset perpendicular to local fault tangent
        tanX, tanY = tangents[ftName]
        nx = -tanY
        ny =  tanX
        sx = points[slave_ids[ftName], 0] + offset_km * nx
        sy = points[slave_ids[ftName], 1] + offset_km * ny
        ax.scatter(sx, sy, s=24, c='red', marker='*', zorder=6,
                   label='slave (offset)', edgecolors='darkred',
                   linewidths=0.3)

        bbox = get_single_fault_bbox(fault_traces, ftName, buffer_km=2.0)
        _set_bbox(ax, bbox)
        ax.set_aspect('equal')
        ax.set_title(f'{ftName} — split nodes', fontsize=12,
                     fontweight='bold')
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.legend(fontsize=8, loc='best')
        fig.tight_layout()
        _save_pdf(fig, f'fig5_{ftName}_split_nodes.pdf')


# ============================================================
# Main
# ============================================================

def main():
    print('Loading mesh data ...')
    points, cells = load_mesh(output_dir)
    totalNodes, totalElems, faultNodeCounts, bounds = load_mesh_info(output_dir)
    fault_traces = load_fault_traces(fault_input_dir, ftNames)
    master_ids, slave_ids = load_split_nodes(output_dir, ftNames,
                                             faultNodeCounts)
    tangents = load_tangent_data(output_dir, ftNames, faultNodeCounts)

    print(f'  {totalNodes} nodes, {totalElems} elements, '
          f'{len(ftNames)} faults')

    print('Figure 1: Overview ...')
    plot_full_domain(points, cells, fault_traces, ftNames, bounds)
    plot_fault_region_zoom(points, cells, fault_traces, ftNames)

    print('Figure 2: Per-fault zoom ...')
    plot_per_fault_zoom(points, cells, fault_traces, ftNames)

    print('Figure 3: Element quality ...')
    ar = compute_aspect_ratio(points, cells)
    plot_quality_full(points, cells, fault_traces, ftNames, bounds, ar)
    plot_quality_zoom(points, cells, fault_traces, ftNames, ar)

    print('Figure 4: Fault overlay ...')
    plot_fault_overlay(points, cells, fault_traces, ftNames)

    print('Figure 5: Split-node verification ...')
    plot_split_nodes(points, cells, fault_traces, ftNames,
                     master_ids, slave_ids, tangents)

    print('Done — all figures saved to', plot_dir)


if __name__ == '__main__':
    main()
