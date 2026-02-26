"""
MOOSE Quick-Start: Visualization Script for All 36 Cases
=========================================================
Generates 2-3 plots per case and saves each PNG into that case's own directory.
For example, Case 01's plot is written to case01-1d-steady-diffusion/case01_diffusion_1d.png.

Requirements: matplotlib, numpy, netCDF4
Usage:
    python visualize_all.py       (run from quickstart-runs/ directory)

All paths are resolved relative to the directory that contains this script.
"""

import csv
import math
import os
import sys
import traceback

import matplotlib

matplotlib.use("Agg")  # Non-interactive backend — write PNG without a display
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers the 3D projection)

try:
    import netCDF4 as nc
except ImportError:
    print("ERROR: netCDF4 is not installed. Run: pip install netCDF4")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Visual style constants
# ---------------------------------------------------------------------------
FIGSIZE_SINGLE = (7, 5)
FIGSIZE_WIDE = (10, 4)
FIGSIZE_SQUARE = (6, 6)
FIGSIZE_LARGE = (11, 8)
DPI = 150
CMAP_SCALAR = "viridis"
CMAP_TEMP = "coolwarm"
FONT_TITLE = 13
FONT_AXIS = 11

plt.rcParams.update(
    {
        "font.size": FONT_AXIS,
        "axes.titlesize": FONT_TITLE,
        "axes.labelsize": FONT_AXIS,
        "figure.dpi": DPI,
    }
)

# ---------------------------------------------------------------------------
# Tracking
# ---------------------------------------------------------------------------
generated_plots: list[str] = []
skipped_cases: list[str] = []


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

CASE_DIRS = {
    "01": "case01-1d-steady-diffusion",
    "02": "case02-2d-steady-diffusion",
    "03": "case03-transient-heat",
    "04": "case04-manufactured-solution",
    "05": "case05-neumann-bc",
    "06": "case06-two-material-domain",
    "07": "case07-nonlinear-diffusion",
    "08": "case08-advection-diffusion",
    "09": "case09-coupled-reaction-diffusion",
    "10": "case10-adaptive-mesh-refinement",
    "11": "case11-adaptive-timestepping",
    "12": "case12-multiapp-coupling",
    "13": "case13-custom-kernel",
    "14": "case14-thermoelasticity",
    "15": "case15-lid-driven-cavity",
    "16": "case16-natural-convection",
    "17": "case17-joule-heating",
    "18": "case18-cahn-hilliard",
    "19": "case19-porous-flow",
    "20": "case20-elastic-wave",
    "21": "case21-bimetallic-strip",
    "22": "case22-charge-relaxation",
    "23": "case23-magnetic-diffusion",
    "24": "case24-drift-diffusion",
    "25": "case25-induction-heating",
    "26": "case26-ehd-pumping",
    "27": "case27-hartmann-flow",
    "28": "case28-twoway-joule-heating",
    "29": "case29-electroconvection",
    "30": "case30-waveguide-cutoff",
    "31": "case31-driven-cavity",
    "32": "case32-dielectric-slab",
    "33": "case33-coupled-resonators",
    "34": "case34-thermal-noise",
    "35": "case35-dispersive-pulse",
    "36": "case36-soliton-pulse",
}


def case_path(casenum: str, filename: str) -> str:
    """Return absolute path for a file in a case directory."""
    return os.path.join(SCRIPT_DIR, CASE_DIRS[casenum], filename)


def save_fig(fig: plt.Figure, casenum: str, name: str) -> None:
    """Save figure into its case directory, record it, and close it."""
    path = case_path(casenum, name)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    generated_plots.append(path)
    print(f"    saved: {CASE_DIRS[casenum]}/{name}")


def file_exists(*paths: str) -> bool:
    """Return True only when every listed path exists."""
    return all(os.path.isfile(p) for p in paths)


def open_exodus(path: str):
    """Open an Exodus II file (netCDF4 Dataset)."""
    return nc.Dataset(path)


def get_nod_var_names(ds) -> list[str]:
    """Decode the list of nodal variable names stored in the Exodus file."""
    nv = ds.variables["name_nod_var"]
    names = []
    for i in range(nv.shape[0]):
        row = nv[i]
        name = b"".join(row.compressed()).decode("utf-8", errors="ignore").strip()
        names.append(name)
    return names


def get_nod_var(ds, index_1based: int, timestep: int = -1) -> np.ndarray:
    """
    Return nodal variable values.

    Parameters
    ----------
    ds            : open netCDF4 Dataset
    index_1based  : variable index starting at 1 (vals_nod_var1, vals_nod_var2 …)
    timestep      : time-step index; -1 means the last available step
    """
    key = f"vals_nod_var{index_1based}"
    data = ds.variables[key][timestep, :]
    return np.array(data, dtype=float)


def get_coords_2d(ds):
    """Return (x, y) coordinate arrays as float64 numpy arrays."""
    x = np.array(ds.variables["coordx"][:], dtype=float)
    y = np.array(ds.variables["coordy"][:], dtype=float)
    return x, y


def get_times(ds) -> np.ndarray:
    """Return time values as a float64 array."""
    return np.array(ds.variables["time_whole"][:], dtype=float)


def read_csv(path: str) -> dict[str, list[float]]:
    """
    Read a MOOSE postprocessor CSV file.

    Returns a dict mapping column name -> list of floats.
    """
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    result: dict[str, list[float]] = {}
    if not rows:
        return result
    for key in rows[0]:
        result[key] = [float(r[key]) for r in rows]
    return result


def find_closest_timestep(times: np.ndarray, target: float) -> int:
    """Return the index of the time value closest to `target`."""
    return int(np.argmin(np.abs(times - target)))


def get_elem_var_names(ds) -> list[str]:
    """Decode the list of element variable names stored in the Exodus file."""
    ev = ds.variables["name_elem_var"]
    return list(nc.chartostring(ev[:]))


def get_elem_var(ds, index_1based: int, block: int = 1, timestep: int = -1) -> np.ndarray:
    """
    Return element variable values for a given element block.

    Parameters
    ----------
    ds            : open netCDF4 Dataset
    index_1based  : variable index starting at 1
    block         : element block number (1-based)
    timestep      : time-step index; -1 means the last available step
    """
    key = f"vals_elem_var{index_1based}eb{block}"
    data = ds.variables[key][timestep, :]
    return np.array(data, dtype=float)


def get_elem_centroids_2d(ds, block: int = 1):
    """
    Compute element centroid coordinates for a 2D mesh.

    Returns (cx, cy) arrays where each entry is the average of the
    element's node coordinates.
    """
    x = np.array(ds.variables["coordx"][:], dtype=float)
    y = np.array(ds.variables["coordy"][:], dtype=float)
    conn_key = f"connect{block}"
    conn = ds.variables[conn_key][:] - 1  # Exodus uses 1-based indexing
    cx = np.mean(x[conn], axis=1)
    cy = np.mean(y[conn], axis=1)
    return cx, cy


def nodes_near_y(y: np.ndarray, y_target: float, tol: float = 0.02) -> np.ndarray:
    """Return integer indices of nodes within `tol` of `y_target`."""
    return np.where(np.abs(y - y_target) < tol)[0]


def add_2d_contour(ax, x, y, values, levels=20, cmap=CMAP_SCALAR, label=""):
    """Draw a tricontourf on ax, add a colorbar, and return the mappable."""
    tcf = ax.tricontourf(x, y, values, levels=levels, cmap=cmap)
    cb = plt.colorbar(tcf, ax=ax)
    if label:
        cb.set_label(label)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    return tcf


# ---------------------------------------------------------------------------
# Case implementations
# ---------------------------------------------------------------------------

def plot_case01():
    """Case 01: 1D Steady Diffusion — u(x) = x"""
    print("Case 01: 1D Steady Diffusion")
    efile = case_path("01", "case01_diffusion_1d_out.e")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case01")
        return

    ds = open_exodus(efile)
    x = np.array(ds.variables["coordx"][:], dtype=float)
    u = get_nod_var(ds, 1, timestep=-1)
    ds.close()

    order = np.argsort(x)
    x_sorted = x[order]
    u_sorted = u[order]

    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
    ax.plot(x_sorted, u_sorted, "b-o", markersize=4, label="MOOSE solution")
    ax.plot(x_sorted, x_sorted, "r--", linewidth=2, label="Exact: u = x")
    ax.set_xlabel("x")
    ax.set_ylabel("u")
    ax.set_title("Case 01: 1D Steady Diffusion — u(x) = x")
    ax.legend()
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    save_fig(fig, "01", "case01_diffusion_1d.png")


def plot_case02():
    """Case 02: 2D Steady Diffusion — 2D contour and 3D surface"""
    print("Case 02: 2D Steady Diffusion")
    efile = case_path("02", "case02_diffusion_2d_out.e")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case02")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    u = get_nod_var(ds, 1, timestep=-1)
    ds.close()

    # Plot 1: 2D contour
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    add_2d_contour(ax, x, y, u, cmap=CMAP_SCALAR, label="u")
    ax.set_title("Case 02: 2D Steady Diffusion — u(x,y)")
    fig.tight_layout()
    save_fig(fig, "02", "case02_contour_2d.png")

    # Plot 2: 3D surface
    fig = plt.figure(figsize=FIGSIZE_SQUARE)
    ax3 = fig.add_subplot(111, projection="3d")
    ax3.plot_trisurf(x, y, u, cmap=CMAP_SCALAR, linewidth=0, antialiased=False)
    ax3.set_xlabel("x")
    ax3.set_ylabel("y")
    ax3.set_zlabel("u")
    ax3.set_title("Case 02: 3D Surface — u(x,y) = x")
    fig.tight_layout()
    save_fig(fig, "02", "case02_surface_3d.png")


def plot_case03():
    """Case 03: Transient Heat Equation"""
    print("Case 03: Transient Heat Equation")
    efile = case_path("03", "case03_heat_transient_out.e")
    cfile = case_path("03", "case03_heat_transient_out.csv")

    if not file_exists(cfile, efile):
        print("    SKIP: output file(s) not found")
        skipped_cases.append("case03")
        return

    data = read_csv(cfile)
    time = data["time"]
    avg_T = data["avg_temperature"]
    max_T = data["max_temperature"]

    # Plot 1: CSV time history
    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
    ax.plot(time, avg_T, "b-", label="avg_temperature")
    ax.plot(time, max_T, "r-", label="max_temperature")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Temperature")
    ax.set_title("Case 03: Temperature Evolution")
    ax.legend()
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    save_fig(fig, "03", "case03_temperature_history.png")

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    times = get_times(ds)

    # Plot 2: 2D heatmap at final timestep
    T_final = get_nod_var(ds, 1, timestep=-1)
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    add_2d_contour(ax, x, y, T_final, cmap=CMAP_TEMP, label="T")
    ax.set_title("Case 03: Temperature at t=0.5s")
    fig.tight_layout()
    save_fig(fig, "03", "case03_temperature_final.png")

    # Plot 3: three timestep snapshots
    n_steps = len(times)
    idx_early = max(0, n_steps // 10)
    idx_mid = n_steps // 2
    idx_late = n_steps - 1
    t_labels = [
        f"t={times[idx_early]:.2f}",
        f"t={times[idx_mid]:.2f}",
        f"t={times[idx_late]:.2f}",
    ]

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    fig.suptitle("Case 03: Temperature Evolution — t=0.05, 0.25, 0.50", fontsize=FONT_TITLE)
    for col, (idx, label) in enumerate(zip([idx_early, idx_mid, idx_late], t_labels)):
        T_snap = get_nod_var(ds, 1, timestep=idx)
        tcf = axes[col].tricontourf(x, y, T_snap, levels=20, cmap=CMAP_TEMP)
        plt.colorbar(tcf, ax=axes[col])
        axes[col].set_title(label)
        axes[col].set_xlabel("x")
        axes[col].set_ylabel("y")
        axes[col].set_aspect("equal")
    ds.close()
    fig.tight_layout()
    save_fig(fig, "03", "case03_temperature_snapshots.png")


def plot_case04():
    """Case 04: Manufactured Solution — numerical, exact, and error"""
    print("Case 04: Manufactured Solution")
    efile = case_path("04", "case04_manufactured_out.e")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case04")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    u_num = get_nod_var(ds, 1, timestep=-1)
    ds.close()

    u_exact = np.sin(math.pi * x) * np.sin(math.pi * y)
    u_error = np.abs(u_num - u_exact)

    vmin = min(u_num.min(), u_exact.min())
    vmax = max(u_num.max(), u_exact.max())

    # Plot 1: numerical solution
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    tcf = ax.tricontourf(x, y, u_num, levels=20, cmap=CMAP_SCALAR,
                         vmin=vmin, vmax=vmax)
    plt.colorbar(tcf, ax=ax, label="u")
    ax.set_title("Case 04: Numerical Solution")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    fig.tight_layout()
    save_fig(fig, "04", "case04_numerical.png")

    # Plot 2: exact solution
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    tcf = ax.tricontourf(x, y, u_exact, levels=20, cmap=CMAP_SCALAR,
                         vmin=vmin, vmax=vmax)
    plt.colorbar(tcf, ax=ax, label="u")
    ax.set_title("Case 04: Exact Solution sin(\u03c0x)sin(\u03c0y)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    fig.tight_layout()
    save_fig(fig, "04", "case04_exact.png")

    # Plot 3: pointwise error
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    tcf = ax.tricontourf(x, y, u_error, levels=20, cmap="magma")
    plt.colorbar(tcf, ax=ax, label="|error|")
    ax.set_title("Case 04: Point-wise Error")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    fig.tight_layout()
    save_fig(fig, "04", "case04_error.png")


def plot_case05():
    """Case 05: Varying Conductivity k(x) = 1 + x"""
    print("Case 05: Varying Conductivity")
    efile = case_path("05", "case05_varying_k_out.e")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case05")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    u = get_nod_var(ds, 1, timestep=-1)
    ds.close()

    # Plot 1: 2D contour
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    add_2d_contour(ax, x, y, u, cmap=CMAP_SCALAR, label="u")
    ax.set_title("Case 05: Solution with k(x) = 1+x")
    fig.tight_layout()
    save_fig(fig, "05", "case05_contour_2d.png")

    # Plot 2: line plot along y = 0.5 vs exact
    idx_line = nodes_near_y(y, 0.5, tol=0.02)
    if len(idx_line) < 2:
        # Fall back to a wider tolerance if the mesh doesn't have nodes near y=0.5
        idx_line = nodes_near_y(y, 0.5, tol=0.06)

    x_line = x[idx_line]
    u_line = u[idx_line]
    order = np.argsort(x_line)
    x_line = x_line[order]
    u_line = u_line[order]

    # Exact: u = ln(1+x) / ln(2)
    x_exact = np.linspace(0.0, 1.0, 200)
    u_exact = np.log(1.0 + x_exact) / np.log(2.0)

    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
    ax.plot(x_line, u_line, "b-o", markersize=4, label="MOOSE (y \u22480.5)")
    ax.plot(x_exact, u_exact, "r--", linewidth=2, label="Exact: ln(1+x)/ln(2)")
    ax.set_xlabel("x")
    ax.set_ylabel("u")
    ax.set_title("Case 05: u along y=0.5 vs Exact Solution")
    ax.legend()
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    save_fig(fig, "05", "case05_line_exact.png")


def plot_case06():
    """Case 06: Two-Material Domain"""
    print("Case 06: Two-Material Domain")
    efile = case_path("06", "case06_two_materials_out.e")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case06")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    u = get_nod_var(ds, 1, timestep=-1)
    ds.close()

    # Plot 1: 2D contour
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    add_2d_contour(ax, x, y, u, cmap=CMAP_SCALAR, label="u")
    ax.axvline(x=0.5, color="white", linestyle="--", linewidth=1.5, label="Interface")
    ax.set_title("Case 06: Two-Material Domain (k=1 left, k=5 right)")
    ax.legend(loc="upper left")
    fig.tight_layout()
    save_fig(fig, "06", "case06_contour_2d.png")

    # Plot 2: line along y = 0.5 showing kink
    idx_line = nodes_near_y(y, 0.5, tol=0.02)
    if len(idx_line) < 2:
        idx_line = nodes_near_y(y, 0.5, tol=0.06)

    x_line = x[idx_line]
    u_line = u[idx_line]
    order = np.argsort(x_line)
    x_line = x_line[order]
    u_line = u_line[order]

    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
    ax.plot(x_line, u_line, "b-o", markersize=4, label="u(x, y\u22480.5)")
    ax.axvline(x=0.5, color="red", linestyle="--", linewidth=1.5, label="k interface")
    ax.set_xlabel("x")
    ax.set_ylabel("u")
    ax.set_title("Case 06: u along y=0.5 — Slope Change at Interface")
    ax.legend()
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    save_fig(fig, "06", "case06_line_interface.png")


def plot_case07():
    """Case 07: Nonlinear Diffusion k(T) = 1 + T"""
    print("Case 07: Nonlinear Diffusion")
    efile = case_path("07", "case07_nonlinear_diffusion_out.e")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case07")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    T = get_nod_var(ds, 1, timestep=-1)
    ds.close()

    # Plot 1: 2D contour
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    add_2d_contour(ax, x, y, T, cmap=CMAP_TEMP, label="T")
    ax.set_title("Case 07: Nonlinear Diffusion k(T) = 1+T")
    fig.tight_layout()
    save_fig(fig, "07", "case07_contour_2d.png")

    # Plot 2: 3D surface
    fig = plt.figure(figsize=FIGSIZE_SQUARE)
    ax3 = fig.add_subplot(111, projection="3d")
    ax3.plot_trisurf(x, y, T, cmap=CMAP_TEMP, linewidth=0, antialiased=False)
    ax3.set_xlabel("x")
    ax3.set_ylabel("y")
    ax3.set_zlabel("T")
    ax3.set_title("Case 07: 3D Temperature Surface")
    fig.tight_layout()
    save_fig(fig, "07", "case07_surface_3d.png")


def plot_case08():
    """Case 08: Advection-Diffusion — Gaussian blob"""
    print("Case 08: Advection-Diffusion")
    efile = case_path("08", "case08_advection_diffusion_out.e")
    cfile = case_path("08", "case08_advection_diffusion_out.csv")
    if not file_exists(efile, cfile):
        print("    SKIP: output file(s) not found")
        skipped_cases.append("case08")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    times = get_times(ds)

    # Plot 1: 2x2 grid of snapshots at t~0, 0.5, 1.0, 2.0
    target_times = [0.0, 0.5, 1.0, 2.0]
    snap_indices = [find_closest_timestep(times, t) for t in target_times]

    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    fig.suptitle("Case 08: Advecting Gaussian Blob", fontsize=FONT_TITLE)
    for panel, (idx, t_tgt) in enumerate(zip(snap_indices, target_times)):
        row, col = divmod(panel, 2)
        c_snap = get_nod_var(ds, 1, timestep=idx)
        tcf = axes[row, col].tricontourf(x, y, c_snap, levels=20, cmap=CMAP_SCALAR)
        plt.colorbar(tcf, ax=axes[row, col], label="c")
        axes[row, col].set_title(f"t \u2248 {t_tgt:.1f}  (step {idx})")
        axes[row, col].set_xlabel("x")
        axes[row, col].set_ylabel("y")
        axes[row, col].set_aspect("equal")
    ds.close()
    fig.tight_layout()
    save_fig(fig, "08", "case08_blob_snapshots.png")

    # Plot 2: total concentration vs time (conservation check)
    data = read_csv(cfile)
    t_csv = data["time"]
    total_c = data["total_c"]

    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
    ax.plot(t_csv, total_c, "b-o", markersize=3)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("total_c")
    ax.set_title("Case 08: Total Concentration vs Time")
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    save_fig(fig, "08", "case08_total_concentration.png")


def plot_case09():
    """Case 09: Coupled System (u, v)"""
    print("Case 09: Coupled System")
    efile = case_path("09", "case09_coupled_system_out.e")
    cfile = case_path("09", "case09_coupled_system_out.csv")
    if not file_exists(efile, cfile):
        print("    SKIP: output file(s) not found")
        skipped_cases.append("case09")
        return

    # Plot 1: CSV — avg_u and avg_v vs time
    data = read_csv(cfile)
    t_csv = data["time"]
    avg_u = data["avg_u"]
    avg_v = data["avg_v"]

    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
    ax.plot(t_csv, avg_u, "b-o", markersize=4, label="avg_u")
    ax.plot(t_csv, avg_v, "r-s", markersize=4, label="avg_v")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Average value")
    ax.set_title("Case 09: Coupled Variables — Average Values")
    ax.legend()
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    save_fig(fig, "09", "case09_coupled_averages.png")

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    u_final = get_nod_var(ds, 1, timestep=-1)
    v_final = get_nod_var(ds, 2, timestep=-1)
    ds.close()

    # Plot 2: field u at final time
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    add_2d_contour(ax, x, y, u_final, cmap=CMAP_SCALAR, label="u")
    ax.set_title("Case 09: Field u at t=2.0")
    fig.tight_layout()
    save_fig(fig, "09", "case09_u_final.png")

    # Plot 3: field v at final time
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    add_2d_contour(ax, x, y, v_final, cmap=CMAP_SCALAR, label="v")
    ax.set_title("Case 09: Field v at t=2.0")
    fig.tight_layout()
    save_fig(fig, "09", "case09_v_final.png")


def plot_case10():
    """Case 10: Adaptive Mesh Refinement"""
    print("Case 10: Adaptive Mesh Refinement")
    # MOOSE writes a series file for AMR; use the last one written (.e-s002)
    efile_series = case_path("10", "case10_adaptive_refinement_out.e-s002")
    efile_base = case_path("10", "case10_adaptive_refinement_out.e")

    # Prefer the series file (last AMR step) but fall back to the base file
    efile = efile_series if file_exists(efile_series) else efile_base
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case10")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    u = get_nod_var(ds, 1, timestep=-1)
    ds.close()

    # Plot 1: 2D contour on refined mesh
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    add_2d_contour(ax, x, y, u, cmap=CMAP_SCALAR, label="u")
    ax.set_title("Case 10: AMR Solution")
    fig.tight_layout()
    save_fig(fig, "10", "case10_amr_solution.png")

    # Plot 2: scatter of node positions colored by u — reveals refinement zones
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    sc = ax.scatter(x, y, c=u, s=10, cmap=CMAP_SCALAR)
    plt.colorbar(sc, ax=ax, label="u")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.set_title("Case 10: Node Distribution After AMR")
    fig.tight_layout()
    save_fig(fig, "10", "case10_node_distribution.png")


def plot_case11():
    """Case 11: Adaptive Time Stepping"""
    print("Case 11: Adaptive Timestepping")
    efile = case_path("11", "case11_adaptive_dt_out.e")
    cfile = case_path("11", "case11_adaptive_dt_out.csv")
    if not file_exists(cfile, efile):
        print("    SKIP: output file(s) not found")
        skipped_cases.append("case11")
        return

    data = read_csv(cfile)
    t_csv = data["time"]
    max_T = data["max_T"]
    avg_T = data["avg_T"]
    dt_pp = data["dt_pp"]

    # Plot 1: max_T and avg_T vs time
    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
    ax.plot(t_csv, max_T, "r-o", markersize=4, label="max_T")
    ax.plot(t_csv, avg_T, "b-s", markersize=4, label="avg_T")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Temperature")
    ax.set_title("Case 11: Temperature with Adaptive Timestepping")
    ax.legend()
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    save_fig(fig, "11", "case11_temperature_history.png")

    # Plot 2: dt vs time (log scale)
    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
    # Guard against zero dt values before log-scale plot
    dt_safe = [max(v, 1e-12) for v in dt_pp]
    ax.semilogy(t_csv, dt_safe, "k-x", markersize=5)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("dt (log scale)")
    ax.set_title("Case 11: Adaptive Timestep Size")
    ax.grid(True, which="both", alpha=0.4)
    fig.tight_layout()
    save_fig(fig, "11", "case11_timestep_size.png")

    # Plot 3: 2D contour at final time
    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    T_final = get_nod_var(ds, 1, timestep=-1)
    ds.close()

    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    add_2d_contour(ax, x, y, T_final, cmap=CMAP_TEMP, label="T")
    ax.set_title("Case 11: Temperature at Steady State")
    fig.tight_layout()
    save_fig(fig, "11", "case11_temperature_final.png")


def plot_case12():
    """Case 12: MultiApp Coupling"""
    print("Case 12: MultiApp Coupling")
    efile_parent = case_path("12", "case12_parent_out.e")
    efile_sub = case_path("12", "case12_parent_out_thermal_sub0.e")
    if not file_exists(efile_parent):
        print("    SKIP: parent output file not found")
        skipped_cases.append("case12")
        return

    # Plot 1: parent temperature field T
    ds_p = open_exodus(efile_parent)
    x_p, y_p = get_coords_2d(ds_p)
    T_parent = get_nod_var(ds_p, 1, timestep=-1)
    ds_p.close()

    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    add_2d_contour(ax, x_p, y_p, T_parent, cmap=CMAP_TEMP, label="T")
    ax.set_title("Case 12: Parent Temperature T")
    fig.tight_layout()
    save_fig(fig, "12", "case12_parent_temperature.png")

    # Plot 2: sub-app phi field (vals_nod_var2 in the sub file)
    if not file_exists(efile_sub):
        print("    NOTE: sub-app output not found — skipping sub plot")
        return

    ds_s = open_exodus(efile_sub)
    x_s, y_s = get_coords_2d(ds_s)
    # The sub file has: vals_nod_var1 = T_from_parent, vals_nod_var2 = phi
    n_vars = ds_s.variables["name_nod_var"].shape[0]
    phi_index = 2 if n_vars >= 2 else 1
    phi = get_nod_var(ds_s, phi_index, timestep=-1)
    ds_s.close()

    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    add_2d_contour(ax, x_s, y_s, phi, cmap=CMAP_SCALAR, label="\u03c6")
    ax.set_title("Case 12: Sub-App Field \u03c6 (driven by T)")
    fig.tight_layout()
    save_fig(fig, "12", "case12_sub_phi.png")


def plot_case13():
    """Case 13: Comprehensive Postprocessor Analysis"""
    print("Case 13: Postprocessor Analysis")
    cfile = case_path("13", "case13_postprocessors_out.csv")
    efile = case_path("13", "case13_postprocessors_out.e")
    if not file_exists(cfile):
        print("    SKIP: CSV output file not found")
        skipped_cases.append("case13")
        return

    data = read_csv(cfile)
    time = data["time"]
    max_temp = data["max_temp"]
    avg_temp = data["avg_temp"]
    total_energy = data["total_energy"]
    T_L2_norm = data["T_L2_norm"]
    dt_size = data["dt_size"]

    # Analytical steady-state average temperature approximation (1-term Fourier)
    k = 2.0
    Q = 5.0
    T_ss_approx = Q * 8.0 / (k * math.pi ** 4)

    # Plot 1: 4-panel comprehensive figure
    fig, axes = plt.subplots(2, 2, figsize=FIGSIZE_LARGE)
    fig.suptitle("Case 13: Comprehensive Postprocessor Analysis", fontsize=FONT_TITLE)

    # Top-left: temperature histories
    axes[0, 0].plot(time, max_temp, "r-o", markersize=3, label="max_temp")
    axes[0, 0].plot(time, avg_temp, "b-s", markersize=3, label="avg_temp")
    axes[0, 0].axhline(
        T_ss_approx, color="gray", linestyle="--", linewidth=1.2,
        label=f"Steady-state avg ({T_ss_approx:.4f})"
    )
    axes[0, 0].set_xlabel("Time (s)")
    axes[0, 0].set_ylabel("Temperature")
    axes[0, 0].set_title("max_temp and avg_temp vs Time")
    axes[0, 0].legend(fontsize=9)
    axes[0, 0].grid(True, alpha=0.4)

    # Top-right: total energy
    axes[0, 1].plot(time, total_energy, "g-^", markersize=3)
    axes[0, 1].set_xlabel("Time (s)")
    axes[0, 1].set_ylabel("\u222bT dV")
    axes[0, 1].set_title("total_energy vs Time")
    axes[0, 1].grid(True, alpha=0.4)

    # Bottom-left: L2 norm
    axes[1, 0].plot(time, T_L2_norm, "m-D", markersize=3)
    axes[1, 0].set_xlabel("Time (s)")
    axes[1, 0].set_ylabel("||T||_L2")
    axes[1, 0].set_title("T_L2_norm vs Time")
    axes[1, 0].grid(True, alpha=0.4)

    # Bottom-right: adaptive dt (log scale)
    dt_safe = [max(v, 1e-12) for v in dt_size]
    axes[1, 1].semilogy(time, dt_safe, "k-x", markersize=3)
    axes[1, 1].set_xlabel("Time (s)")
    axes[1, 1].set_ylabel("dt (log scale)")
    axes[1, 1].set_title("dt_size vs Time")
    axes[1, 1].grid(True, which="both", alpha=0.4)

    fig.tight_layout()
    save_fig(fig, "13", "case13_postprocessors.png")

    # Bonus: 2D temperature field at final time (if Exodus present)
    if file_exists(efile):
        ds = open_exodus(efile)
        x, y = get_coords_2d(ds)
        T_final = get_nod_var(ds, 1, timestep=-1)
        ds.close()

        fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
        add_2d_contour(ax, x, y, T_final, cmap=CMAP_TEMP, label="T")
        ax.set_title("Case 13: Temperature at Final Time")
        fig.tight_layout()
        save_fig(fig, "13", "case13_temperature_final.png")


# ---------------------------------------------------------------------------
# Cases 14-21: Advanced Multi-Physics
# ---------------------------------------------------------------------------

def plot_case14():
    """Case 14: Thermoelasticity — Temperature, displacement, and von Mises stress"""
    print("Case 14: Thermoelasticity")
    efile = case_path("14", "case14_thermoelasticity_out.e")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case14")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    names = get_nod_var_names(ds)
    elem_names = get_elem_var_names(ds)
    cx, cy = get_elem_centroids_2d(ds)

    # Nodal variables
    T_idx = names.index("T") + 1
    disp_x_idx = names.index("disp_x") + 1
    disp_y_idx = names.index("disp_y") + 1

    T = get_nod_var(ds, T_idx, timestep=-1)
    disp_x = get_nod_var(ds, disp_x_idx, timestep=-1)
    disp_y = get_nod_var(ds, disp_y_idx, timestep=-1)

    # Element variable (von Mises stress lives on elements, not nodes)
    vm_idx = elem_names.index("vonmises_stress") + 1
    vonmises = get_elem_var(ds, vm_idx, timestep=-1)
    ds.close()

    # Plot 1: 2x2 panel — T, disp_x, disp_y, vonmises_stress
    fig, axes = plt.subplots(2, 2, figsize=FIGSIZE_LARGE)
    fig.suptitle("Case 14: Thermoelasticity — Heated Plate", fontsize=FONT_TITLE)

    # Top row: nodal variables (T, disp_x) — use tricontourf on nodes
    for col_idx, (vals, cmap, label) in enumerate([
        (T, CMAP_TEMP, "Temperature (K)"),
        (disp_x, CMAP_SCALAR, "disp_x (m)"),
    ]):
        tcf = axes[0, col_idx].tricontourf(x, y, vals, levels=20, cmap=cmap)
        plt.colorbar(tcf, ax=axes[0, col_idx], label=label)
        axes[0, col_idx].set_xlabel("x")
        axes[0, col_idx].set_ylabel("y")
        axes[0, col_idx].set_aspect("equal")
        axes[0, col_idx].set_title(label.split("(")[0].strip())

    # Bottom-left: disp_y (nodal)
    tcf = axes[1, 0].tricontourf(x, y, disp_y, levels=20, cmap=CMAP_SCALAR)
    plt.colorbar(tcf, ax=axes[1, 0], label="disp_y (m)")
    axes[1, 0].set_xlabel("x")
    axes[1, 0].set_ylabel("y")
    axes[1, 0].set_aspect("equal")
    axes[1, 0].set_title("disp_y")

    # Bottom-right: von Mises stress (element variable — use scatter on centroids)
    sc = axes[1, 1].scatter(cx, cy, c=vonmises, s=8, cmap="inferno", edgecolors="none")
    plt.colorbar(sc, ax=axes[1, 1], label="von Mises (Pa)")
    axes[1, 1].set_xlabel("x")
    axes[1, 1].set_ylabel("y")
    axes[1, 1].set_aspect("equal")
    axes[1, 1].set_title("von Mises Stress")

    fig.tight_layout()
    save_fig(fig, "14", "case14_thermoelasticity.png")


def plot_case15():
    """Case 15: Lid-Driven Cavity — velocity and pressure fields"""
    print("Case 15: Lid-Driven Cavity")
    efile = case_path("15", "case15_lid_driven_cavity_out.e")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case15")
        return

    ds = open_exodus(efile)
    cx, cy = get_elem_centroids_2d(ds)
    names = get_elem_var_names(ds)

    # FV variables are element-centered
    vel_x_name = "vel_x" if "vel_x" in names else "superficial_vel_x"
    vel_y_name = "vel_y" if "vel_y" in names else "superficial_vel_y"
    p_name = "pressure"

    vel_x_idx = names.index(vel_x_name) + 1
    vel_y_idx = names.index(vel_y_name) + 1
    p_idx = names.index(p_name) + 1

    vx = get_elem_var(ds, vel_x_idx, timestep=-1)
    vy = get_elem_var(ds, vel_y_idx, timestep=-1)
    p = get_elem_var(ds, p_idx, timestep=-1)
    ds.close()

    vel_mag = np.sqrt(vx ** 2 + vy ** 2)

    # Plot 1: velocity magnitude contour (element centroids)
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    tcf = ax.tricontourf(cx, cy, vel_mag, levels=20, cmap="plasma")
    plt.colorbar(tcf, ax=ax, label="|velocity|")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.set_title("Case 15: Lid-Driven Cavity — Velocity Magnitude (Re=100)")
    fig.tight_layout()
    save_fig(fig, "15", "case15_velocity_magnitude.png")

    # Plot 2: pressure field
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    tcf = ax.tricontourf(cx, cy, p, levels=20, cmap="RdBu_r")
    plt.colorbar(tcf, ax=ax, label="Pressure")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.set_title("Case 15: Pressure Field")
    fig.tight_layout()
    save_fig(fig, "15", "case15_pressure.png")

    # Plot 3: quiver plot
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    n = len(cx)
    step = max(1, n // 400)
    ax.quiver(cx[::step], cy[::step], vx[::step], vy[::step],
              vel_mag[::step], cmap="plasma", scale=15, width=0.003)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.set_title("Case 15: Velocity Vectors")
    fig.tight_layout()
    save_fig(fig, "15", "case15_velocity_vectors.png")


def plot_case16():
    """Case 16: Natural Convection — temperature and velocity in heated cavity"""
    print("Case 16: Natural Convection")
    efile = case_path("16", "case16_natural_convection_out.e")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case16")
        return

    ds = open_exodus(efile)
    cx, cy = get_elem_centroids_2d(ds)
    names = get_elem_var_names(ds)

    T_name = "T_fluid" if "T_fluid" in names else "temperature"
    vel_x_name = "vel_x" if "vel_x" in names else "superficial_vel_x"
    vel_y_name = "vel_y" if "vel_y" in names else "superficial_vel_y"

    T_idx = names.index(T_name) + 1
    vx_idx = names.index(vel_x_name) + 1
    vy_idx = names.index(vel_y_name) + 1

    T = get_elem_var(ds, T_idx, timestep=-1)
    vx = get_elem_var(ds, vx_idx, timestep=-1)
    vy = get_elem_var(ds, vy_idx, timestep=-1)
    ds.close()

    vel_mag = np.sqrt(vx ** 2 + vy ** 2)

    # Plot 1: temperature field
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    tcf = ax.tricontourf(cx, cy, T, levels=20, cmap=CMAP_TEMP)
    plt.colorbar(tcf, ax=ax, label="Temperature")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.set_title("Case 16: Natural Convection — Temperature (Ra=10$^4$)")
    fig.tight_layout()
    save_fig(fig, "16", "case16_temperature.png")

    # Plot 2: velocity magnitude
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    tcf = ax.tricontourf(cx, cy, vel_mag, levels=20, cmap="plasma")
    plt.colorbar(tcf, ax=ax, label="|velocity|")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.set_title("Case 16: Velocity Magnitude")
    fig.tight_layout()
    save_fig(fig, "16", "case16_velocity_magnitude.png")


def plot_case17():
    """Case 17: Joule Heating — voltage and temperature evolution"""
    print("Case 17: Joule Heating")
    efile = case_path("17", "case17_joule_heating_out.e")
    cfile = case_path("17", "case17_joule_heating_out.csv")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case17")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    names = get_nod_var_names(ds)
    times = get_times(ds)

    V_idx = names.index("V") + 1
    T_idx = names.index("T") + 1

    V = get_nod_var(ds, V_idx, timestep=-1)
    T_final = get_nod_var(ds, T_idx, timestep=-1)

    # Plot 1: side-by-side V and T at final time
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Case 17: Joule Heating — Final State", fontsize=FONT_TITLE)

    tcf1 = axes[0].tricontourf(x, y, V, levels=20, cmap="RdYlBu_r")
    plt.colorbar(tcf1, ax=axes[0], label="V (Volts)")
    axes[0].set_title("Electric Potential V")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    axes[0].set_aspect("equal")

    tcf2 = axes[1].tricontourf(x, y, T_final, levels=20, cmap=CMAP_TEMP)
    plt.colorbar(tcf2, ax=axes[1], label="T (K)")
    axes[1].set_title("Temperature T")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("y")
    axes[1].set_aspect("equal")

    ds.close()
    fig.tight_layout()
    save_fig(fig, "17", "case17_joule_heating.png")

    # Plot 2: temperature history from CSV
    if file_exists(cfile):
        data = read_csv(cfile)
        t_csv = data["time"]
        max_T = data.get("max_T", data.get("max_temperature", []))
        avg_T = data.get("avg_T", data.get("avg_temperature", []))

        if max_T and avg_T:
            fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
            ax.plot(t_csv, max_T, "r-o", markersize=3, label="max_T")
            ax.plot(t_csv, avg_T, "b-s", markersize=3, label="avg_T")
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Temperature (K)")
            ax.set_title("Case 17: Temperature Rise from Joule Heating")
            ax.legend()
            ax.grid(True, alpha=0.4)
            fig.tight_layout()
            save_fig(fig, "17", "case17_temperature_history.png")


def plot_case18():
    """Case 18: Cahn-Hilliard — phase separation snapshots"""
    print("Case 18: Cahn-Hilliard Spinodal Decomposition")
    efile = case_path("18", "case18_cahn_hilliard_out.e")
    cfile = case_path("18", "case18_cahn_hilliard_out.csv")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case18")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    names = get_nod_var_names(ds)
    times = get_times(ds)

    c_idx = names.index("c") + 1

    # Plot 1: 2x2 snapshots at different times
    n_steps = len(times)
    target_frac = [0.0, 0.15, 0.5, 1.0]
    snap_indices = [max(0, min(n_steps - 1, int(f * (n_steps - 1)))) for f in target_frac]

    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    fig.suptitle("Case 18: Cahn-Hilliard Spinodal Decomposition", fontsize=FONT_TITLE)
    for panel, idx in enumerate(snap_indices):
        row, col = divmod(panel, 2)
        c_snap = get_nod_var(ds, c_idx, timestep=idx)
        tcf = axes[row, col].tricontourf(x, y, c_snap, levels=20, cmap="PiYG")
        plt.colorbar(tcf, ax=axes[row, col], label="c")
        axes[row, col].set_title(f"t = {times[idx]:.1f}")
        axes[row, col].set_xlabel("x")
        axes[row, col].set_ylabel("y")
        axes[row, col].set_aspect("equal")
    ds.close()
    fig.tight_layout()
    save_fig(fig, "18", "case18_phase_separation.png")

    # Plot 2: free energy and avg_c from CSV
    if file_exists(cfile):
        data = read_csv(cfile)
        t_csv = data["time"]
        avg_c = data.get("avg_c", [])
        bulk_E = data.get("bulk_energy", data.get("total_free_energy", []))

        if avg_c and bulk_E:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=FIGSIZE_WIDE)
            fig.suptitle("Case 18: Conservation and Energy", fontsize=FONT_TITLE)

            ax1.plot(t_csv, avg_c, "b-", linewidth=1.5)
            ax1.set_xlabel("Time")
            ax1.set_ylabel("Average c")
            ax1.set_title("Mass Conservation")
            ax1.grid(True, alpha=0.4)

            ax2.plot(t_csv, bulk_E, "r-", linewidth=1.5)
            ax2.set_xlabel("Time")
            ax2.set_ylabel("Bulk Free Energy")
            ax2.set_title("Free Energy Decay")
            ax2.grid(True, alpha=0.4)

            fig.tight_layout()
            save_fig(fig, "18", "case18_energy_conservation.png")


def plot_case19():
    """Case 19: Porous Flow — pressure and temperature fields"""
    print("Case 19: Porous Flow with Heat Transport")
    efile = case_path("19", "case19_porous_flow_out.e")
    cfile = case_path("19", "case19_porous_flow_out.csv")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case19")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    names = get_nod_var_names(ds)
    times = get_times(ds)

    p_idx = names.index("porepressure") + 1
    T_idx = names.index("temperature") + 1

    p_final = get_nod_var(ds, p_idx, timestep=-1)
    T_final = get_nod_var(ds, T_idx, timestep=-1)

    # Plot 1: side-by-side pressure and temperature
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    fig.suptitle("Case 19: Porous Flow — Final State", fontsize=FONT_TITLE)

    tcf1 = axes[0].tricontourf(x, y, p_final, levels=20, cmap="Blues")
    plt.colorbar(tcf1, ax=axes[0], label="Pressure (Pa)")
    axes[0].set_title("Porepressure")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    axes[0].set_aspect("equal")

    tcf2 = axes[1].tricontourf(x, y, T_final, levels=20, cmap=CMAP_TEMP)
    plt.colorbar(tcf2, ax=axes[1], label="T (K)")
    axes[1].set_title("Temperature")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("y")
    axes[1].set_aspect("equal")

    # Plot 2: temperature snapshots (thermal plume evolution)
    n_steps = len(times)
    if n_steps >= 4:
        target_frac = [0.0, 0.33, 0.66, 1.0]
        snap_indices = [max(0, min(n_steps - 1, int(f * (n_steps - 1)))) for f in target_frac]

        fig2, axes2 = plt.subplots(1, 4, figsize=(16, 3.5))
        fig2.suptitle("Case 19: Thermal Plume Evolution", fontsize=FONT_TITLE)
        for col, idx in enumerate(snap_indices):
            T_snap = get_nod_var(ds, T_idx, timestep=idx)
            tcf = axes2[col].tricontourf(x, y, T_snap, levels=20, cmap=CMAP_TEMP,
                                          vmin=295, vmax=355)
            plt.colorbar(tcf, ax=axes2[col])
            axes2[col].set_title(f"t = {times[idx]:.0f} s")
            axes2[col].set_xlabel("x")
            axes2[col].set_ylabel("y")
            axes2[col].set_aspect("equal")
        fig2.tight_layout()
        save_fig(fig2, "19", "case19_thermal_plume.png")

    ds.close()
    fig.tight_layout()
    save_fig(fig, "19", "case19_porous_flow.png")

    # Plot 3: CSV time history
    if file_exists(cfile):
        data = read_csv(cfile)
        t_csv = data["time"]
        max_T = data.get("max_T", [])
        avg_T = data.get("avg_T", [])

        if max_T and avg_T:
            fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
            ax.plot(t_csv, max_T, "r-o", markersize=3, label="max_T")
            ax.plot(t_csv, avg_T, "b-s", markersize=3, label="avg_T")
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Temperature (K)")
            ax.set_title("Case 19: Temperature History")
            ax.legend()
            ax.grid(True, alpha=0.4)
            fig.tight_layout()
            save_fig(fig, "19", "case19_temperature_history.png")


def plot_case20():
    """Case 20: Elastic Wave — displacement and stress snapshots"""
    print("Case 20: Elastic Wave Propagation")
    # The named [exodus] sub-block produces a file without "_out" in the name:
    #   case20_elastic_wave_exodus.e  (not case20_elastic_wave_out.e)
    efile = case_path("20", "case20_elastic_wave_exodus.e")
    cfile = case_path("20", "case20_elastic_wave_out.csv")

    if not file_exists(efile):
        # Fallback to standard naming in case output config changes
        alt = case_path("20", "case20_elastic_wave_out.e")
        if file_exists(alt):
            efile = alt
        else:
            print("    SKIP: output file not found")
            skipped_cases.append("case20")
            return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    names = get_nod_var_names(ds)
    times = get_times(ds)

    dx_idx = names.index("disp_x") + 1

    # Plot 1: displacement snapshots at different times
    n_steps = len(times)
    if n_steps < 4:
        snap_indices = list(range(n_steps))
    else:
        target_frac = [0.1, 0.3, 0.6, 0.9]
        snap_indices = [max(0, min(n_steps - 1, int(f * (n_steps - 1)))) for f in target_frac]

    fig, axes = plt.subplots(len(snap_indices), 1, figsize=(12, 3 * len(snap_indices)))
    if len(snap_indices) == 1:
        axes = [axes]
    fig.suptitle("Case 20: Elastic Wave — disp_x Along Bar", fontsize=FONT_TITLE)
    for row, idx in enumerate(snap_indices):
        dx_snap = get_nod_var(ds, dx_idx, timestep=idx)
        tcf = axes[row].tricontourf(x, y, dx_snap, levels=20, cmap="RdBu_r")
        plt.colorbar(tcf, ax=axes[row], label="disp_x (m)")
        axes[row].set_title(f"t = {times[idx]:.5f} s")
        axes[row].set_xlabel("x (m)")
        axes[row].set_ylabel("y (m)")
    ds.close()
    fig.tight_layout()
    save_fig(fig, "20", "case20_wave_snapshots.png")

    # Plot 2: tip displacement vs time from CSV
    if file_exists(cfile):
        data = read_csv(cfile)
        t_csv = data["time"]
        tip_dx = data.get("disp_x_right", data.get("tip_disp_x", []))

        if tip_dx:
            fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
            ax.plot(t_csv, tip_dx, "b-", linewidth=1)
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("disp_x at right end (m)")
            ax.set_title("Case 20: Tip Displacement — Wave Arrival and Reflection")
            ax.grid(True, alpha=0.4)
            fig.tight_layout()
            save_fig(fig, "20", "case20_tip_displacement.png")


def plot_case21():
    """Case 21: Bimetallic Strip — deformation and stress"""
    print("Case 21: Bimetallic Strip")
    efile = case_path("21", "case21_bimetallic_strip_out.e")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case21")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    names = get_nod_var_names(ds)
    elem_names = get_elem_var_names(ds)

    dx_idx = names.index("disp_x") + 1
    dy_idx = names.index("disp_y") + 1

    disp_x = get_nod_var(ds, dx_idx, timestep=-1)
    disp_y = get_nod_var(ds, dy_idx, timestep=-1)

    # von Mises stress is an element variable — read from element centroids
    vm_idx = elem_names.index("vonmises_stress") + 1
    # Case 21 has two element blocks (steel, aluminum); read both and concatenate
    vm_parts = []
    cx_parts = []
    cy_parts = []
    for blk in range(1, 3):
        key = f"vals_elem_var{vm_idx}eb{blk}"
        if key in ds.variables:
            vm_parts.append(get_elem_var(ds, vm_idx, block=blk, timestep=-1))
            blk_cx, blk_cy = get_elem_centroids_2d(ds, block=blk)
            cx_parts.append(blk_cx)
            cy_parts.append(blk_cy)
    vonmises = np.concatenate(vm_parts)
    cx = np.concatenate(cx_parts)
    cy = np.concatenate(cy_parts)
    ds.close()

    # Plot 1: deformed shape colored by disp_y
    # Scale displacements for visibility
    disp_scale = 5.0
    x_def = x + disp_x * disp_scale
    y_def = y + disp_y * disp_scale

    fig, axes = plt.subplots(2, 1, figsize=(12, 6))
    fig.suptitle("Case 21: Bimetallic Strip — Thermal Bending", fontsize=FONT_TITLE)

    # Top: original + deformed shape
    axes[0].scatter(x, y, c="lightgray", s=2, label="Original")
    sc = axes[0].scatter(x_def, y_def, c=disp_y, s=3, cmap="RdBu_r")
    plt.colorbar(sc, ax=axes[0], label="disp_y (m)")
    axes[0].set_xlabel("x (m)")
    axes[0].set_ylabel("y (m)")
    axes[0].set_title(f"Deformed Shape ({disp_scale}x magnification)")
    axes[0].set_aspect("equal")
    axes[0].legend(loc="upper right", markerscale=5)

    # Bottom: von Mises stress on element centroids
    sc2 = axes[1].scatter(cx, cy, c=vonmises, s=8, cmap="inferno", edgecolors="none")
    plt.colorbar(sc2, ax=axes[1], label="von Mises (Pa)")
    axes[1].axhline(y=0.5, color="white", linestyle="--", linewidth=1, alpha=0.7,
                     label="Material interface")
    axes[1].set_xlabel("x (m)")
    axes[1].set_ylabel("y (m)")
    axes[1].set_title("von Mises Stress Distribution")
    axes[1].set_aspect("equal")
    axes[1].legend(loc="upper right")

    fig.tight_layout()
    save_fig(fig, "21", "case21_bimetallic_strip.png")


def plot_case22():
    """Case 22: Charge Relaxation — exponential decay of free charge"""
    print("Case 22: Charge Relaxation in an Ohmic Medium")
    efile = case_path("22", "case22_charge_relaxation_out.e")
    cfile = case_path("22", "case22_charge_relaxation_out.csv")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case22")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    names = get_nod_var_names(ds)
    times = get_times(ds)

    rho_idx = names.index("rho_e") + 1
    phi_idx = names.index("phi") + 1

    # Plot 1: rho_e snapshots at 4 times
    n_steps = len(times)
    target_frac = [0.0, 0.1, 0.3, 1.0]
    snap_indices = [max(0, min(n_steps - 1, int(f * (n_steps - 1)))) for f in target_frac]

    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    fig.suptitle("Case 22: Charge Relaxation — Charge Density", fontsize=FONT_TITLE)
    for panel, idx in enumerate(snap_indices):
        row, col = divmod(panel, 2)
        rho = get_nod_var(ds, rho_idx, timestep=idx)
        tcf = axes[row, col].tricontourf(x, y, rho, levels=20, cmap="YlOrRd")
        plt.colorbar(tcf, ax=axes[row, col], label=r"$\rho_e$")
        axes[row, col].set_title(f"t = {times[idx]:.3f}")
        axes[row, col].set_xlabel("x")
        axes[row, col].set_ylabel("y")
        axes[row, col].set_aspect("equal")
    ds.close()
    fig.tight_layout()
    save_fig(fig, "22", "case22_charge_relaxation.png")

    # Plot 2: time history from CSV
    if file_exists(cfile):
        data = read_csv(cfile)
        t_csv = data["time"]
        max_rho = data.get("max_rho_e", data.get("max_rho", []))
        avg_rho = data.get("avg_rho_e", data.get("avg_rho", []))

        if max_rho and avg_rho:
            fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
            ax.plot(t_csv, max_rho, "r-o", markersize=3, label="max rho_e")
            ax.plot(t_csv, avg_rho, "b-s", markersize=3, label="avg rho_e")
            # Analytical: exponential decay with tau = eps/sigma = 0.1
            t_an = np.linspace(0, max(t_csv), 200)
            ax.plot(t_an, np.exp(-t_an / 0.1), "k--", linewidth=1,
                    label=r"Analytical: $e^{-t/\tau}$, $\tau$=0.1")
            ax.set_xlabel("Time (s)")
            ax.set_ylabel(r"Charge density $\rho_e$")
            ax.set_title("Case 22: Charge Decay — Exponential Relaxation")
            ax.legend()
            ax.grid(True, alpha=0.4)
            fig.tight_layout()
            save_fig(fig, "22", "case22_charge_history.png")


def plot_case23():
    """Case 23: Magnetic Diffusion — erfc profile penetration"""
    print("Case 23: Magnetic Diffusion into a Conducting Slab")
    efile = case_path("23", "case23_magnetic_diffusion_out.e")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case23")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    names = get_nod_var_names(ds)
    times = get_times(ds)

    B_idx = names.index("B") + 1

    # Plot 1: B(x) profiles at several times (quasi-1D, pick centerline)
    y_mid = (y.max() + y.min()) / 2.0
    center = nodes_near_y(y, y_mid, tol=0.04)
    x_center = x[center]
    order = np.argsort(x_center)
    x_sorted = x_center[order]

    # Select ~5 snapshots across the simulation
    n_steps = len(times)
    snap_fracs = [0.05, 0.15, 0.3, 0.6, 1.0]
    snap_idx = [max(0, min(n_steps - 1, int(f * (n_steps - 1)))) for f in snap_fracs]

    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
    from scipy.special import erfc as _erfc
    D_m = 0.01
    for idx in snap_idx:
        B_snap = get_nod_var(ds, B_idx, timestep=idx)
        B_line = B_snap[center][order]
        t_val = times[idx]
        ax.plot(x_sorted, B_line, "-o", markersize=2, label=f"MOOSE t={t_val:.1f}")
        # Analytical erfc profile
        if t_val > 0:
            x_an = np.linspace(0, x_sorted.max(), 200)
            B_an = _erfc(x_an / (2.0 * np.sqrt(D_m * t_val)))
            ax.plot(x_an, B_an, "k--", linewidth=0.8, alpha=0.5)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("B (T)")
    ax.set_title("Case 23: Magnetic Field Diffusion Profiles")
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.4)
    ds.close()
    fig.tight_layout()
    save_fig(fig, "23", "case23_magnetic_diffusion.png")


def plot_case24():
    """Case 24: Drift-Diffusion — charge front propagation"""
    print("Case 24: Charge Drift-Diffusion Between Parallel Plates")
    efile = case_path("24", "case24_drift_diffusion_out.e")
    cfile = case_path("24", "case24_drift_diffusion_out.csv")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case24")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    names = get_nod_var_names(ds)
    times = get_times(ds)

    rho_idx = names.index("rho_e") + 1
    phi_idx = names.index("phi") + 1

    # Plot 1: rho_e(x) profiles at several times (quasi-1D centerline)
    y_mid = (y.max() + y.min()) / 2.0
    center = nodes_near_y(y, y_mid, tol=0.02)
    x_center = x[center]
    order = np.argsort(x_center)
    x_sorted = x_center[order]

    n_steps = len(times)
    snap_fracs = [0.05, 0.15, 0.3, 0.5, 0.8, 1.0]
    snap_idx = [max(0, min(n_steps - 1, int(f * (n_steps - 1)))) for f in snap_fracs]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Case 24: Charge Drift-Diffusion", fontsize=FONT_TITLE)

    for idx in snap_idx:
        rho = get_nod_var(ds, rho_idx, timestep=idx)
        axes[0].plot(x_sorted, rho[center][order], "-", linewidth=1.5,
                     label=f"t={times[idx]:.2f}")

    axes[0].set_xlabel("x (m)")
    axes[0].set_ylabel(r"$\rho_e$ (C/m$^3$)")
    axes[0].set_title("Charge Density Profiles")
    axes[0].legend(fontsize=8)
    axes[0].grid(True, alpha=0.4)

    # Right panel: phi(x) at same times
    for idx in snap_idx:
        phi = get_nod_var(ds, phi_idx, timestep=idx)
        axes[1].plot(x_sorted, phi[center][order], "-", linewidth=1.5,
                     label=f"t={times[idx]:.2f}")

    axes[1].set_xlabel("x (m)")
    axes[1].set_ylabel(r"$\phi$ (V)")
    axes[1].set_title("Electric Potential")
    axes[1].legend(fontsize=8)
    axes[1].grid(True, alpha=0.4)
    ds.close()
    fig.tight_layout()
    save_fig(fig, "24", "case24_drift_diffusion.png")

    # Plot 2: time history from CSV
    if file_exists(cfile):
        data = read_csv(cfile)
        t_csv = data["time"]
        avg_rho = data.get("avg_rho", [])
        max_rho = data.get("max_rho", [])

        if avg_rho and max_rho:
            fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
            ax.plot(t_csv, avg_rho, "b-s", markersize=3, label="avg rho_e")
            ax.plot(t_csv, max_rho, "r-o", markersize=3, label="max rho_e")
            ax.set_xlabel("Time (s)")
            ax.set_ylabel(r"$\rho_e$")
            ax.set_title("Case 24: Charge Density vs Time")
            ax.legend()
            ax.grid(True, alpha=0.4)
            fig.tight_layout()
            save_fig(fig, "24", "case24_charge_history.png")


def plot_case25():
    """Case 25: Induction Heating — skin-depth heating"""
    print("Case 25: Induction Heating (Magnetic Diffusion + Heat)")
    efile = case_path("25", "case25_induction_heating_out.e")
    cfile = case_path("25", "case25_induction_heating_out.csv")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case25")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    names = get_nod_var_names(ds)
    times = get_times(ds)

    B_idx = names.index("B") + 1
    T_idx = names.index("T") + 1

    # Plot 1: side-by-side B and T at final time (quasi-1D centerline)
    y_mid = (y.max() + y.min()) / 2.0
    center = nodes_near_y(y, y_mid, tol=0.04)
    x_center = x[center]
    order = np.argsort(x_center)
    x_sorted = x_center[order]

    B_final = get_nod_var(ds, B_idx, timestep=-1)[center][order]
    T_final = get_nod_var(ds, T_idx, timestep=-1)[center][order]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Case 25: Induction Heating — Final State", fontsize=FONT_TITLE)

    axes[0].plot(x_sorted, B_final, "b-o", markersize=2)
    axes[0].set_xlabel("x (m)")
    axes[0].set_ylabel("B (T)")
    axes[0].set_title("Magnetic Field (skin depth)")
    axes[0].grid(True, alpha=0.4)

    axes[1].plot(x_sorted, T_final, "r-o", markersize=2)
    axes[1].set_xlabel("x (m)")
    axes[1].set_ylabel("T (K)")
    axes[1].set_title("Temperature (surface heating)")
    axes[1].grid(True, alpha=0.4)
    ds.close()
    fig.tight_layout()
    save_fig(fig, "25", "case25_induction_heating.png")

    # Plot 2: temperature history from CSV
    if file_exists(cfile):
        data = read_csv(cfile)
        t_csv = data["time"]
        max_T = data.get("max_T", [])
        avg_T = data.get("avg_T", [])

        if max_T and avg_T:
            fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
            ax.plot(t_csv, max_T, "r-o", markersize=3, label="max T")
            ax.plot(t_csv, avg_T, "b-s", markersize=3, label="avg T")
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Temperature (K)")
            ax.set_title("Case 25: Temperature History — Induction Heating")
            ax.legend()
            ax.grid(True, alpha=0.4)
            fig.tight_layout()
            save_fig(fig, "25", "case25_temperature_history.png")


def plot_case26():
    """Case 26: EHD Pumping — Coulomb force driven recirculation"""
    print("Case 26: EHD Pumping")
    efile = case_path("26", "case26_ehd_pumping_out.e")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case26")
        return

    ds = open_exodus(efile)
    cx, cy = get_elem_centroids_2d(ds)
    names = get_elem_var_names(ds)

    vel_x_name = "vel_x" if "vel_x" in names else "superficial_vel_x"
    vel_y_name = "vel_y" if "vel_y" in names else "superficial_vel_y"
    p_name = "pressure"

    vx_idx = names.index(vel_x_name) + 1
    vy_idx = names.index(vel_y_name) + 1
    p_idx = names.index(p_name) + 1

    vx = get_elem_var(ds, vx_idx, timestep=-1)
    vy = get_elem_var(ds, vy_idx, timestep=-1)
    p = get_elem_var(ds, p_idx, timestep=-1)
    ds.close()

    vel_mag = np.sqrt(vx ** 2 + vy ** 2)

    # Plot 1: velocity magnitude + vectors
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Case 26: EHD Pumping — Coulomb-Force-Driven Flow", fontsize=FONT_TITLE)

    tcf = axes[0].tricontourf(cx, cy, vel_mag, levels=20, cmap="plasma")
    plt.colorbar(tcf, ax=axes[0], label="|velocity|")
    n = len(cx)
    step = max(1, n // 300)
    axes[0].quiver(cx[::step], cy[::step], vx[::step], vy[::step],
                   color="white", scale=max(vel_mag) * 8, width=0.003, alpha=0.8)
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    axes[0].set_aspect("equal")
    axes[0].set_title("Velocity Magnitude + Vectors")

    tcf2 = axes[1].tricontourf(cx, cy, p, levels=20, cmap="RdBu_r")
    plt.colorbar(tcf2, ax=axes[1], label="Pressure")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("y")
    axes[1].set_aspect("equal")
    axes[1].set_title("Pressure Field")

    fig.tight_layout()
    save_fig(fig, "26", "case26_ehd_pumping.png")


def plot_case27():
    """Case 27: Hartmann Flow — Lorentz-braked channel flow"""
    print("Case 27: MHD Hartmann Flow")
    efile = case_path("27", "case27_hartmann_flow_out.e")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case27")
        return

    ds = open_exodus(efile)
    cx, cy = get_elem_centroids_2d(ds)
    names = get_elem_var_names(ds)

    vel_x_name = "vel_x" if "vel_x" in names else "superficial_vel_x"
    vel_y_name = "vel_y" if "vel_y" in names else "superficial_vel_y"
    p_name = "pressure"

    vx_idx = names.index(vel_x_name) + 1
    vy_idx = names.index(vel_y_name) + 1
    p_idx = names.index(p_name) + 1

    vx = get_elem_var(ds, vx_idx, timestep=-1)
    vy = get_elem_var(ds, vy_idx, timestep=-1)
    p = get_elem_var(ds, p_idx, timestep=-1)
    ds.close()

    # Plot 1: vel_x(y) profile at the channel midpoint (x ~ 2.5)
    # Extract a vertical slice at x ~ 2.5
    x_mid = 2.5
    tol_x = 0.15
    mid_mask = np.abs(cx - x_mid) < tol_x
    cy_mid = cy[mid_mask]
    vx_mid = vx[mid_mask]
    order = np.argsort(cy_mid)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Case 27: MHD Hartmann Flow (Ha = 5)", fontsize=FONT_TITLE)

    axes[0].plot(cy_mid[order], vx_mid[order], "b-o", markersize=3, label="MOOSE")
    # Analytical Hartmann profile
    Ha = 5.0
    y_an = np.linspace(0, 1, 200)
    # The actual profile shape depends on the pressure gradient that
    # develops; we normalise to match the MOOSE peak.
    v_max_moose = vx_mid.max()
    v_shape = 1.0 - np.cosh(Ha * (y_an - 0.5)) / np.cosh(Ha / 2.0)
    v_shape_max = v_shape.max()
    if v_shape_max > 0:
        v_analytical = v_max_moose * v_shape / v_shape_max
        axes[0].plot(y_an, v_analytical, "r--", linewidth=1.5, label="Hartmann (analytical)")
    axes[0].set_xlabel("y (channel width)")
    axes[0].set_ylabel("v_x")
    axes[0].set_title("Velocity Profile (mid-channel)")
    axes[0].legend()
    axes[0].grid(True, alpha=0.4)

    # Right panel: 2D velocity field
    vel_mag = np.sqrt(vx ** 2 + vy ** 2)
    tcf = axes[1].tricontourf(cx, cy, vel_mag, levels=20, cmap="plasma")
    plt.colorbar(tcf, ax=axes[1], label="|velocity|")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("y")
    axes[1].set_title("Velocity Magnitude (full domain)")

    fig.tight_layout()
    save_fig(fig, "27", "case27_hartmann_flow.png")


def plot_case28():
    """Case 28: Two-Way Joule Heating — T-dependent sigma"""
    print("Case 28: Two-Way Joule Heating")
    efile = case_path("28", "case28_twoway_joule_heating_out.e")
    cfile = case_path("28", "case28_twoway_joule_heating_out.csv")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case28")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    names = get_nod_var_names(ds)

    V_idx = names.index("V") + 1
    T_idx = names.index("T") + 1

    V = get_nod_var(ds, V_idx, timestep=-1)
    T_final = get_nod_var(ds, T_idx, timestep=-1)

    # Plot 1: side-by-side V and T at final time
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Case 28: Two-Way Joule Heating — Final State", fontsize=FONT_TITLE)

    tcf1 = axes[0].tricontourf(x, y, V, levels=20, cmap="RdYlBu_r")
    plt.colorbar(tcf1, ax=axes[0], label="V (Volts)")
    axes[0].set_title("Electric Potential V")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    axes[0].set_aspect("equal")

    tcf2 = axes[1].tricontourf(x, y, T_final, levels=20, cmap=CMAP_TEMP)
    plt.colorbar(tcf2, ax=axes[1], label="T (K)")
    axes[1].set_title("Temperature T")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("y")
    axes[1].set_aspect("equal")

    ds.close()
    fig.tight_layout()
    save_fig(fig, "28", "case28_twoway_joule_heating.png")

    # Plot 2: temperature history from CSV
    if file_exists(cfile):
        data = read_csv(cfile)
        t_csv = data["time"]
        max_T = data.get("max_T", [])
        avg_T = data.get("avg_T", [])

        if max_T and avg_T:
            fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
            ax.plot(t_csv, max_T, "r-o", markersize=3, label="max T")
            ax.plot(t_csv, avg_T, "b-s", markersize=3, label="avg T")
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Temperature (K)")
            ax.set_title("Case 28: Temperature Rise — T-Dependent Conductivity")
            ax.legend()
            ax.grid(True, alpha=0.4)
            fig.tight_layout()
            save_fig(fig, "28", "case28_temperature_history.png")


def plot_case29():
    """Case 29: Electroconvection — EHD-enhanced natural convection"""
    print("Case 29: Electroconvection")
    efile = case_path("29", "case29_electroconvection_out.e")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case29")
        return

    ds = open_exodus(efile)
    cx, cy = get_elem_centroids_2d(ds)
    names = get_elem_var_names(ds)

    vel_x_name = "vel_x" if "vel_x" in names else "superficial_vel_x"
    vel_y_name = "vel_y" if "vel_y" in names else "superficial_vel_y"
    T_name = "T_fluid"
    p_name = "pressure"

    vx_idx = names.index(vel_x_name) + 1
    vy_idx = names.index(vel_y_name) + 1
    T_idx = names.index(T_name) + 1

    vx = get_elem_var(ds, vx_idx, timestep=-1)
    vy = get_elem_var(ds, vy_idx, timestep=-1)
    T = get_elem_var(ds, T_idx, timestep=-1)
    ds.close()

    vel_mag = np.sqrt(vx ** 2 + vy ** 2)

    # Plot: temperature + velocity vectors, and velocity magnitude
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Case 29: Electroconvection — EHD-Enhanced (Fe=5)", fontsize=FONT_TITLE)

    tcf1 = axes[0].tricontourf(cx, cy, T, levels=20, cmap=CMAP_TEMP)
    plt.colorbar(tcf1, ax=axes[0], label="T_fluid")
    n = len(cx)
    step = max(1, n // 300)
    axes[0].quiver(cx[::step], cy[::step], vx[::step], vy[::step],
                   color="black", scale=max(vel_mag) * 8, width=0.003, alpha=0.6)
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    axes[0].set_aspect("equal")
    axes[0].set_title("Temperature + Velocity Vectors")

    tcf2 = axes[1].tricontourf(cx, cy, vel_mag, levels=20, cmap="plasma")
    plt.colorbar(tcf2, ax=axes[1], label="|velocity|")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("y")
    axes[1].set_aspect("equal")
    axes[1].set_title("Velocity Magnitude")

    fig.tight_layout()
    save_fig(fig, "29", "case29_electroconvection.png")


def plot_case30():
    """Case 30: Rectangular Waveguide Cutoff Frequencies"""
    print("Case 30: Rectangular Waveguide Cutoff")
    efile = case_path("30", "case30_waveguide_cutoff_out.e")
    cfile = case_path("30", "case30_waveguide_cutoff_out_eigenvalues_0002.csv")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case30")
        return

    # Plot mode shape from Exodus
    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    names = get_nod_var_names(ds)
    psi_idx = names.index("psi") + 1
    psi = get_nod_var(ds, psi_idx, timestep=-1)
    ds.close()

    # Analytical eigenvalues for a=2, b=1
    analytical = {
        "TM\u2081\u2081": (1*math.pi/2)**2 + (1*math.pi/1)**2,
        "TM\u2082\u2081": (2*math.pi/2)**2 + (1*math.pi/1)**2,
        "TM\u2083\u2081": (3*math.pi/2)**2 + (1*math.pi/1)**2,
    }

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    fig.suptitle("Case 30: Rectangular Waveguide TM Cutoff Modes (Haus Ch 2)", fontsize=FONT_TITLE)

    # Left: mode shape (first eigenmode = psi at final)
    tcf = axes[0].tricontourf(x, y, psi, levels=20, cmap="RdBu_r")
    plt.colorbar(tcf, ax=axes[0], label="\u03c8")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    axes[0].set_aspect("equal")
    axes[0].set_title("Fundamental TM\u2081\u2081 Mode Shape")

    # Right: eigenvalue comparison bar chart
    if file_exists(cfile):
        csv_data = read_csv(cfile)
        # Eigenvalues CSV has column "eigenvalues_0" or similar
        eig_key = [k for k in csv_data if "eigen" in k.lower()]
        if eig_key:
            computed = csv_data[eig_key[0]][:3]
            labels = list(analytical.keys())
            ana_vals = list(analytical.values())
            x_pos = np.arange(len(labels))
            w = 0.35
            axes[1].bar(x_pos - w/2, ana_vals[:len(labels)], w, label="Analytical", color="steelblue")
            axes[1].bar(x_pos + w/2, computed[:len(labels)], w, label="MOOSE", color="coral")
            axes[1].set_xticks(x_pos)
            axes[1].set_xticklabels(labels)
            axes[1].set_ylabel("k_c\u00b2")
            axes[1].legend()
            axes[1].set_title("Eigenvalue Comparison")
            axes[1].grid(True, alpha=0.3)
        else:
            axes[1].text(0.5, 0.5, "Eigenvalue CSV\nnot parsed", ha="center", va="center", transform=axes[1].transAxes)
    else:
        axes[1].text(0.5, 0.5, "Eigenvalue CSV\nnot found", ha="center", va="center", transform=axes[1].transAxes)

    fig.tight_layout()
    save_fig(fig, "30", "case30_waveguide_cutoff.png")


def plot_case31():
    """Case 31: Driven Resonant Cavity"""
    print("Case 31: Driven Resonant Cavity")
    efile = case_path("31", "case31_driven_cavity_out.e")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case31")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    names = get_nod_var_names(ds)
    E_idx = names.index("E") + 1
    E = get_nod_var(ds, E_idx, timestep=-1)
    ds.close()

    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
    tcf = ax.tricontourf(x, y, E, levels=30, cmap="RdBu_r")
    plt.colorbar(tcf, ax=ax, label="E field")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.set_title("Case 31: Driven Cavity \u2014 Near Resonance (k\u00b2=12.3)")
    fig.tight_layout()
    save_fig(fig, "31", "case31_driven_cavity.png")


def plot_case32():
    """Case 32: Dielectric Slab Reflection"""
    print("Case 32: Dielectric Slab Reflection")
    efile = case_path("32", "case32_dielectric_slab_out.e")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case32")
        return

    ds = open_exodus(efile)
    x = np.array(ds.variables["coordx"][:], dtype=float)
    names = get_nod_var_names(ds)
    Er_idx = names.index("E_real") + 1
    Ei_idx = names.index("E_imag") + 1
    Er = get_nod_var(ds, Er_idx, timestep=-1)
    Ei = get_nod_var(ds, Ei_idx, timestep=-1)
    ds.close()

    E_mag = np.sqrt(Er**2 + Ei**2)
    order = np.argsort(x)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    fig.suptitle("Case 32: EM Wave Reflection from Dielectric Slab (Haus Ch 1)", fontsize=FONT_TITLE)

    # Left: real and imaginary components
    axes[0].plot(x[order], Er[order], "b-", label="E_real", linewidth=1)
    axes[0].plot(x[order], Ei[order], "r--", label="E_imag", linewidth=1)
    axes[0].axvspan(0, 25, alpha=0.1, color="gray", label="Slab (\u03b5\u1d63=4)")
    axes[0].set_xlabel("x (m)")
    axes[0].set_ylabel("E field")
    axes[0].legend()
    axes[0].set_title("Real and Imaginary Components")
    axes[0].grid(True, alpha=0.3)

    # Right: magnitude
    axes[1].plot(x[order], E_mag[order], "k-", linewidth=1.5)
    axes[1].axvspan(0, 25, alpha=0.1, color="gray", label="Slab (\u03b5\u1d63=4)")
    axes[1].set_xlabel("x (m)")
    axes[1].set_ylabel("|E|")
    axes[1].legend()
    axes[1].set_title("Field Magnitude")
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    save_fig(fig, "32", "case32_dielectric_slab.png")


def plot_case33():
    """Case 33: Coupled Resonator Beating"""
    print("Case 33: Coupled Resonator Beating")
    cfile = case_path("33", "case33_coupled_resonators_out.csv")
    if not file_exists(cfile):
        print("    SKIP: output file not found")
        skipped_cases.append("case33")
        return

    csv_data = read_csv(cfile)
    t = csv_data["time"]
    avg_u = csv_data["avg_u"]
    avg_v = csv_data["avg_v"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    fig.suptitle("Case 33: Coupled Resonator Beating (Haus Ch 3)", fontsize=FONT_TITLE)

    # Left: time series
    axes[0].plot(t, avg_u, "b-", label="\u27e8u\u27e9 (mode 1)", linewidth=1.2)
    axes[0].plot(t, avg_v, "r--", label="\u27e8v\u27e9 (mode 2)", linewidth=1.2)
    axes[0].set_xlabel("Time (s)")
    axes[0].set_ylabel("Average amplitude")
    axes[0].legend()
    axes[0].set_title("Mode Amplitudes \u2014 Beating")
    axes[0].grid(True, alpha=0.3)

    # Right: envelope showing exponential decay
    total = [u**2 + v**2 for u, v in zip(avg_u, avg_v)]
    axes[1].plot(t, total, "k-", linewidth=1.2, label="\u27e8u\u27e9\u00b2 + \u27e8v\u27e9\u00b2")
    # Overlay exponential envelope
    t_arr = np.array(t)
    if len(t_arr) > 1 and total[0] > 0:
        envelope = total[0] * np.exp(-2 * 0.5 * t_arr)  # decay at 2*gamma
        axes[1].plot(t_arr, envelope, "g--", linewidth=1, label="e^{-2\u03b3t} envelope")
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Energy (arb.)")
    axes[1].legend()
    axes[1].set_title("Total Energy Decay")
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    save_fig(fig, "33", "case33_coupled_resonators.png")


def plot_case34():
    """Case 34: Thermal Noise Relaxation"""
    print("Case 34: Thermal Noise Relaxation")
    efile = case_path("34", "case34_thermal_noise_out.e")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case34")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    times = get_times(ds)
    names = get_nod_var_names(ds)
    T_idx = names.index("T") + 1

    snap_times = [0.0, 0.1, 0.5, 2.0]
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.5))
    fig.suptitle("Case 34: Thermal Noise Relaxation (Haus Ch 5)", fontsize=FONT_TITLE)

    for i, st in enumerate(snap_times):
        ti = find_closest_timestep(times, st)
        T = get_nod_var(ds, T_idx, timestep=ti)
        tcf = axes[i].tricontourf(x, y, T, levels=20, cmap=CMAP_TEMP, vmin=0.3, vmax=0.7)
        axes[i].set_aspect("equal")
        axes[i].set_title(f"t = {times[ti]:.2f}")
        if i == 0:
            axes[i].set_ylabel("y")
        axes[i].set_xlabel("x")
    plt.colorbar(tcf, ax=axes[-1], label="T")
    ds.close()

    fig.tight_layout()
    save_fig(fig, "34", "case34_thermal_noise.png")


def plot_case35():
    """Case 35: Dispersive Pulse Broadening"""
    print("Case 35: Dispersive Pulse Broadening")
    efile = case_path("35", "case35_dispersive_pulse_out.e")
    cfile = case_path("35", "case35_dispersive_pulse_out.csv")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case35")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    times = get_times(ds)
    names = get_nod_var_names(ds)
    A_idx = names.index("A") + 1

    # Extract midline (y near 0.1)
    mid = nodes_near_y(y, 0.1, tol=0.06)

    snap_times = [0.0, 2.0, 4.0, 6.0]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    fig.suptitle("Case 35: Dispersive Pulse Broadening (Haus Ch 4)", fontsize=FONT_TITLE)

    colors = ["blue", "green", "orange", "red"]
    for i, st in enumerate(snap_times):
        ti = find_closest_timestep(times, st)
        A = get_nod_var(ds, A_idx, timestep=ti)
        xm = x[mid]
        Am = A[mid]
        order = np.argsort(xm)
        axes[0].plot(xm[order], Am[order], color=colors[i], linewidth=1.2,
                     label=f"t={times[ti]:.1f}")
    ds.close()

    axes[0].set_xlabel("x")
    axes[0].set_ylabel("A")
    axes[0].legend()
    axes[0].set_title("Pulse Profiles at Different Times")
    axes[0].grid(True, alpha=0.3)

    # Right: CSV postprocessor time series
    if file_exists(cfile):
        csv_data = read_csv(cfile)
        t = csv_data["time"]
        axes[1].plot(t, csv_data.get("max_A", [0]*len(t)), "b-", label="max(A)")
        axes[1].plot(t, csv_data.get("total_A", [0]*len(t)), "r--", label="\u222bA dV")
        axes[1].set_xlabel("Time")
        axes[1].set_ylabel("Value")
        axes[1].legend()
        axes[1].set_title("Peak Amplitude and Total Mass")
        axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    save_fig(fig, "35", "case35_dispersive_pulse.png")


def plot_case36():
    """Case 36: Soliton Pulse Propagation"""
    print("Case 36: Soliton Pulse Propagation")
    efile = case_path("36", "case36_soliton_pulse_out.e")
    cfile = case_path("36", "case36_soliton_pulse_out.csv")
    if not file_exists(efile):
        print("    SKIP: output file not found")
        skipped_cases.append("case36")
        return

    ds = open_exodus(efile)
    x, y = get_coords_2d(ds)
    times = get_times(ds)
    names = get_nod_var_names(ds)
    A_idx = names.index("A") + 1

    mid = nodes_near_y(y, 0.1, tol=0.06)

    snap_times = [0.0, 3.0, 6.0, 9.0]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    fig.suptitle("Case 36: Soliton Pulse Propagation (Haus Ch 10)", fontsize=FONT_TITLE)

    colors = ["blue", "green", "orange", "red"]
    for i, st in enumerate(snap_times):
        ti = find_closest_timestep(times, st)
        A = get_nod_var(ds, A_idx, timestep=ti)
        xm = x[mid]
        Am = A[mid]
        order = np.argsort(xm)
        axes[0].plot(xm[order], Am[order], color=colors[i], linewidth=1.2,
                     label=f"t={times[ti]:.1f}")

    ds.close()

    axes[0].set_xlabel("x")
    axes[0].set_ylabel("A")
    axes[0].legend()
    axes[0].set_title("Soliton Profiles (\u03b1=0.1)")
    axes[0].grid(True, alpha=0.3)

    # Right: max_A time series (should stay constant for soliton)
    if file_exists(cfile):
        csv_data = read_csv(cfile)
        t = csv_data["time"]
        axes[1].plot(t, csv_data.get("max_A", [0]*len(t)), "b-", linewidth=1.5,
                     label="max(A)")
        axes[1].axhline(y=1.0, color="gray", linestyle="--", alpha=0.5, label="A\u2080=1.0")
        axes[1].set_xlabel("Time")
        axes[1].set_ylabel("max(A)")
        axes[1].legend()
        axes[1].set_title("Peak Amplitude vs Time")
        axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    save_fig(fig, "36", "case36_soliton_pulse.png")


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

CASE_FUNCTIONS = [
    ("01", plot_case01),
    ("02", plot_case02),
    ("03", plot_case03),
    ("04", plot_case04),
    ("05", plot_case05),
    ("06", plot_case06),
    ("07", plot_case07),
    ("08", plot_case08),
    ("09", plot_case09),
    ("10", plot_case10),
    ("11", plot_case11),
    ("12", plot_case12),
    ("13", plot_case13),
    ("14", plot_case14),
    ("15", plot_case15),
    ("16", plot_case16),
    ("17", plot_case17),
    ("18", plot_case18),
    ("19", plot_case19),
    ("20", plot_case20),
    ("21", plot_case21),
    ("22", plot_case22),
    ("23", plot_case23),
    ("24", plot_case24),
    ("25", plot_case25),
    ("26", plot_case26),
    ("27", plot_case27),
    ("28", plot_case28),
    ("29", plot_case29),
    ("30", plot_case30),
    ("31", plot_case31),
    ("32", plot_case32),
    ("33", plot_case33),
    ("34", plot_case34),
    ("35", plot_case35),
    ("36", plot_case36),
]


if __name__ == "__main__":
    print("=" * 60)
    print("MOOSE Quick-Start Visualization — All 36 Cases")
    print("=" * 60)
    print(f"Script directory: {SCRIPT_DIR}")
    print("Plots will be saved into each case subdirectory.\n")

    # Run each case
    for case_id, func in CASE_FUNCTIONS:
        try:
            func()
        except Exception:
            print(f"    ERROR in case {case_id}:")
            traceback.print_exc()
            skipped_cases.append(f"case{case_id} (exception)")
        print()

    # Summary
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"\nPlots generated ({len(generated_plots)} total):")
    for path in generated_plots:
        print(f"  {os.path.relpath(path, SCRIPT_DIR)}")

    if skipped_cases:
        print(f"\nSkipped ({len(skipped_cases)}):")
        for name in skipped_cases:
            print(f"  {name}")
    else:
        print("\nAll cases processed successfully.")

    print("\nDone.")
