"""
MOOSE Quick-Start: Visualization Script for All 13 Cases
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
]


if __name__ == "__main__":
    print("=" * 60)
    print("MOOSE Quick-Start Visualization — All 13 Cases")
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
