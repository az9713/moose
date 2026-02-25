"""
Plot postprocessor history from Case 13.

Usage:
    python plot_case13.py

Requires: matplotlib, numpy (standard scientific Python stack)
Install:  pip install matplotlib numpy
"""

import csv
import math
import sys

# ---- Read the CSV file ----
filename = "case13_postprocessors_out.csv"
try:
    with open(filename, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
except FileNotFoundError:
    print(f"ERROR: '{filename}' not found.")
    print("Run the MOOSE input file first: ./moose_test-opt -i case13_postprocessors.i")
    sys.exit(1)

# ---- Extract columns ----
time         = [float(r["time"])         for r in rows]
max_temp     = [float(r["max_temp"])     for r in rows]
avg_temp     = [float(r["avg_temp"])     for r in rows]
total_energy = [float(r["total_energy"]) for r in rows]
T_L2_norm    = [float(r["T_L2_norm"])   for r in rows]
dt_size      = [float(r["dt_size"])     for r in rows]

# ---- Analytical steady-state check ----
# For -k*div(grad T) = Q with T=0 on boundary of [0,1]^2,
# the steady-state average temperature is Q / (k * 8*pi^2) * 16
# (first eigenfunction approximation: Q*8/(k*pi^4))
k = 2.0
Q = 5.0
T_ss_approx = Q * 8.0 / (k * math.pi**4)
print(f"Analytical steady-state avg T (1-term approximation): {T_ss_approx:.4f}")
print(f"Final simulated avg T: {avg_temp[-1]:.4f}")

# ---- Print a table ----
print(f"\n{'time':>8} {'max_temp':>10} {'avg_temp':>10} {'energy':>12} {'dt':>10}")
print("-" * 55)
for i in range(0, len(time), max(1, len(time)//12)):
    print(f"{time[i]:8.4f} {max_temp[i]:10.5f} {avg_temp[i]:10.5f} "
          f"{total_energy[i]:12.6f} {dt_size[i]:10.6f}")

# ---- Matplotlib plots (optional, comment out if not installed) ----
try:
    import matplotlib
    matplotlib.use("Agg")   # non-interactive backend (writes PNG)
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    fig.suptitle("MOOSE Case 13: Transient Heat Equation", fontsize=14)

    # Temperature histories
    axes[0, 0].plot(time, max_temp, "r-o", markersize=3, label="Max T")
    axes[0, 0].plot(time, avg_temp, "b-s", markersize=3, label="Avg T")
    axes[0, 0].axhline(T_ss_approx, color="gray", linestyle="--", label="Steady-state avg (analytical)")
    axes[0, 0].set_xlabel("Time (s)")
    axes[0, 0].set_ylabel("Temperature")
    axes[0, 0].set_title("Temperature vs. Time")
    axes[0, 0].legend()
    axes[0, 0].grid(True)

    # Total energy
    axes[0, 1].plot(time, total_energy, "g-^", markersize=3)
    axes[0, 1].set_xlabel("Time (s)")
    axes[0, 1].set_ylabel("int(T) dV")
    axes[0, 1].set_title("Total Thermal Energy")
    axes[0, 1].grid(True)

    # L2 norm
    axes[1, 0].plot(time, T_L2_norm, "m-D", markersize=3)
    axes[1, 0].set_xlabel("Time (s)")
    axes[1, 0].set_ylabel("||T||_L2")
    axes[1, 0].set_title("L2 Norm of Temperature")
    axes[1, 0].grid(True)

    # Adaptive timestep size
    axes[1, 1].semilogy(time, dt_size, "k-x", markersize=3)
    axes[1, 1].set_xlabel("Time (s)")
    axes[1, 1].set_ylabel("dt (log scale)")
    axes[1, 1].set_title("Adaptive Timestep Size")
    axes[1, 1].grid(True)

    plt.tight_layout()
    outpng = "case13_results.png"
    plt.savefig(outpng, dpi=120)
    print(f"\nPlot saved to: {outpng}")

except ImportError:
    print("\nmatplotlib not found — skipping plots.")
    print("Install with: pip install matplotlib")
