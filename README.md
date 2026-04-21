# EE414 Assignment III — Numerical Methods for ODEs

**Author:** Fortune Olose
**Student ID:** 21455956
**Course:** EE414

This repository contains the submission materials for EE414 Assignment III. The
assignment solves a first-order linear ODE arising from an RL-circuit initial
value problem using three numerical methods: Euler, Heun's RK2, and classical
RK4. Results are compared against the analytic solution for step sizes
h = 40, 20, and 10 ps over t ∈ [0, 240] ps.

## Problem

RL circuit with R = 4 MΩ, L = 580 µH, giving time constant τ = L/R = 145 ps.

IVP:

    dy/dt = (1 - y) / τ,    y(0) = 0

Analytic solution:

    y(t) = 1 - exp(-t/τ)

## Repository Contents

### Assignment brief

| File | Description |
|---|---|
| `Assignment 3 instructions (1).pdf` | Original assignment specification |

### Written report (Q1–Q5)

| File | Description |
|---|---|
| `EE414_Assignment_Q1_Q5.docx` | Full written solution, editable Word document |
| `EE414_Assignment_Q1_Q5_preview.pdf` | PDF export of the report for viewing |

### MATLAB source (numerical solvers)

| File | Description |
|---|---|
| `EE414_ODE_solvers.m` | Combined driver script: runs Euler, RK2, and RK4 at h = 40, 20, 10 ps; prints comparison tables and generates the convergence plot. Contains local functions `euler`, `rk2`, `rk4`. |
| `EE414_Q3_Euler.m` | Q3 — standalone Euler's method, tabulates ŷ(tᵢ), y(tᵢ), and absolute error Eₜ for each step size |
| `EE414_Q3_Euler_1.m` … `EE414_Q3_Euler_4.m` | Iterative revisions of the Q3 Euler script (history of formatting/output refinements) |
| `EE414_Q4_RK2.m` | Q4 — Heun's method (RK2); writes per-step tables to `Q4_RK2_Tables.xlsx` |
| `EE414_Q5_RK4.m` | Q5 — classical 4th-order Runge-Kutta; writes per-step tables to `Q5_RK4_Tables.xlsx` |

### Results (generated tables)

| File | Description |
|---|---|
| `Q3_Euler_Tables.xlsx` | Euler results, one sheet per step size (h = 40, 20, 10 ps) |
| `Q4_RK2_Tables.xlsx` | RK2 (Heun) results, one sheet per step size |
| `Q5_RK4_Tables.xlsx` | RK4 results, one sheet per step size |

Each sheet columns: `i`, `t_i` (ps), `yhat_ti`, `y_ti` (exact), `E_t` (absolute error).

### Figures

| File | Description |
|---|---|
| `convergence.png` | Log-log plot of global error |Eₜ| at t = 240 ps vs step size h, showing the expected orders of accuracy: O(h) for Euler, O(h²) for RK2, O(h⁴) for RK4 |

## Reproducing the results

1. Open MATLAB in this directory.
2. Run `EE414_ODE_solvers.m` to print the combined comparison tables and
   regenerate the convergence plot.
3. Run `EE414_Q3_Euler.m`, `EE414_Q4_RK2.m`, `EE414_Q5_RK4.m` individually to
   regenerate the per-question Excel tables.

All time values are in picoseconds.
