# Engine Cycle Simulation (MATLAB)

This project simulates the in-cylinder cycle of a spark-ignition engine using a **single-zone thermodynamic model** with a **Wiebe heat release law**.  
It iteratively adjusts intake pressure, exhaust temperature, and residual burn fraction until consistency criteria are met.

## Features
- Crank-angle resolved volume, pressure, temperature, and heat release.
- **Wiebe function** for combustion (`a`, `n`, `ts`, `td`).
- Iterative convergence on:
  - Intake pressure `Pi`
  - Exhaust temperature `Te`
  - Burned fraction parameter `f1`
- Numerical integration of indicated work (trapz over P–V curve).
- Automatic plots:
  - Cylinder pressure vs crank angle
  - Temperature vs crank angle (full-cycle)
  - P–V diagram
- Printed results: work, efficiency, Pi, Pe, Te, f1, error metrics.

## File
- `Engine_Cycle_Simulation.m` → single complete script (no external dependencies, runs as-is in MATLAB).

## Workflow
1. Open `engine_cycle.m` in MATLAB.
2. Adjust the input parameters:
   - Fuel properties (`qc`, `A`, `F`, `as`, `Mf`).
   - Engine geometry (`VD`, `rc`, `N`, rod ratio).
   - Combustion parameters (`ts`, `td`, `a`, `nWie`).
3. Run the script.
4. The while-loop iteratively updates `Pi`, `Te`, and `f1` until all error metrics (`r1, r2, r3`) fall below tolerance `ae` or `maxIter` is reached.
5. Results are printed in the console and plots are displayed.

---

## Parameters you may tune
- `ts` → Start of combustion [°CA]
- `td` → Combustion duration [°CA]
- `a` , `nWie` → Wiebe law shape parameters
- `rc` → Compression ratio
- `N` → Engine speed [rpm]
- `Pi`, `Pe` → Intake and exhaust pressures [bar]
- `Te` → Initial exhaust gas temperature guess [K]
- `f1` → Initial burned fraction parameter

## Outputs
- **Printed values:**
  - Iterations until convergence
  - Indicated work per cycle [J]
  - Indicated efficiency [%]
  - Converged `Pi`, `Pe`, `Te`, `f1`
  - Final error values (`r1`, `r2`, `r3`)
- **Plots:**
  - Pressure vs crank angle
  - Temperature vs crank angle
  - P–V diagram

## Notes
- Units are SI, except `Pi` and `Pe` which are in bar on input.
- Default parameters are placeholders — tune for real engine data.
- Model is **single-zone** with constant specific heat ratio γ (no heat losses, no variable γ).
- Convergence uses signed, damped updates to ensure stability.

## Limitations & Extensions
- No heat-transfer model (can be added with Woschni correlation).
- Constant γ — could be replaced by temperature-dependent cp/cv.
- Single Wiebe — can be extended to double-Wiebe for premixed + diffusion burn.
- No pumping/volumetric efficiency — `bmep` is treated as IMEP target proxy.
