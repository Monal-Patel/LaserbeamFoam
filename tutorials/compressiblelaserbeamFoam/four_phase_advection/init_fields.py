#!/usr/bin/env python3
"""Initialize four-phase Zalesak disk alpha fields for OpenFOAM.

Replicates the geometry from amrsolver's prob.h:
- Disk centered at (0.5, 0.5), radius 0.30, rotated 30 degrees
- 3 sectors (120 deg each) for liquid/vapor/solid
- 3 outer radial slots with rounded ends
- 3 sinusoidal notch indentations
- Gas fills the remaining domain

Phase mapping (amrsolver -> OpenFOAM):
  phase 0 -> liquid  (rho=998.2)
  phase 1 -> vapor   (rho=800.0)
  phase 2 -> solid   (rho=1200.0)
  phase 3 -> gas     (rho=1.2, ambient)
"""

import numpy as np
import os

# Geometry parameters (matching amrsolver prob.h)
CX, CY = 0.5, 0.5
RADIUS = 0.30
SLOT_W = 0.05
SLOT_H = 0.10
NOTCH_AMP = 0.025
ROTATION_DEG = 30.0

# Mesh parameters
NX, NY, NZ = 100, 100, 4
X_LO, X_HI = 0.0, 1.0
Y_LO, Y_HI = 0.0, 1.0


def compute_alpha(x, y):
    """Compute phase volume fractions at cell center (x, y).

    Returns (alpha_liquid, alpha_vapor, alpha_solid).
    alpha_gas = 1 - sum.
    """
    theta = np.radians(ROTATION_DEG)
    ct, st = np.cos(theta), np.sin(theta)

    xr = (x - CX) * ct + (y - CY) * st
    yr = -(x - CX) * st + (y - CY) * ct

    r = np.sqrt(xr * xr + yr * yr)
    phi = np.arctan2(yr, xr)
    in_disk = r <= RADIUS

    notch_end = (RADIUS - SLOT_H) * 0.5
    base_angles = [-np.pi / 2, -np.pi / 2 + 2 * np.pi / 3, -np.pi / 2 + 4 * np.pi / 3]

    in_outer_slots = False
    in_notch = [False, False, False]

    for s in range(3):
        ca, sa = np.cos(base_angles[s]), np.sin(base_angles[s])
        xl = xr * ca + yr * sa
        yl = -xr * sa + yr * ca

        in_rect = (abs(yl) <= 0.5 * SLOT_W) and (xl >= RADIUS - SLOT_H)
        ccx = RADIUS - SLOT_H
        in_round = (xl < ccx) and ((xl - ccx) ** 2 + yl ** 2 <= 0.25 * SLOT_W * SLOT_W)
        in_outer_slots = in_outer_slots or in_rect or in_round

        in_notch[s] = (
            (xl >= 0.0)
            and (xl <= notch_end)
            and (abs(yl) <= NOTCH_AMP * np.sin(np.pi * xl / notch_end) if notch_end > 0 else False)
        )

    in_any_notch = in_notch[0] or in_notch[1] or in_notch[2]
    solid_disk = in_disk and (not in_outer_slots) and (not in_any_notch)

    alpha0, alpha1, alpha2 = 0.0, 0.0, 0.0

    if solid_disk:
        sp = (phi + np.pi / 2 + 4 * np.pi) % (2 * np.pi)
        if sp < 2 * np.pi / 3:
            alpha0 = 1.0
        elif sp < 4 * np.pi / 3:
            alpha1 = 1.0
        else:
            alpha2 = 1.0

    if in_notch[0]:
        alpha0, alpha1, alpha2 = 1.0, 0.0, 0.0
    if in_notch[1]:
        alpha0, alpha1, alpha2 = 0.0, 1.0, 0.0
    if in_notch[2]:
        alpha0, alpha1, alpha2 = 0.0, 0.0, 1.0

    return alpha0, alpha1, alpha2


def write_of_field(filepath, field_name, values, nx, ny, nz):
    """Write OpenFOAM volScalarField in ascii format."""
    total = nx * ny * nz
    with open(filepath, "w") as f:
        f.write(f"FoamFile {{ version 2.0; format ascii; class volScalarField; location \"0\"; object {field_name}; }}\n")
        f.write("dimensions [0 0 0 0 0 0 0];\n")
        f.write(f"internalField nonuniform List<scalar>\n{total}\n(\n")
        for v in values:
            f.write(f"{v}\n")
        f.write(");\n")
        f.write("boundaryField { periodic_x_left { type cyclic; } periodic_x_right { type cyclic; } periodic_y_bottom { type cyclic; } periodic_y_top { type cyclic; } periodic_z_back { type cyclic; } periodic_z_front { type cyclic; }}\n")


def main():
    dx = (X_HI - X_LO) / NX
    dy = (Y_HI - Y_LO) / NY

    alpha_liquid = []
    alpha_vapor = []
    alpha_solid = []
    alpha_gas = []

    # OpenFOAM cell ordering: z fastest, then y, then x
    # For blockMesh hex (0 1 2 3 4 5 6 7) with (nx ny nz):
    # index = k + nz * (j + ny * i) ... but actually OpenFOAM uses (i + nx*(j + ny*k))
    # Standard OpenFOAM ordering: x varies fastest, then y, then z
    for k in range(NZ):
        for j in range(NY):
            for i in range(NX):
                x = X_LO + (i + 0.5) * dx
                y = Y_LO + (j + 0.5) * dy

                a0, a1, a2 = compute_alpha(x, y)
                a3 = 1.0 - a0 - a1 - a2

                alpha_liquid.append(a0)
                alpha_vapor.append(a1)
                alpha_solid.append(a2)
                alpha_gas.append(a3)

    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "0")
    write_of_field(os.path.join(out_dir, "alpha.liquid"), "alpha.liquid", alpha_liquid, NX, NY, NZ)
    write_of_field(os.path.join(out_dir, "alpha.vapor"), "alpha.vapor", alpha_vapor, NX, NY, NZ)
    write_of_field(os.path.join(out_dir, "alpha.solid"), "alpha.solid", alpha_solid, NX, NY, NZ)
    write_of_field(os.path.join(out_dir, "alpha.gas"), "alpha.gas", alpha_gas, NX, NY, NZ)

    n_liq = sum(1 for v in alpha_liquid if v > 0.5)
    n_vap = sum(1 for v in alpha_vapor if v > 0.5)
    n_sol = sum(1 for v in alpha_solid if v > 0.5)
    n_gas = sum(1 for v in alpha_gas if v > 0.5)
    total = NX * NY * NZ
    print(f"Initialized {total} cells:")
    print(f"  liquid: {n_liq} cells ({100*n_liq/total:.1f}%)")
    print(f"  vapor:  {n_vap} cells ({100*n_vap/total:.1f}%)")
    print(f"  solid:  {n_sol} cells ({100*n_sol/total:.1f}%)")
    print(f"  gas:    {n_gas} cells ({100*n_gas/total:.1f}%)")


if __name__ == "__main__":
    main()
