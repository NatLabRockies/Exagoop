#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_particles.py
Generate MPM particle file (mpm_particles.dat) and inputs_axialbar for ExaGOOP.

Positional arguments (original ordering kept):
  1  no_of_cell_in_x         (int >0)
  2  buffery                 (int >=0)
  3  periodic                (0 or 1)
  4  np_per_cell_x           (int >0)
  5  order_scheme            (int)
  6  alpha_pic_flip          (float)
  7  stress_update_scheme    (int)
  8  CFL                     (float)  
 9  output_tag              (string) used for Solution/<tag> prefixes

Example:
  python3 generate_particles.py 25 3 1 1 3 1 1 0.1 mytag
"""
from __future__ import annotations
import argparse
import sys
import numpy as np
import tempfile
import shutil
import os
from typing import List, Tuple


def die(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(1)


def warn(msg: str) -> None:
    print(f"Warning: {msg}", file=sys.stderr)


def write_atomic_with_count(lines: List[str], out_filename: str) -> None:
    """
    Write lines to out_filename with the first line equal to the count.
    Uses a temporary file and atomic move.
    """
    out_dir = os.path.dirname(os.path.abspath(out_filename)) or "."
    os.makedirs(out_dir, exist_ok=True)
    tmp = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=out_dir, delete=False) as tmpf:
            tmp = tmpf.name
            tmpf.write(f"{len(lines)}\n")
            for ln in lines:
                tmpf.write(ln)
        shutil.move(tmp, out_filename)
    except Exception as e:
        if tmp and os.path.exists(tmp):
            try:
                os.remove(tmp)
            except Exception:
                pass
        die(f"Failed to write '{out_filename}': {e}")


def generate_particles_and_return(ncells_x: int,
                                  buffery: int,
                                  periodic: int,
                                  np_per_cell_x: int,
                                  order_scheme: int,
                                  alpha_pic_flip: float,
                                  stress_update_scheme: int,
                                  CFL: float,
                                  output_tag: str,
                                  out_particles: str = "mpm_particles.dat"
                                  ) -> Tuple[int, float]:
    """
    Generate particle lines and write mpm_particles.dat atomically.
    Returns (npart, dx1).
    """
    if ncells_x <= 0:
        die("no_of_cell_in_x must be > 0")
    if buffery < 0:
        die("buffery must be >= 0")
    if np_per_cell_x <= 0:
        die("np_per_cell_x must be > 0")
    if periodic not in (0, 1):
        warn("periodic should be 0 or 1; treating other values as 0")
        periodic = 0

    # domain geometry (kept from original script)
    L = 25.0
    bufferx = 4
    bufferz = buffery
    dx1 = L / float(ncells_x)

    blo = np.array([0.0, -(buffery + 0.5) * dx1, -(bufferz + 0.5) * dx1], dtype=float)
    bhi = np.array([L + bufferx * dx1, (buffery + 0.5) * dx1, (bufferz + 0.5) * dx1], dtype=float)
    ncells = np.array([ncells_x + bufferx, 2 * buffery + 1, 2 * bufferz + 1], dtype=int)
    dx = (bhi - blo) / ncells

    if not np.allclose(dx[0], dx):
        warn(f"mesh sizes differ: dx = {dx}")

    # sampling window for placing particles (same as original)
    xmin, xmax = 0.0, L
    ymin, ymax = -dx1 / 2.0, dx1 / 2.0
    zmin, zmax = -dx1 / 2.0, dx1 / 2.0

    # particle geometry and physical properties (original defaults)
    vol_cell = float(dx[0] * dx[1] * dx[2])
    vol_particle = vol_cell / float(np_per_cell_x * 1 * 1)
    if vol_particle <= 0.0 or not np.isfinite(vol_particle):
        die("Computed particle volume is non-positive or invalid.")
    rad = (3.0 / 4.0 * vol_particle / np.pi) ** (1.0 / 3.0)

    dens = 1.0
    phase = 0
    E = 100.0
    nu = 0.0
    v0 = 0.1
    n = 1
    beta_n = (2 * n - 1) / 2.0 * np.pi / L

    particle_lines: List[str] = []

    # generate particle lines; keep formatting stable so 'phase' prints as an integer token
    for k in range(int(ncells[2])):
        for j in range(int(ncells[1])):
            for i in range(int(ncells[0])):
                c_cx = blo[0] + i * dx[0]
                c_cy = blo[1] + j * dx[1]
                c_cz = blo[2] + k * dx[2]
                if (xmin <= c_cx < xmax) and (ymin <= c_cy < ymax) and (zmin <= c_cz < zmax):
                    for ii in range(int(np_per_cell_x)):
                        cell_cx = c_cx + (2 * ii + 1) * dx[0] / (2.0 * np_per_cell_x)
                        velx = v0 * np.sin(beta_n * cell_cx)
                        vely = 0.0
                        velz = 0.0
                        # explicit formatting keeps tokens separated and phase as integer
                        line = "{phase:d} {cx:.6e} {cy:.6e} {cz:.6e} {rad:.6e} {dens:.6e} {vx:.6e} {vy:.6e} {vz:.6e} {flag:d} {E:.6e} {nu:.6e}\n".format(
                            phase=int(phase),
                            cx=cell_cx,
                            cy=0.0,
                            cz=0.0,
                            rad=rad,
                            dens=dens,
                            vx=velx,
                            vy=vely,
                            vz=velz,
                            flag=0,
                            E=E,
                            nu=nu
                        )
                        particle_lines.append(line)

    npart = len(particle_lines)
    if npart == 0:
        warn("Zero particles generated; check domain/parameters.")

    # write atomically (header + particle lines)
    write_atomic_with_count(particle_lines, out_particles)
    print(f"WROTE: {out_particles} with {npart} particles")
    return npart, dx1


def write_inputs_file(ncells_x: int,
                      buffery: int,
                      periodic: int,
                      np_per_cell_x: int,
                      order_scheme: int,
                      alpha_pic_flip: float,
                      stress_update_scheme: int,
                      CFL: float,                      
                      output_tag: str,
                      dx1: float,
                      out_filename: str = "inputs_axialbar") -> None:
    """
    Write inputs_axialbar file using provided parameters. Creates Solution/<output_tag> directory if needed.
    """
    bufferx = 4
    xmin = 0.0
    xmax = 25.0 + bufferx * dx1
    ymin = -dx1 * (buffery + 0.5)
    ymax = dx1 * (buffery + 0.5)
    zmin = -dx1 * (buffery + 0.5)
    zmax = dx1 * (buffery + 0.5)

    sol_dir = os.path.join(".", "Solution", output_tag)
    try:
        os.makedirs(sol_dir, exist_ok=True)
    except Exception as e:
        warn(f"Could not create Solution directory '{sol_dir}': {e}")

    try:
        with open(out_filename, "w", encoding="utf-8") as f:
            f.write("#geometry parameters\n")
            f.write(f"mpm.prob_lo = {xmin} {ymin} {zmin}\t\t\t#Lower corner of physical domain\n")
            f.write(f"mpm.prob_hi = {xmax} {ymax} {zmax}\t\t\t#Upper corner of physical domain\n")
            f.write(f"mpm.ncells  = {ncells_x + bufferx} {2 * buffery + 1} {2 * buffery + 1}\n")
            f.write(f"mpm.max_grid_size = {ncells_x + bufferx + 1}\n")
            f.write(f"mpm.is_it_periodic = 0  {periodic}  {periodic}\n")

            f.write("\n\n#AMR Parameters\n")
            f.write("#restart_checkfile = \"\"\n")

            f.write("\n\n#Input files\n")
            f.write("mpm.use_autogen=0\n")
            f.write("mpm.mincoords_autogen=0.0 0.0 0.0\n")
            f.write("mpm.maxcoords_autogen=1.0 1.0 1.0\n")
            f.write("mpm.vel_autogen=0.0 0.0 0.0\n")
            f.write("mpm.constmodel_autogen=0\n")
            f.write("mpm.dens_autogen=1.0\n")
            f.write("mpm.E_autogen=1e6\n")
            f.write("mpm.nu_autogen=0.3\n")
            f.write("mpm.bulkmod_autogen=2e6\n")
            f.write("mpm.Gama_pres_autogen=7\n")
            f.write("mpm.visc_autogen=0.001\n")
            f.write("mpm.multi_part_per_cell_autogen=1\n")
            f.write("mpm.particle_file=\"mpm_particles.dat\"\n")

            f.write("\n\n#File output parameters\n")
            f.write(f"mpm.prefix_particlefilename=\"./Solution/{output_tag}/plt\"\n")
            f.write(f"mpm.prefix_gridfilename=\"./Solution/{output_tag}/nplt\"\n")
            f.write(f"mpm.prefix_densityfilename=\"./Solution/{output_tag}/dens\"\n")
            f.write(f"mpm.prefix_checkpointfilename=\"./Solution/{output_tag}/chk\"\n")
            f.write("mpm.num_of_digits_in_filenames=6\n")

            f.write("\n\n#Simulation run parameters\n")
            f.write("mpm.final_time=50.0\n")
            f.write("mpm.max_steps=5000000\n")
            f.write("mpm.screen_output_time = 0.001\n")
            f.write("mpm.write_output_time=0.5\n")
            f.write("mpm.num_redist = 1\n")

            f.write("\n\n#Timestepping parameters\n")
            f.write("mpm.fixed_timestep = 0\n")
            f.write("mpm.timestep = 1.0e-5\n")
            f.write(f"mpm.CFL={CFL}\n")
            f.write("mpm.dt_min_limit=1e-12\n")
            f.write("mpm.dt_max_limit=1e+00\n")

            f.write("\n\n#Numerical schemes\n")
            f.write(f"mpm.order_scheme={order_scheme}\n")
            f.write(f"mpm.alpha_pic_flip = {alpha_pic_flip}\n")
            f.write(f"mpm.stress_update_scheme= {stress_update_scheme}\n")
            f.write("mpm.mass_tolerance = 1e-18\n")

            f.write("\n\n#Physics parameters\n")
            f.write("mpm.gravity = 0.0 0.0 0.0\n")
            f.write("mpm.applied_strainrate_time=0.0\n")
            f.write("mpm.applied_strainrate=0.0\n")
            f.write("mpm.external_loads=0\n")
            f.write("mpm.force_slab_lo= 0.0 0.0 0.0\n")
            f.write("mpm.force_slab_hi= 1.0 1.0 1.0\n")
            f.write("mpm.extforce = 0.0 0.0 0.0\n")

            f.write("\n\n#Diagnostics and Test\n")
            f.write("mpm.print_diagnostics= 1\n")
            f.write("mpm.is_standard_test= 1\n")
            f.write("mpm.test_number= 1\n")
            f.write("mpm.axial_bar_E= 100\n")
            f.write("mpm.axial_bar_rho= 1\n")
            f.write("mpm.axial_bar_L= 25.0\n")
            f.write("mpm.axial_bar_modenumber= 1\n")
            f.write("mpm.axial_bar_v0= 0.1\n")

            f.write("\n\n#Boundary conditions\n")
            f.write("mpm.bc_lower=1 0 0\n")
            f.write("mpm.bc_upper=2 0 0\n")

    except Exception as e:
        die(f"Failed to write inputs file '{out_filename}': {e}")

    print(f"WROTE: {out_filename}")


def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate MPM particle file and inputs_axialbar.in")
    p.add_argument("no_of_cell_in_x", type=int)
    p.add_argument("buffery", type=int)
    p.add_argument("periodic", type=int, choices=[0, 1])
    p.add_argument("np_per_cell_x", type=int)
    p.add_argument("order_scheme", type=int)
    p.add_argument("alpha_pic_flip", type=float)
    p.add_argument("stress_update_scheme", type=int)
    p.add_argument("CFL", type=float)    
    p.add_argument("output_tag", type=str)
    p.add_argument("--debug", action="store_true", help="Print parsed args and exit")
    return p.parse_args()


def main() -> None:
    args = parse_cli()
    if args.debug:
        print("DEBUG ARGS:", args)

    npart, dx1 = generate_particles_and_return(
        ncells_x=args.no_of_cell_in_x,
        buffery=args.buffery,
        periodic=args.periodic,
        np_per_cell_x=args.np_per_cell_x,
        order_scheme=args.order_scheme,
        alpha_pic_flip=args.alpha_pic_flip,
        stress_update_scheme=args.stress_update_scheme,
        CFL=args.CFL,
        output_tag=args.output_tag,
        out_particles="mpm_particles.dat"
    )

    write_inputs_file(
        ncells_x=args.no_of_cell_in_x,
        buffery=args.buffery,
        periodic=args.periodic,
        np_per_cell_x=args.np_per_cell_x,
        order_scheme=args.order_scheme,
        alpha_pic_flip=args.alpha_pic_flip,
        stress_update_scheme=args.stress_update_scheme,
        CFL=args.CFL,        
        output_tag=args.output_tag,
        dx1=dx1,
        out_filename="inputs_axialbar.dat"
    )

    print("All done.")


if __name__ == "__main__":
    main()

