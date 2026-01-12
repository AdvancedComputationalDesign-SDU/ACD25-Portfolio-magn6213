"""
Assignment 4: Agent-Based Model for Surface Panelization

Author: Magnus Støvring

GhPython Script: Heightmap Surface
Purpose:
    Generate a Perlin-like noise heightmap and build a Rhino NURBS surface
    from the resulting grid of points.

Inputs:
    Width   : int, number of samples in X direction
    Height  : int, number of samples in Y direction
    Scale   : int, roughness scale (larger = smoother)
    Seed    : int, random seed
    Amp     : float, amplitude scaling for Z heights
Outputs:
    Surface : Rhino NurbsSurface built from heightmap
    Pts     : list of Point3d grid points
"""

import Rhino.Geometry as rg
import numpy as np

# --- Noise function ---
def perlin_noise_vectorized(height, width, scale=30, seed=0):
    rng = np.random.default_rng(seed)

    grid_y = height // scale + 2
    grid_x = width // scale + 2
    grid = rng.random((grid_y, grid_x))

    y = np.arange(height)
    x = np.arange(width)
    X, Y = np.meshgrid(x, y)

    gx = X // scale
    gy = Y // scale
    tx = (X % scale) / scale
    ty = (Y % scale) / scale

    v00 = grid[gy, gx]
    v10 = grid[gy, gx + 1]
    v01 = grid[gy + 1, gx]
    v11 = grid[gy + 1, gx + 1]

    i1 = v00 + tx * (v10 - v00)
    i2 = v01 + tx * (v11 - v01)
    noise = i1 + ty * (i2 - i1)

    return noise

# --- Main program ---
if Width < 2 or Height < 2:
    Surface = None
    Pts = []
else:
    noise = perlin_noise_vectorized(Height, Width, scale=Scale, seed=Seed)

    # Normalize to [0,1]
    noise = (noise - noise.min()) / (noise.max() - noise.min())

    # Scale by amplitude
    H = Amp * noise

    # Build grid of points
    pt_grid = []
    for i in range(Width):
        row = []
        for j in range(Height):
            x = float(i)
            y = float(j)
            z = float(H[j, i])  # note: meshgrid indexing
            row.append(rg.Point3d(x, y, z))
        pt_grid.append(row)

    # Flatten points for surface constructor
    flat_pts = [p for row in pt_grid for p in row]

    # Create NURBS surface from point grid
    try:
        Surface = rg.NurbsSurface.CreateFromPoints(flat_pts, Width, Height, 3, 3)
    except:
        Surface = None

    Pts = flat_pts
