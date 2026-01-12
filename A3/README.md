---
layout: default
title: Project Documentation
parent: "A3: Parametric Structural Canopy"
nav_order: 2
nav_exclude: false
search_exclude: false
---

# Assignment 3: Parametric Structural Canopy

## Design Variation 1
![alt text](images/myarchitectai_3du3dwqfm_sd.jpg)

### Pseudo-Code
```
## 1. Initialize Environment & Parameters
- Import required libraries:
  - Rhino.Geometry
  - rhinoscriptsyntax
  - scriptcontext
  - NumPy
  - random, math, heapq
- Set global random seed (`Seed`)
- Define global parameters:
  - Surface resolution (`ResU`, `ResV`)
  - Noise parameters (`Scale`, `Amp`, `Height`)
  - Tessellation parameters (`N`, `MinDist`, `PipeRadius`)
  - Point selection parameters (`SelectedCount`, `MinSpacing`)
  - Branching parameters (`Levels`, `Step`, `Thickness`)

---

## 2. Perlin Noise Surface Generation

### 2.1 Define Perlin Noise Function
- Initialize NumPy random generator using `Seed`
- Generate random gradient vectors at grid points
- Create scaled coordinate grid
- Compute integer cell corner indices
- Compute local coordinates within each cell
- Define smooth fade function
- Compute interpolation weights
- Evaluate gradient dot products at four corners
- Perform bilinear interpolation
- Normalize noise values to `[0, 1]`
- Return Perlin noise heightmap

---

### 2.2 Generate Heightmap
- Call Perlin noise function with `(ResU, ResV, Scale, Seed)`
- Apply amplitude scaling
- Apply vertical height offset

---

### 2.3 Build Point Grid
- Initialize empty 2D grid
- For each `(i, j)` in resolution:
  - Assign `(x, y)` from grid indices
  - Assign `z` from heightmap
  - Create `Point3d(x, y, z)`
- Flatten grid into ordered point list

---

### 2.4 Create Surface
- Construct NURBS surface from flattened point list
- Use degree 3 in both U and V directions
- Output surface (`Surface`)

---

## 3. Curvature-Weighted Voronoi Tessellation

---

### 3.1 Generate Curvature-Weighted Sites
- Sample surface UV domain on a regular grid
- Evaluate surface mean curvature at each UV
- Use curvature magnitude as sampling weight
- Normalize weights
- Randomly select `N` sites:
  - Enforce minimum spatial distance constraint
- Store selected UV coordinates

---

### 3.2 Convert Sites to 3D
- Map UV site coordinates to 3D surface points
- Store as `Sites3d`

---

### 3.3 Construct Voronoi Cells in UV Space
- Define UV bounding box
- For each site:
  - Initialize polygon as bounding box
  - For each other site:
    - Compute perpendicular bisector
    - Clip polygon using half-plane test
  - Store resulting cell polygon

---

### 3.4 Map Voronoi Cells to 3D Geometry
- Convert UV polygon vertices to surface points
- Close each polygon
- Convert to NURBS curves
- Store Voronoi boundary curves

---

### 3.5 Pipe Voronoi Curves
- For each valid Voronoi curve:
  - Create circular pipe with given radius
  - Use document tolerances
- Store resulting Brep pipes

---

### 3.6 Select Structural Anchor Points
- Compute k-nearest neighbor distance sum for each site
- Sort sites by total distance
- Select evenly distributed subset with minimum spacing
- Output selected points (`SelectedPoints3d`)

---

## 4. Unique Point Extraction & Local Sampling

### 4.1 Extract Unique Curve Points
- Convert Voronoi curves to polylines
- Bucket points using spatial hashing
- Remove duplicates within tolerance
- Store unique points

---

### 4.2 Random Closest-Point Selection
- For each selected site:
  - Compute distances to all unique points
  - Find closest candidates
  - Randomly select fixed number of points
- Output attractor points for branching system

---

## 5. Branching Vertical Support Structure

---

### 5.1 Define Root Points
- Define normalized UV root positions
- Map UV roots to surface
- Project roots vertically to ground plane

---

### 5.2 Prepare Canopy Geometry
- Convert surface to mesh
- Use mesh for line-intersection trimming

---

### 5.3 Assign Attractors to Roots
- For each attractor:
  - Assign to nearest root
- Group attractors per root

---

### 5.4 Branch Direction Computation
- Compute vector toward attractor cluster
- Normalize direction
- Apply height-dependent random noise
- Return final growth direction

---

### 5.5 Pipe Geometry Creation
- Compute branch radius based on recursion level
- Convert line segments into cylindrical Breps
- Store pipe geometry

---

### 5.6 Recursive Binary Branching
- If no attractors or level exhausted → stop
- If final level:
  - Connect directly to attractor
- If one attractor:
  - Grow single branch toward cluster
- If multiple attractors:
  - Split into two groups by maximum separation
  - Recursively grow branches for each group
- Trim branches against canopy mesh

---

### 5.7 Tree Growth
- Grow vertical trunk from root
- If trunk hits canopy → stop
- Otherwise:
  - Start recursive binary branching

---

## 6. Final Execution
- Assign attractors to root points
- Grow branching structure from each root
- Output:
  - Original surface
  - Voronoi pipe mesh
  - Branching support structure Breps
```
---

### Technical explanation
```
The system is structured as a multi-stage parametric pipeline for generating a structurally differentiated canopy driven by surface geometry, curvature analysis, and branching support logic. The workflow is composed of three connected GhPython scripts, each operating at a different level of abstraction but sharing a common geometric domain.

Procedural Surface Generation (Perlin Noise Canopy)

The first stage generates the base canopy surface using a procedural heightmap derived from Perlin noise. A custom 2D Perlin noise function produces a smooth scalar field sampled across a regular UV grid.

Random gradient vectors are assigned per grid node using a seeded random number generator, ensuring repeatability. Noise values are interpolated using a quintic fade function to maintain smooth transitions between grid cells. The resulting noise field is normalized, scaled, and offset to control surface amplitude and base height.

The heightmap is mapped onto a regular XY point grid, where each point’s Z-value corresponds to the noise-derived height. This structured point lattice is used to construct a degree-3 NURBS surface. The resulting surface defines the global geometry of the canopy and acts as the spatial and analytical foundation for all subsequent computations.

Curvature-Weighted Voronoi Structural Network

The second stage extracts a structural network directly from the canopy surface using curvature-driven sampling and explicit Voronoi construction in parameter space.

2.1 Curvature-Weighted Site Generation

Candidate UV locations are sampled across the surface domain at a fixed resolution. At each sample point, the surface mean curvature is evaluated. The absolute curvature value is used as a probability weight, biasing site selection toward regions of higher geometric intensity.

A stochastic sampling process selects a specified number of sites while enforcing a minimum 3D distance constraint between them. This prevents clustering and ensures an even yet curvature-aware distribution of structural nodes across the surface.

2.2 Voronoi Diagram Construction in UV Space

A Voronoi diagram is constructed manually in UV space using geometric half-plane clipping rather than relying on built-in Voronoi solvers.

For each site, a bounding polygon corresponding to the surface UV domain is iteratively clipped against perpendicular bisectors formed between the site and every other site. Each clipping operation preserves only the half-plane closer to the current site, resulting in a convex Voronoi cell.

Once computed, the UV-space Voronoi polygons are mapped back onto the 3D surface by evaluating surface points at their UV coordinates.

2.3 Structural Curve and Pipe Generation

The edges of each Voronoi cell are converted into closed 3D polylines and rebuilt as NURBS curves. These curves are thickened into pipe breps using Rhino’s pipe construction, producing a continuous structural network that follows the curvature-informed surface partitioning.

2.4 Node Selection for Structural Supports

Surface points corresponding to Voronoi vertices are collected and deduplicated using spatial hashing. A subset of points is selected based on nearest-neighbor distance summation, favoring points located in denser regions of the Voronoi network.

For each selected point, a randomized subset of nearby Voronoi vertices is chosen. These points serve as attractors for the vertical branching support system in the next stage.

Branching Vertical Support Growth

The third stage generates branching vertical supports that connect the ground plane to the canopy, guided by attractor points derived from the Voronoi network.

3.1 Root Placement and Canopy Intersection

Root points are defined by projecting selected surface UV locations vertically down to the ground plane. The canopy surface is meshed to allow efficient ray-based intersection tests.

Each support begins with a vertical trunk that grows upward from the root point. Growth is trimmed automatically upon intersection with the canopy mesh, ensuring precise structural connection.

3.2 Attractor-Based Binary Branching

Attractor points are assigned to the nearest root, forming localized growth domains. Each support grows recursively using a binary branching algorithm:

Attractors are split into two spatial clusters based on maximum separation.

Branch directions are computed as normalized vectors toward attractor clusters.

Height-dependent noise is applied to branch direction vectors, increasing variation toward upper levels.

Branch segments grow incrementally and are trimmed when intersecting the canopy. At the final recursion level, branches snap directly to attractor points, ensuring geometric connectivity.

3.3 Structural Thickness Logic

Branch thickness decreases linearly with recursion depth. Trunks are thickest at the base, while terminal branches taper toward a minimum radius. Each segment is generated as a cylindrical brep aligned to the growth direction.
```

---

<table>
  <tr>
    <td align="center">
      <img src="images/V1-Top.jpg" width="400"><br>
      <strong>Variation 1 Top View</strong>
    </td>
    <td align="center">
      <img src="images/V1-Front.jpg" width="400"><br>
      <strong>Variation 1 Front View</strong>
    </td>
  </tr>
</table>

<table>
  <tr>
    <td align="center">
      <img src="images/V1-P.jpg" width="400"><br>
      <strong>Variation 1 Perspective View</strong>
    </td>
    <td align="center">
      <img src="images/V1-2P.jpg" width="400"><br>
      <strong>Variation 1 Two Point Perspective View</strong>
    </td>
  </tr>
</table>