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
This design has a wavy, organic surface shaped by Perlin noise, giving it a natural, terrain‑like feel. The tessellation is done with Voronoi cells, which break the surface into irregular, puzzle‑like patches. Vertical branching supports rise up underneath, connecting to the canopy in a way that feels almost like roots meeting a forest floor.

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
### Results
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




---

## Design Variation 2
![alt text](images/myarchitectai_zcqz1jisk_sd.jpg)
Here the canopy surface is defined by a radial sinus function, creating a rhythmic, circular ripple pattern. The tessellation is a simple rectangular grid, which contrasts with the flowing surface by adding order and clarity. The same branching supports grow upward, but against this more regular grid they feel like organic elements weaving into a structured framework.


### Pseudo-Code
```
## 1. Initialize Environment & Parameters
- Import required libraries:
  - Rhino.Geometry as rg
  - NumPy
  - random, heapq
  - scriptcontext as sc
- Set global random seed (`Seed`)
- Define surface generation parameters:
  - Resolution (`ResU`, `ResV`)
  - Radial heightmap parameters (`Amp`, `Frequency`, `Height`)
- Define tessellation and point selection parameters:
  - Grid size (`N`)
  - Minimum spacing (`MinDist`)
  - Pipe radius (`PipeRadius`)
  - Number of selected anchor points (`SelectedCount`)
- Define branching support parameters:
  - Levels (`Levels`)
  - Branch step size (`Step`)
  - Base branch thickness (`Thickness`)
  - Trunk initial height (`up_first`)
  - Noise magnitude (`noise`)
  - Enable noise (`use_noise`)

---

## 2. Radial Heightmap Surface Generation
### 2.1 Define radial heightmap function
- Create normalized grid from -1 to 1 for U and V directions
- Compute radial distance from center
- Apply sinusoidal function: `Z = amplitude * sin(frequency * R * 2π)`
- Normalize Z values to [0,1] and add height offset
- Return heightmap array `Z`

### 2.2 Generate heightmap
- Call `radial_heightmap(ResU, ResV, Amp, Frequency, Height)`

### 2.3 Build point grid
- Initialize empty 2D list `pt_grid`
- For each `(i, j)` in resolution:
  - Assign x, y from indices
  - Assign z from heightmap
  - Create `Point3d(x, y, z)`
  - Append to row
- Flatten 2D grid into `flat_pts`

### 2.4 Create NURBS surface
- Construct surface from `flat_pts` using degree 3 in U and V
- Output surface `Surface`

---

## 3. Safe Rectangular Tessellation & Anchor Point Selection
### 3.1 Tessellate surface grid
- Retrieve U and V domains of surface
- Limit subdivisions: `Nu = min(max(2,N),40)` and `Nv = min(max(2,N),40)`
- Generate `u_vals` and `v_vals` lists
- Evaluate surface points at all `(u,v)` pairs
- Flatten 2D grid into `Sites3d`

### 3.2 Build rectangular Voronoi edges
- For each cell in the grid:
  - Identify four corner points
  - Create four edges (polyline) per cell
  - Convert edges to NURBS curves
- Store all curves in `VoronoiCurves3d`

### 3.3 Pipe Voronoi curves
- For each valid curve:
  - Create Brep pipe with `PipeRadius`
  - Respect document tolerance and angle tolerances
- Store resulting Breps in `VoronoiPipeMeshes`

### 3.4 Select anchor points based on nearest neighbors
- Compute sum of distances to `k=10` nearest neighbors for each site
- Sort points by distance sum
- Select `SelectedCount` points with minimum spacing `MinDist`
- Output selected anchor points `SelectedPoints3d`

---

## 4. Unique Point Extraction & Local Sampling
### 4.1 Extract unique points from curves
- Convert curves to polylines
- Assign spatial hash keys to points using tolerance `tol`
- Keep only unique points in `UniquePoints`

### 4.2 Random closest-point selection
- For each selected anchor point:
  - Compute distances to all `UniquePoints`
  - Find closest 30 points
  - Randomly select 6 points from closest
- Output `SelectedPoints` for branching system

---

## 5. Branching Vertical Support Structures
### 5.1 Define root points
- Specify normalized UV positions for roots
- Map UVs to 3D surface points
- Project roots vertically to ground plane (Z=0)

### 5.2 Prepare canopy geometry
- Convert surface to mesh for line intersection
- Store mesh as `Canopy`

### 5.3 Assign attractors to roots
- For each attractor:
  - Find nearest root
  - Append attractor to that root’s list
- Output dictionary `attr_map` mapping roots → attractors

### 5.4 Compute branch growth direction
- Compute vector toward attractor cluster
- Normalize vector
- Apply height-dependent random noise if enabled
- Return final growth vector

### 5.5 Compute branch radius and create pipe
- Interpolate radius linearly between base `Thickness` and tip `0.05` based on recursion level
- Create cylinder Brep from line segment and radius
- Return Brep

### 5.6 Recursive binary branching
- Base case: stop if no attractors or level = 0
- If last level:
  - Snap directly to attractor and create pipe
- If one attractor:
  - Grow single branch toward cluster
- If multiple attractors:
  - Split into two groups by maximum separation
  - Recursively grow each branch
- Trim all lines against canopy mesh

### 5.7 Grow tree from root points
- Grow vertical trunk from root
- If trunk hits canopy → stop
- Otherwise:
  - Begin recursive binary branching

### 5.8 Main execution
- Assign attractors to roots (`attr_map`)
- For each root:
  - Call `grow_tree(root, attr_map[root], Levels)`
- Store all branch Breps in `PipeBreps`

---

## 6. Outputs
- Original surface (`Surface`)
- Voronoi pipe mesh (`VoronoiPipeMeshes`)
- Branching vertical support pipes (`PipeBreps`)

```
---

### Technical explanation
```
Procedural Surface Generation (Radial Heightmap Canopy)

The first stage generates the base canopy surface using a radial heightmap. The heightmap is defined as a sinusoidal function of the radial distance from the center of the UV domain, producing a smooth scalar field sampled on a regular UV grid.

A normalized 2D grid is constructed with ResU × ResV resolution.

Each grid point is evaluated with a sinusoidal radial function Z = amplitude * sin(frequency * R * 2π).

The resulting height values are normalized to [0,1] and offset by a user-defined height parameter.

The heightmap values are mapped onto a regular XY point lattice, where the Z-coordinate corresponds to the computed height.

This structured point lattice forms the basis for constructing a degree-3 NURBS surface, defining the global canopy geometry and providing a spatial framework for all subsequent stages.

The surface produced in this stage is deterministic given the input parameters and serves as the domain for tessellation, node selection, and branching support growth.

Curvature-Inspired Voronoi Structural Network

The second stage extracts a structural network from the canopy surface using tessellation and local density-based point selection.

Grid-Based Site Sampling

The surface UV domain is subdivided into a rectangular grid of Nu × Nv points.

Grid points are evaluated on the 3D surface to form candidate structural sites (Sites3d).

A nearest-neighbor distance summation metric is computed for each site, representing local density.

Sites with lower total distances are prioritized, enforcing minimum spacing constraints to ensure structural distribution without overcrowding.

A subset of points (SelectedPoints3d) is chosen as anchor nodes for the subsequent branching supports.

Rectangular Tessellation and Voronoi Curves

The UV grid is used to define rectangular cells, each producing four edges.

Edges are converted to polylines and then to NURBS curves, forming the preliminary Voronoi-like structural network.

Pipes are generated along these curves using a fixed radius (PipeRadius), creating a continuous 3D structural mesh.

Voronoi vertices are deduplicated via spatial hashing to identify unique node points suitable for attractor selection.

Localized Attractor Selection

For each anchor point, the 30 nearest unique Voronoi vertices are identified.

A randomized subset of these nearby vertices is selected to form attractor points, which drive the branching vertical support growth in the next stage.

This ensures the attractors are geometrically informed by both local network density and spatial proximity, producing a natural distribution for vertical supports.

Branching Vertical Support Growth

The third stage generates branching vertical supports from ground plane roots to the canopy, guided by attractor points derived from the Voronoi network.

Root Placement and Canopy Intersection

Roots are defined by projecting normalized UV positions of selected points vertically to the ground plane (Z=0).

The canopy surface is converted to a mesh to facilitate line-mesh intersection computations.

Initial growth begins with a vertical trunk extending upward from each root.

Branch growth is trimmed upon canopy intersection, ensuring a precise connection between supports and surface geometry.

Attractor-Based Recursive Binary Branching

Each attractor is assigned to the nearest root, forming localized growth domains.

The branching system grows recursively using a binary splitting strategy:

Attractors are partitioned into two clusters based on maximum pairwise separation.

Branch directions are computed as normalized vectors toward cluster centroids.

Height-dependent random noise is applied to branch vectors, increasing directional variation toward higher recursion levels.

Branches extend incrementally, generating cylindrical Breps along each segment.

At the final recursion level, branches snap directly to attractor points, guaranteeing geometric connectivity.

Branch Thickness and Structural Logic

Branch radius decreases linearly with recursion depth: thickest at the trunk base and tapering toward terminal branches.

Each branch segment is aligned along its computed growth vector, forming smooth, structurally coherent cylinders.

The final output is a hierarchical branching support system, visually and functionally integrated with the canopy surface.
```

---
### Results
<table>
  <tr>
    <td align="center">
      <img src="images/V2-Top.jpg" width="400"><br>
      <strong>Variation 2 Top View</strong>
    </td>
    <td align="center">
      <img src="images/V2-Front.jpg" width="400"><br>
      <strong>Variation 2 Front View</strong>
    </td>
  </tr>
</table>

<table>
  <tr>
    <td align="center">
      <img src="images/V2-P.jpg" width="400"><br>
      <strong>Variation 2 Perspective View</strong>
    </td>
    <td align="center">
      <img src="images/V2-2P.jpg" width="400"><br>
      <strong>Variation 2 Two Point Perspective View</strong>
    </td>
  </tr>
</table>

---

## Design Variation 3
![alt text](images/myarchitectai_lxt60yx8e_sd.jpg)
This variation uses scattered Gaussian bumps to create a surface with soft, hill‑like protrusions. The tessellation is a diagrid, forming a network of diagonal lines that emphasize geometry and stability. The branching supports rise into this canopy, reinforcing the sense of a lightweight but intricate lattice structure.

### Pseudo-Code
```
1. Initialize Environment & Parameters

Import required libraries:

Rhino.Geometry as rg  
rhinoscriptsyntax as rs  
numpy as np  
random, heapq, math  
scriptcontext as sc  

Set global random seed (Seed)

Define Gaussian bump parameters:

Resolution (ResU, ResV)  
Amplitude (Amp)  
Number of bump centers (NumCenters)  
Spread (Spread)  
Height offset (Height)  

Define tessellation and point selection parameters:

Grid divisions (N)  
Minimum spacing (MinSpacing)  
Pipe radius (PipeRadius)  
Number of selected anchor points (SelectedCount)  

Define branching support parameters:

Levels (Levels)  
Step size (Step)  
Base branch thickness (Thickness)  
Trunk initial height (up_first)  
Noise magnitude (noise)  
Enable noise (use_noise)

2. Gaussian Bump Surface Generation

2.1 Define Gaussian bump heightmap function

Initialize random number generator with seed  
Generate random 2D bump centers in normalized [0,1] domain  
Create normalized meshgrid (xx, yy) for U and V  

For each bump center (cx, cy):  
    Compute squared distance dist2 = (xx - cx)^2 + (yy - cy)^2  
    Add Gaussian: Z += exp(-dist2 / spread)  

Normalize Z to [0,1], scale by amplitude, and add height offset  
Return heightmap array Z  

2.2 Generate heightmap

Call gaussian_bumps(ResU, ResV, Amp, NumCenters, Spread, Seed, Height)

2.3 Build point grid

Initialize empty 2D list pt_grid  

For each (i, j) in resolution:  
    Assign x, y from indices  
    Assign z from heightmap  
    Create Point3d(x, y, z)  
    Append to row  

Flatten 2D grid into flat_pts  

2.4 Create NURBS surface

Construct NURBS surface from flat_pts using degree 3 in U and V  
Output surface Surface  

3. Triangular / Diagrid Tessellation & Anchor Point Selection

3.1 Build surface grid

Retrieve surface U and V domains  
Limit subdivisions: Nu = max(2, N), Nv = max(2, N)  
Generate u_vals and v_vals lists  
Evaluate surface points at all (u, v) pairs  
Store 2D grid_points and flatten to Sites3d  

3.2 Create base grid edges

For each row in grid_points:  
    Create polyline along row → convert to NURBS curve  

For each column in grid_points:  
    Create polyline along column → convert to NURBS curve  

Store all curves in TriCurves3d  

3.3 Add diagonals for triangles/diagrid

For each quad (i, j) in grid:  
    Diagonal 1: p00 → p11 → polyline → NURBS  
    Diagonal 2: p10 → p01 → polyline → NURBS  
    Append to TriCurves3d  

3.4 Pipe all curves

For each valid curve in TriCurves3d:  
    Create Brep pipe using PipeRadius and document tolerances  
    Append to VoronoiPipeMeshes  

Output curves as VoronoiCurves3d  

3.5 Select anchor points based on nearest neighbors

For each site, compute sum of distances to k=10 nearest neighbors  
Sort points by distance sum  
Select SelectedCount points ensuring minimum spacing MinSpacing  
Output SelectedPoints3d  

4. Unique Point Extraction & Random Closest-Point Sampling

4.1 Extract unique points from curves

Convert curves to polylines  
Compute spatial hash keys with tolerance tol  
Keep only unique points → UniquePoints  

4.2 Random selection of closest points

For each selected anchor point:  
    Compute distances to all UniquePoints  
    Find closest 30 points  
    Randomly select 6 points  

Output SelectedPoints for branching system  

5. Branching Vertical Support Structures

5.1 Define root points

Specify normalized UV positions for roots  
Map UVs to 3D surface points  
Project roots vertically to ground plane (Z=0) → roots  

5.2 Prepare canopy geometry

Convert surface to mesh → Canopy  

5.3 Attractor assignment to roots

For each attractor point:  
    Find nearest root  
    Append attractor to that root’s list  

Output dictionary attr_map mapping roots → attractors  

5.4 Compute cluster growth direction

For given start_pt and attractor cluster:  
    Compute vector toward attractor cluster  
    Normalize vector  
    Apply height-dependent random noise if use_noise=True  

Return final growth vector  

5.5 Branch radius & pipe creation

Compute radius linearly from Thickness at base to 0.05 at tip  
Create cylinder Brep along line segment with computed radius  
Return Brep  

5.6 Trim line to canopy

Check intersection of line with mesh canopy  
If intersection exists → shorten line at hit point  
Return trimmed endpoint and hit status  

5.7 Recursive binary branching

Base case: stop if no attractors or level = 0  
Level 1: snap directly to single attractor and create pipe  

Single attractor:  
    Grow branch toward cluster  

Multiple attractors:  
    Split into two groups by maximum separation  
    Recursively grow each branch  

Trim all lines against canopy mesh  

5.8 Grow tree from root points

Grow vertical trunk from root to up_first height  
Trim trunk against canopy  

If trunk hits canopy → stop  
Otherwise → begin recursive binary branching  

Store all branch Breps in PipeBreps  

5.9 Main execution

Assign attractors to roots → attr_map  

For each root:  
    Call grow_tree(root, attr_map[root], Levels)  

6. Outputs

Original NURBS surface → Surface  
Triangular/diagrid pipe mesh → VoronoiPipeMeshes  
Branching vertical support pipes → PipeBreps

```
---

### Technical explanation
```
Three parameter groups define the system’s behaviour:

Gaussian bump surface parameters control the resolution and morphology of the heightmap.

Tessellation parameters govern the density and spacing of the diagrid structure.

Branching parameters define the recursive growth behaviour of the vertical support system.

2. Gaussian Bump Surface Generation
A procedural heightmap is created by superimposing Gaussian functions across a normalized domain. Randomly distributed centres define the bumps, and each contributes a smooth radial elevation. The combined field is normalized, scaled, and offset to produce a continuous terrain. This scalar field is sampled into a grid of points, which are then reconstructed into a degree‑3 NURBS surface. The result is a smooth, differentiable surface suitable for tessellation and structural analysis.

3. Triangular / Diagrid Tessellation & Anchor Point Selection
The surface domain is subdivided into a grid of UV coordinates, yielding a structured point set. From this, polylines are generated along rows and columns, then diagonals are added to form a triangular or diagrid pattern. These curves are thickened into Brep pipes, producing a tessellated canopy structure. Anchor points are selected by ranking sites based on their proximity to neighbours, ensuring both distribution and minimum spacing. This step identifies structurally meaningful attractors for the branching system.

4. Unique Point Extraction & Random Closest‑Point Sampling
All tessellation curves are converted to polylines, and duplicates are removed using spatial hashing. This yields a clean set of unique points. For each anchor point, the nearest neighbours are identified, and a random subset is selected. This introduces controlled variability into the branching system, preventing deterministic repetition and ensuring organic distribution of attractors.

5. Branching Vertical Support Structures
Roots are defined in UV space, mapped to the surface, and projected vertically to the ground plane. The canopy surface is meshed to allow efficient intersection tests. Attractors are assigned to roots, forming clusters that guide growth. Branches grow recursively toward attractor clusters, with optional noise introducing organic variation. Thickness decreases linearly from trunk to tip, and each segment is piped into a cylindrical Brep. Branches are trimmed against the canopy mesh to ensure structural integration. Recursive binary branching produces bifurcated tree‑like supports, combining deterministic geometry with stochastic variation.

6. Outputs
The system produces three primary outputs:

The Gaussian bump NURBS surface.

The tessellated diagrid pipe mesh.

The branching vertical support structures.

Together, these outputs form a complete architectural system that integrates terrain generation, tessellation, and biologically inspired branching supports.
```

---
### Results
<table>
  <tr>
    <td align="center">
      <img src="images/V3-Top.jpg" width="400"><br>
      <strong>Variation 3 Top View</strong>
    </td>
    <td align="center">
      <img src="images/V3-Front.jpg" width="400"><br>
      <strong>Variation 3 Front View</strong>
    </td>
  </tr>
</table>

<table>
  <tr>
    <td align="center">
      <img src="images/V3-P.jpg" width="400"><br>
      <strong>Variation 3 Perspective View</strong>
    </td>
    <td align="center">
      <img src="images/V3-2P.jpg" width="400"><br>
      <strong>Variation 3 Two Point Perspective View</strong>
    </td>
  </tr>
</table>


## Surfaces
The three canopy surfaces demonstrate distinct parametric approaches to surface modulation. The Perlin noise surface is generated through a stochastic gradient‑based function, producing irregular but continuous undulations across the domain. This method highlights the capacity of noise functions to introduce controlled randomness while maintaining overall continuity. The radial sinus surface applies a periodic function in polar coordinates, resulting in concentric oscillations that emphasize symmetry and regularity. This approach illustrates how trigonometric functions can be used to impose ordered variation on a surface. The Gaussian bump surface is constructed by superimposing multiple Gaussian kernels distributed across the domain, creating localized height variations with smooth decay. This technique demonstrates how kernel‑based functions can generate surfaces with discrete, controllable features.
![alt text](images/V1-Surface.jpg)
![alt text](images/V2-Surface.jpg)
![alt text](images/V3-Surface.jpg)