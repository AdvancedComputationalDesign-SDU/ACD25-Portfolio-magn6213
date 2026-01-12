## Magnus Støvring - ACD-E25 Portfolio

## Overview
This portfolio documents four computational design studies developed for Advanced Computational Design. Across the series, I investigate emergent geometric organisation through algorithmic modelling, recursive logic, field‑driven geometry, and agent‑based systems, iterating through prototypes, parameter sweeps, and visual evaluation. The work examines how local rules, scalar fields, and computational constraints can be composed into coherent spatial outcomes.

## Assignments

### A1: NumPy Array Manipulation for 2D Pattern Generation
![alt text](A1/images/perlin_output_seed_5_Scale_5_Colormap.png)
In A1, I construct a Perlin‑like noise generator using fully vectorized NumPy operations. The system builds a coarse random grid, interpolates across it bilinearly, and normalizes the result into a continuous heightmap. This heightmap becomes the basis for pixel‑level pattern generation.

Two visualization modes are implemented:

- Manual RGB mode, where the noise drives a single colour channel for controlled, minimalistic outputs.

- Colormap mode, where the noise is mapped through magma_r to produce terrain‑like gradients.

The study focuses on how array logic, interpolation, and controlled randomness can produce structured yet expressive image families. The toggle‑ready code structure supports rapid iteration across seeds, scales, and colour mappings, enabling systematic exploration of pattern behaviour.

### A2: Exploring Fractals through Recursive Geometric Patterns
![alt text](A2/images/DLA_seed42_p1000_att500_a1-1.5.png)
A2 investigates diffusion‑limited aggregation (DLA) through a recursive particle‑walking algorithm. Particles spawn near the boundary of a circular domain and perform random walks until they either attach to the growing cluster or exceed a recursion limit.

Key features of the system include:

- Recursive walker function that encodes the entire diffusion process as a self‑calling rule.

- Adaptive attachment radius, modulated by distance to an attractor point, producing subtle directional bias.

- Dynamic spawn radius, expanding as the cluster grows.

- Line‑based visualization, where each new particle connects to its nearest neighbour with a colour and thickness gradient.

The resulting structure exhibits classic fractal branching behaviour, with local randomness accumulating into global hierarchy. The recursion‑based implementation highlights how simple rules can generate complex, organic geometries.

### A3: Parametric Structural Canopy
![alt text](A3/images/myarchitectai_3du3dwqfm_sd.jpg)
In A3, I design a series of canopy systems driven by different scalar fields and tessellation strategies, implemented in Grasshopper and GhPython. Each variation explores how surface logic, panelization, and structural support can be composed into coherent spatial assemblies.

#### Variation 1 — Perlin‑Noise Terrain + Voronoi Tessellation
A wavy, organic surface generated from Perlin noise produces a terrain‑like canopy. Voronoi tessellation breaks the surface into irregular, puzzle‑like patches. Vertical branching supports rise from below, connecting to the canopy like roots meeting a forest floor, reinforcing the naturalistic character of the geometry.

#### Variation 2 — Radial Sinus Surface + Rectangular Grid
A radial sinus function creates rhythmic circular ripples across the canopy. A simple rectangular grid tessellation overlays this surface, introducing order and clarity. The branching supports weave upward through this regular framework, creating a dialogue between organic growth and geometric structure.

#### Variation 3 — Gaussian Bumps + Diagrid
Scattered Gaussian bumps generate soft hill‑like protrusions. A diagrid tessellation emphasizes stability and directional flow, producing a lightweight yet intricate lattice. The branching supports integrate into this network, reinforcing the sense of a geometric shell supported by organic struts.

Across all three variations, the study examines how field‑driven surfaces, tessellation logics, and structural behaviours interact to produce coherent architectural systems.

### A4: Agent-Based Modeling for Surface Panelization
![alt text](A4/images/Banner.gif)
A4 develops an agent‑based panelization system in Grasshopper using three custom GhPython components: a surface builder, an agent builder, and an agent simulator. The simulation operates on a Rhino surface derived from my A1 heightmap, establishing continuity across the course.

Agents are initialized across the surface and move through time according to multiple geometric signals:

Curvature‑Driven Behaviour
Agents compute the gradient of absolute Gaussian curvature and move toward regions where curvature changes most rapidly. This produces flow lines that align with geometric features of the surface.

Spatial Influences
- Separation: agents avoid neighbours to maintain spacing.

- Repulsion: agents steer away from user‑defined repulsor points.

- Boundary proximity: agents slow or freeze near UV boundaries to prevent escape.

Agent states persist across iterations using scriptcontext.sticky, enabling interactive, real‑time simulation.

The final panelization is generated by feeding agent positions into Grasshopper’s Delaunay Mesh component. The resulting topology reflects the emergent distribution of agents, producing a panel system shaped not by a single global rule but by the interplay of multiple local behaviours.

