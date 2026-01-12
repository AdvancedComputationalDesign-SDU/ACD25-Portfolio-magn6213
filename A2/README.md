---
layout: default
title: Project Documentation
parent: "A2: Exploring Fractals through Recursive Geometric Patterns"
nav_order: 2
nav_exclude: false
search_exclude: false
---

# Assignment 2: Exploring Fractals through Recursive Geometric Patterns

[View on GitHub]({{ site.github.repository_url }})

![Example Fractal](images/branching.png)

## Table of Contents

- [Pseudo-Code](#pseudo-code)
- [Technical Explanation](#technical-explanation)
- [Geometric Influences](#geometric-influences)
- [Appearance Mapping](#appearance-mapping)
- [Experiments](#experiments)
- [Challenges and Solutions](#challenges-and-solutions)
- [References](#references)

---

## Pseudo‑Code

### 1. Initialize Simulation Parameters
- Set `SEED` and initialize both Python’s `random` and NumPy RNG.
- Define:
  - `circle_radius`
  - `num_particles`
  - `step_size`
  - `max_attempts`
- Define attractor behaviour:
  - `attractor_point = [circle_radius, circle_radius]`
  - `min_attach`, `max_attach`
- Initialize the cluster with:
  - `aggregated_particles = [(0, 0)]`

### 2. Set Up Plotting Environment
- Create Matplotlib figure with black background.
- Fix axis limits to the circular domain.
- Load colormap (`plasma`).

### 3. Helper Functions
- **is_inside_circle(x, y, radius)**  
  Returns True if `(x, y)` lies inside the circular boundary.

- **random_walk(x, y, step_size)**  
  Moves the particle one step in a random direction.

- **is_near_aggregated(x, y, aggregated_particles, step_size)**  
  - Compute distance to attractor.
  - Convert to closeness factor.
  - Compute dynamic attachment radius:  
    `attach_radius = min_attach + closeness * (max_attach - min_attach)`
  - Check if particle is within attachment threshold of any aggregated particle.

- **compute_spawn_radius(aggregated_particles)**  
  - Find maximum squared distance of the cluster.
  - Spawn new particles slightly outside the current boundary.

### 4. Recursive Random‑Walk Function
**walk_particle(x, y, attempts)**  
- Base cases:
  - If `attempts >= max_attempts`: return None
  - If particle leaves circle: return None
  - If particle is close enough to attach: return `(x, y)`
- Recursive step:
  - Move particle using `random_walk`
  - Call `walk_particle` again with incremented attempt count

### 5. Main Aggregation Loop
For each particle in `num_particles`:
- Compute spawn radius.
- Spawn particle at random angle.
- Call `walk_particle(x, y, 0)`.
- If result is valid:
  - Append to `aggregated_particles`.

### 6. Plotting the Final Structure
- Convert aggregated particles to NumPy array.
- For each particle (except the first):
  - Compute distances to all previous particles.
  - Find nearest neighbour.
  - Draw line segment between them.
  - Line width decreases from `max_width` to `min_width`.
  - Line colour based on distance from origin.

### 7. Measure Runtime
- Record start and end times.
- Compute total runtime.

### 8. Add Legend Text
- Display simulation parameters and runtime on the figure.

### 9. Save and Show Figure
- Save final image with filename encoding:
  - seed  
  - particle count  
  - max attempts  
  - attachment parameters  
- Display the figure.


---

## Technical Explanation

This implementation simulates diffusion‑limited aggregation (DLA) using a recursive random‑walk approach. The simulation begins with a single seed particle at the origin. Each new particle is spawned just outside the current cluster boundary and performs a random walk until it either attaches to the cluster or exceeds the maximum number of allowed steps.

Attachment is influenced by an attractor point located in the upper‑right corner of the domain. The particle’s distance to this attractor determines a dynamic attachment radius: particles closer to the attractor have a larger effective radius and therefore a higher chance of sticking. This introduces directional bias and produces asymmetric growth patterns.

The random walk is implemented recursively. Each step checks three conditions:
1. whether the particle has exceeded the maximum number of attempts,
2. whether it has left the circular simulation boundary,
3. whether it is close enough to attach to the existing cluster.

If none of these conditions are met, the particle takes another random step and the function calls itself again.

Once all particles have been aggregated, the structure is visualized by connecting each particle to its nearest neighbour. Line thickness decreases with particle index, creating thicker “trunk” branches and thinner offshoots. Line colour is determined by distance from the origin using a plasma colormap, highlighting the geometric structure of the cluster.

The result is a fractal‑like branching pattern characteristic of DLA systems. The combination of random walk, recursive movement, and attractor‑biased attachment produces complex emergent behaviour from simple local rules.


---

## Geometric Influences

1. Attractor Point:
- The attractor point is a fixed location in the simulation (in this case the top right corner) toward which particle growth is biased.
- Influence Computation:
  - For each particle, compute its Euclidean distance to the attractor point.
  - Calculate a "closeness" factor: closer particles have a higher probability of attaching.
  - Scale the attachment radius using this closeness:
      attach_radius = min_attach + closeness * (max_attach - min_attach)
- Effect on Growth:
  - Particles near the attractor attach more easily, producing asymmetric, biased growth.
  - This modulates the cluster formation dynamically, directing the fractal toward the attractor while still allowing random branching.

2. Circle Radius (Boundary Limit):
- The circle radius defines the maximum allowed extent of the particle cluster.
- Influence Computation:
  - Check each particle’s position: if it moves outside the circle, terminate the random walk.
  - Acts as a clipping region, preventing growth beyond a set boundary.
- Effect on Growth:
  - Limits the overall size of the fractal, keeping the structure contained within a visually meaningful area.
  - Ensures that the fractal’s density and shape are modulated by a geometric boundary, influencing both local and global growth patterns.

Both influences are applied during the particle’s random walk and attachment check:
- The attractor dynamically modulates attachment probability.
- The circle radius constrains the allowable space for growth.
Together, they produce a fractal pattern that is both structured and self-limiting.

---

## Appearance Mapping

- Line Color:
  - The color of each line segment is mapped to the distance of the particle from the origin
    (or optionally the attractor point).
  - Computation:
    - distance = np.hypot(x, y)  # Euclidean distance from origin
    - color = cmap(distance / circle_radius)
  - Justification:
    - Mapping color to distance highlights the spatial structure of the cluster.
    - Provides a visual gradient that emphasizes the cluster's growth toward the attractor.

- Line Width:
  - The line width decreases as more particles are added to the cluster.
  - Computation:
    - line_width = max_width - (max_width - min_width) * (i / total_particles)
  - Justification:
    - Thicker lines represent early “trunk” branches, while thinner lines show later, smaller offshoots.
    - Mimics natural growth patterns (like plant branching) and enhances depth perception.

- Overall Effect:
  - These mappings turn abstract particle positions into an intuitive visualization.
  - Observers can easily identify main growth directions, branching hierarchy, and the influence of attractor points.
  - Provides both aesthetic appeal and meaningful interpretation of the fractal pattern.

---

## Experiments

### Changing random seed
The random seed, dictates the way the pattern is formed, and thus changing the random seed gives very different patterns.
<table>
  <tr>
    <td align="center">
      <img src="images/DLA_seed40_p1000_att900_a1-1.5.png" width="220"><br>
      <strong>Seed 40</strong>
    </td>
    <td align="center">
      <img src="images/DLA_seed41_p1000_att900_a1-1.5.png" width="220"><br>
      <strong>Seed 41</strong>
    </td>
    <td align="center">
      <img src="images/DLA_seed42_p1000_att900_a1-1.5.png" width="220"><br>
      <strong>Seed 42</strong>
    </td>
    <td align="center">
      <img src="images/DLA_seed43_p1000_att900_a1-1.5.png" width="220"><br>
      <strong>Seed 43</strong>
    </td>
  </tr>
</table>


### Changing number of particles
The number of particles dictates how many points will try to attach as part of the pattern but in this case it can be seen that the difference between 2.500 and 5.000 is minimal.

<table>
  <tr>
    <td align="center">
      <img src="images/DLA_seed42_p500_att900_a1-1.5.png" width="220"><br>
      <strong>500 Particles</strong>
    </td>
    <td align="center">
      <img src="images/DLA_seed42_p1000_att900_a1-1.5.png" width="220"><br>
      <strong>1000 Particles</strong>
    </td>
    <td align="center">
      <img src="images/DLA_seed42_p2500_att900_a1-1.5.png" width="220"><br>
      <strong>2500 Particles</strong>
    </td>
    <td align="center">
      <img src="images/DLA_seed42_p5000_att900_a1-1.5.png" width="220"><br>
      <strong>5000 Particles</strong>
    </td>
  </tr>
</table>



### Changing maximum attempts
The maximum attempts dictates the amount of tries a point has to be attached to the existing structure, so changing the number changes the entire way the pattern is formed.

<table>
  <tr>
    <td align="center">
      <img src="images/DLA_seed42_p1000_att250_a1-1.5.png" width="220"><br>
      <strong>250 Attempts</strong>
    </td>
    <td align="center">
      <img src="images/DLA_seed42_p1000_att500_a1-1.5.png" width="220"><br>
      <strong>500 Attempts</strong>
    </td>
    <td align="center">
      <img src="images/DLA_seed42_p1000_att750_a1-1.5.png" width="220"><br>
      <strong>750 Attempts</strong>
    </td>
    <td align="center">
      <img src="images/DLA_seed42_p1000_att950_a1-1.5.png" width="220"><br>
      <strong>950 Attempts</strong>
    </td>
  </tr>
</table>



### Changing max_attach
Increasing the max_attach represents increasing the attraction force of the attraction point, in this case the upper right corner. As it can be seen this "pushes" the branches more and more towards the corner.

<table>
  <tr>
    <td align="center">
      <img src="images/DLA_seed42_p1000_att900_a1-1.5.png" width="220"><br>
      <strong>max_attach = 1.5</strong>
    </td>
    <td align="center">
      <img src="images/DLA_seed42_p1000_att900_a1-2.png" width="220"><br>
      <strong>max_attach = 2</strong>
    </td>
    <td align="center">
      <img src="images/DLA_seed42_p1000_att900_a1-4.png" width="220"><br>
      <strong>max_attach = 4</strong>
    </td>
    <td align="center">
      <img src="images/DLA_seed42_p1000_att900_a1-6.png" width="220"><br>
      <strong>max_attach = 6</strong>
    </td>
  </tr>
</table>



---

## Challenges and Solutions

- **Challenge**: Recursive random walks can be slow when many particles fail to attach.
  - **Solution**: Introduced clear base‑case checks (max attempts, boundary exit, attachment proximity) to terminate unproductive walks early, reducing unnecessary recursion depth.

- **Challenge**: Balancing attractor‑biased growth with natural diffusion‑limited randomness.
  - **Solution**: Computed a dynamic attachment radius based on distance to the attractor. This preserves stochastic branching while still guiding the cluster toward the attractor.

- **Challenge**: Preventing particles from wandering indefinitely and keeping the cluster within a meaningful region.
  - **Solution**: Enforced a circular boundary using `is_inside_circle()`, immediately terminating particles that leave the domain.

- **Challenge**: Efficiently checking whether a particle is close enough to attach to the cluster.
  - **Solution**: Used NumPy vectorized distance calculations inside `is_near_aggregated()`, avoiding Python loops and enabling fast proximity checks even for large clusters.

- **Challenge**: Visualizing the final structure without visual clutter or uniform line thickness.
  - **Solution**: Connected each particle only to its nearest neighbour, varied line width based on particle index, and mapped color to radial distance, producing a readable and expressive fractal visualization.



---

## References

- **Model of Diffusion**: [https://phas.ubc.ca/~berciu/TEACHING/PHYS349/DLA.pdf](https://phas.ubc.ca/~berciu/TEACHING/PHYS349/DLA.pdf)
- **Information on Diffusion Limited Aggregation**: [https://en.wikipedia.org/wiki/Diffusion-limited_aggregation](https://en.wikipedia.org/wiki/Diffusion-limited_aggregation)
- **Youtube tutorial for inspiration**: [https://www.youtube.com/watch?v=BD6s_gCAiAU](https://www.youtube.com/watch?v=BD6s_gCAiAU)

A combination of ChatGPT and Copilot has been used during the project, primarily for troubleshooting related to the fractal generation aswell as generation of code for the README.md file
---