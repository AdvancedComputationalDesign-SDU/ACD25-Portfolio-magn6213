"""
Assignment 2: Fractal Generator

Author: Magnus Støvring

Description:
This script generates fractal patterns using recursive functions and geometric transformations.
"""


import random
import numpy as np
import matplotlib.pyplot as plt
import time 

# Start timer #
start_time = time.time()

# Randomness #
SEED = 42
random.seed(SEED)
np.random.seed(SEED)

# Parameters #
circle_radius = 50
num_particles = 500
step_size = 1
max_attempts = 900

# Attractor point #
attractor_point = np.array([circle_radius, circle_radius])
min_attach = 1
max_attach = 1.5

# Starting point #
aggregated_particles = [(0, 0)]

# Settings for matplotlib #
plt.ion()
fig, ax = plt.subplots()
ax.set_xlim(-circle_radius, circle_radius)
ax.set_ylim(-circle_radius, circle_radius)
ax.set_aspect('equal')
ax.set_facecolor('black')
fig.patch.set_facecolor('black')

# Colormap #
cmap = plt.get_cmap('plasma')

# Helper functions #
def is_inside_circle(x, y, radius):
    return x**2 + y**2 < radius**2

def random_walk(x, y, step_size):
    angle = np.random.uniform(0, 2*np.pi)
    return x + np.cos(angle)*step_size, y + np.sin(angle)*step_size

def is_near_aggregated(x, y, aggregated_particles, step_size):
    positions = np.array(aggregated_particles)

    dist_to_attractor = np.linalg.norm(attractor_point - np.array([x, y]))
    max_dist = circle_radius * np.sqrt(2)
    closeness = 1 - np.clip(dist_to_attractor / max_dist, 0, 1)
    attach_radius = min_attach + closeness * (max_attach - min_attach)
    threshold = attach_radius**2 * 1.3

    dx = positions[:, 0] - x
    dy = positions[:, 1] - y
    dist_sq = dx**2 + dy**2

    return np.any(dist_sq <= threshold)

def compute_spawn_radius(aggregated_particles, margin=10):
    max_dist = max(px**2 + py**2 for px, py in aggregated_particles)
    return np.sqrt(max_dist) + margin

# ==============================
# Recursive diffusion function
# ==============================
def walk_particle(x, y, attempts):
    # Base cases
    if attempts >= max_attempts:
        return None

    if not is_inside_circle(x, y, circle_radius):
        return None

    if is_near_aggregated(x, y, aggregated_particles, step_size):
        return (x, y)

    # Recursive step
    x_new, y_new = random_walk(x, y, step_size)
    return walk_particle(x_new, y_new, attempts + 1)

# ==============================
# Main aggregation loop
# ==============================
positions = np.array(aggregated_particles)

for i in range(num_particles):
    spawn_radius = compute_spawn_radius(aggregated_particles)
    angle = np.random.uniform(0, 2*np.pi)
    x, y = spawn_radius*np.cos(angle), spawn_radius*np.sin(angle)

    result = walk_particle(x, y, 0)

    if result is not None:
        aggregated_particles.append(result)
        positions = np.vstack([positions, result])

# Plot aggregated particles as lines #
max_width = 1.5
min_width = 0.25

all_positions = np.array(aggregated_particles)

for i in range(1, len(aggregated_particles)):
    x1, y1 = all_positions[i]

    prev_positions = all_positions[:i]
    dx = prev_positions[:, 0] - x1
    dy = prev_positions[:, 1] - y1
    dist_sq = dx**2 + dy**2

    nearest_index = np.argmin(dist_sq)
    x0, y0 = prev_positions[nearest_index]

    line_width = max_width - (max_width - min_width) * (i / len(aggregated_particles))
    ax.plot(
        [x0, x1], [y0, y1],
        color=cmap(np.hypot(x1, y1) / circle_radius),
        linewidth=line_width
    )

# End timer #
end_time = time.time()
print(f"Runtime: {end_time - start_time:.2f} seconds")

# Legend #
legend_text = (
    f"Seed: {SEED}\n"
    f"Particles: {num_particles}\n"
    f"Attempts: {max_attempts}\n"
    f"Min attach: {min_attach}\n"
    f"Max attach: {max_attach}\n"
    f"Runtime: {end_time - start_time:.2f} s"
)

fig.text(
    0.5, 0.15, legend_text,
    fontsize=5, color='white',
    va='bottom', ha='center',
    bbox=dict(facecolor='black', alpha=0.7, boxstyle='round,pad=0.3')
)

# Save figure #
plt.savefig(
    f"A2/images/DLA_seed{SEED}_p{num_particles}_att{max_attempts}_a{min_attach}-{max_attach}.png",
    dpi=300,
    bbox_inches='tight'
)
plt.ioff()
plt.show()
