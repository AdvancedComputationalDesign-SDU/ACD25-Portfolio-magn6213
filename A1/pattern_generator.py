import numpy as np
import matplotlib.pyplot as plt

## --- Random seed --- ##
seed = 5
np.random.seed(seed)

## --- Roughness --- ##
scale = 50

## --- Vectorized Perlin-like noise function --- ##
def perlin_noise_vectorized(height, width, scale=scale):
    grid_y = height // scale + 2
    grid_x = width // scale + 2
    grid = np.random.rand(grid_y, grid_x)

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

## --- Canvas size --- ##
height = 1000
width = 1000

## --- Generate Perlin noise --- ##
noise = perlin_noise_vectorized(height, width)

## --- Normalise noise to [0, 1] --- ##
noise = (noise - noise.min()) / (noise.max() - noise.min())

# ============================================================
#   VERSION 1: MANUAL RGB — RED CHANNEL ONLY (TOGGLE ON/OFF)
# ============================================================
""" rgb_image = np.zeros((height, width, 3))
rgb_image[..., 0] = noise   # Red channel = noise
# Green and blue remain zero

active_title = "Scale=50 , Red Channel" """


# ============================================================
#   VERSION 2: COLORMAP VERSION (TOGGLE ON/OFF)
# ============================================================
rgb_image = plt.cm.get_cmap("magma_r")(noise)  # RGBA output

active_title = "Scale=50 , Colormap"

# ============================================================
#   DISPLAY
# ============================================================

fig = plt.figure(figsize=(10, 6))
plt.imshow(rgb_image)
plt.title(active_title)
plt.axis("off")

plt.show()

## --- Save image --- ##
save_path = r"C:\Users\magn6\Documents\GitHub\ACD25-Portfolio-magn6213\A1\images"
filename = f"{save_path}\\perlin_output_seed_{seed}_Scale_{scale}_Colormap.png"

fig.savefig(filename, dpi=300, bbox_inches="tight", pad_inches=0)
print(f"Saved: {filename}")


