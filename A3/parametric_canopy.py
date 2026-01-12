"""
Assignment 3: Parametric Structural Canopy — Pseudocode Scaffold

Author: Magnus Støvring

"""

# -------------------------------
# Variation 1) Perlin Noise Surface Generation
# -------------------------------
import Rhino.Geometry as rg
import numpy as np


## PERLIN NOISE ##
def perlin_noise(width, height, scale=10, seed=0):
    rng = np.random.default_rng(seed)

    # Random gradient vectors at each grid point #
    gx = rng.random((width+1, height+1)) * 2 - 1
    gy = rng.random((width+1, height+1)) * 2 - 1
    grad = np.stack((gx, gy), axis=-1)

    # Pixel coordinates (scaled domain) #
    xs = np.linspace(0, scale, width, endpoint=False)
    ys = np.linspace(0, scale, height, endpoint=False)
    xx, yy = np.meshgrid(xs, ys, indexing="ij")

    # Corner indices #
    x0 = xx.astype(int)
    y0 = yy.astype(int)
    x1 = x0 + 1
    y1 = y0 + 1

    # Local coordinates inside each cell #
    dx = xx - x0
    dy = yy - y0

    # Fade function #
    def fade(t):
        return 6*t**5 - 15*t**4 + 10*t**3

    tx = fade(dx)
    ty = fade(dy)

    # Dot product helper #
    def dot(ix, iy, x, y):
        g = grad[ix % grad.shape[0], iy % grad.shape[1]]
        return g[...,0]*x + g[...,1]*y

    # Dot at each of the 4 corners #
    n00 = dot(x0, y0, dx, dy)
    n10 = dot(x1, y0, dx-1, dy)
    n01 = dot(x0, y1, dx, dy-1)
    n11 = dot(x1, y1, dx-1, dy-1)

    # Bilinear interpolation using fade weights #
    nx0 = n00*(1-tx) + n10*tx
    nx1 = n01*(1-tx) + n11*tx
    nxy = nx0*(1-ty) + nx1*ty

    # Normalize to [0,1] #
    nxy = (nxy - nxy.min()) / (nxy.max() - nxy.min())

    return nxy



## BUILD THE UV GRID ##
u = np.linspace(0, 1, ResU)
v = np.linspace(0, 1, ResV)
uu, vv = np.meshgrid(u, v, indexing="ij")


## GENERATE HEIGHTMAP FROM PERLIN NOISE ##
H = perlin_noise(ResU, ResV, scale=Scale, seed=Seed)
H = Amp * H + Height

## BUILD POINT GRID USING PERLIN HEIGHTMAP ##
pt_grid = []
for i in range(ResU):
    row = []
    for j in range(ResV):
        x = float(i)
        y = float(j)
        z = float(H[i, j])
        row.append(rg.Point3d(x, y, z))
    pt_grid.append(row)

# Flatten for NURBS constructor #
flat_pts = [p for row in pt_grid for p in row]


## MAKE NURBS SURFACE FROM POINT GRID ##
Surface = rg.NurbsSurface.CreateFromPoints(flat_pts, ResU, ResV, 3, 3)









# -------------------------------
# Variation 2) Radial Sinus Heightmap
# -------------------------------
import Rhino.Geometry as rg
import numpy as np

## RADIAL HEIGHTMAP ##
def radial_heightmap(ResU, ResV, amplitude=1.0, frequency=1.0, height=0.0):
    xs = np.linspace(-1, 1, ResU)
    ys = np.linspace(-1, 1, ResV)
    xx, yy = np.meshgrid(xs, ys, indexing="ij")

    R = np.sqrt(xx**2 + yy**2)
    Z = amplitude * np.sin(frequency * R * np.pi * 2)

    # Normalize to [0,1]
    Z = (Z - Z.min()) / (Z.max() - Z.min())

    # Apply height offset
    Z = Z + height

    return Z

## HEIGHTMAP ##
H = radial_heightmap(ResU, ResV, Amp, Frequency, Height)

## POINT GRID ##
pt_grid = []
for i in range(ResU):
    row = []
    for j in range(ResV):
        row.append(rg.Point3d(float(i), float(j), float(H[i, j])))
    pt_grid.append(row)

flat_pts = [p for row in pt_grid for p in row]

## SURFACE ##
Surface = rg.NurbsSurface.CreateFromPoints(flat_pts, ResU, ResV, 3, 3)





# -------------------------------
# Variation 3) Gaussian Bump Heightmap
# -------------------------------
import Rhino.Geometry as rg
import numpy as np

## GAUSSIAN BUMP FIELD ##
def gaussian_bumps(ResU, ResV, amplitude=1.0, num_centers=3, spread=0.2, seed=0, height=0.0):
    rng = np.random.default_rng(seed)

    # Random bump centers in normalized domain
    centers = rng.random((num_centers, 2))

    xs = np.linspace(0, 1, ResU)
    ys = np.linspace(0, 1, ResV)
    xx, yy = np.meshgrid(xs, ys, indexing="ij")

    Z = np.zeros((ResU, ResV))

    for cx, cy in centers:
        dist2 = (xx - cx)**2 + (yy - cy)**2
        Z += np.exp(-dist2 / spread)

    # Normalize
    Z = (Z - Z.min()) / (Z.max() - Z.min())

    # Apply amplitude
    Z = Z * amplitude

    # Apply height offset
    Z = Z + height

    return Z

## HEIGHTMAP ##
H = gaussian_bumps(ResU, ResV, Amp, NumCenters, Spread, Seed, Height)

## POINT GRID ##
pt_grid = []
for i in range(ResU):
    row = []
    for j in range(ResV):
        row.append(rg.Point3d(float(i), float(j), float(H[i, j])))
    pt_grid.append(row)

flat_pts = [p for row in pt_grid for p in row]

## SURFACE ##
Surface = rg.NurbsSurface.CreateFromPoints(flat_pts, ResU, ResV, 3, 3)








# -------------------------------
# Variation 1) Voronoi Surface Tesselation
# -------------------------------
import rhinoscriptsyntax as rs
import Rhino
import random
import math
import heapq
import scriptcontext as sc

EPS = 1e-9
random.seed(Seed)


# CURVATURE-WEIGHTED SITE GENERATION WITH MIN DISTANCE

def generate_curvature_weighted_sites(surface, count, grid=100, min_dist=0.0):
    du = surface.Domain(0)
    dv = surface.Domain(1)
    umin, umax = du.T0, du.T1
    vmin, vmax = dv.T0, dv.T1
    uv_list = []
    weights = []
    for i in range(grid):
        u = umin + (i / float(grid - 1)) * (umax - umin)
        for j in range(grid):
            v = vmin + (j / float(grid - 1)) * (vmax - vmin)
            curv_obj = surface.CurvatureAt(u, v)
            mp = abs(curv_obj.Mean) if curv_obj else 0.0
            mp = max(mp, 0.001)
            uv_list.append((u, v))
            weights.append(mp)
    total = sum(weights)
    weights = [w / total for w in weights]
    sites = []
    attempts = 0
    max_attempts = count * 200
    while len(sites) < count and attempts < max_attempts:
        candidate = random.choices(uv_list, weights=weights, k=1)[0]
        too_close = False
        for s in sites:
            p1 = surface.PointAt(candidate[0], candidate[1])
            p2 = surface.PointAt(s[0], s[1])
            if p1.DistanceTo(p2) < min_dist:
                too_close = True
                break
        if not too_close:
            sites.append(candidate)
        attempts += 1
    return sites


# NEAREST-NEIGHBOR DISTANCE SUM

def sum_nearest_distances(points, k=10):
    n = len(points)
    sums = []
    for i, pt in enumerate(points):
        distances = [pt.DistanceTo(points[j]) for j in range(n) if j != i]
        nearest_k = heapq.nsmallest(k, distances)
        sums.append(sum(nearest_k))
    return sums

def select_low_distance_points(points, total_distances, count, min_dist):
    pts_with_dist = list(zip(points, total_distances))
    pts_with_dist.sort(key=lambda x: x[1])
    selected = []
    for pt, _ in pts_with_dist:
        if len(selected) >= count:
            break
        too_close = False
        for s in selected:
            if pt.DistanceTo(s) < min_dist:
                too_close = True
                break
        if not too_close:
            selected.append(pt)
    return selected


# VORONOI BY BISECTOR LINES

def perp_bisector_line(a, b):
    sx, sy = a
    ox, oy = b
    mx, my = (sx + ox) / 2.0, (sy + oy) / 2.0
    dx, dy = ox - sx, oy - sy
    n = math.hypot(dx, dy)
    if n < EPS: return None
    A, B = dx / n, dy / n
    C = -(A * mx + B * my)
    return A, B, C

def same_side(pt, A, B, C, keep_pt):
    x, y = pt
    kx, ky = keep_pt
    s1 = A * x + B * y + C
    s2 = A * kx + B * ky + C
    if abs(s1) < EPS:
        return True
    return (s1 >= 0) == (s2 >= 0)

def line_seg_intersection(A, B, C, p1, p2):
    x1, y1 = p1
    x2, y2 = p2
    dx, dy = x2 - x1, y2 - y1
    denom = A * dx + B * dy
    if abs(denom) < EPS: return None
    t = -(A * x1 + B * y1 + C) / denom
    if 0 <= t <= 1:
        return (x1 + t*dx, y1 + t*dy)
    return None

def clip_poly(poly, A, B, C, keep_pt):
    if not poly:
        return []
    out = []
    m = len(poly)
    for i in range(m):
        P = poly[i]
        Q = poly[(i + 1) % m]
        Pin = same_side(P, A, B, C, keep_pt)
        Qin = same_side(Q, A, B, C, keep_pt)
        if Pin and Qin:
            out.append(Q)
        elif Pin and not Qin:
            ip = line_seg_intersection(A, B, C, P, Q)
            if ip: out.append(ip)
        elif not Pin and Qin:
            ip = line_seg_intersection(A, B, C, P, Q)
            if ip: out.append(ip)
            out.append(Q)
    return out

def voronoi_on_uv(sites, domain_u, domain_v):
    umin, umax = domain_u
    vmin, vmax = domain_v
    box = [(umin, vmax), (umax, vmax), (umax, vmin), (umin, vmin)]
    cells = []
    for i, site in enumerate(sites):
        cell = box[:]
        for j, other in enumerate(sites):
            if i == j:
                continue
            bis = perp_bisector_line(site, other)
            if bis is None:
                continue
            A, B, C = bis
            cell = clip_poly(cell, A, B, C, site)
            if not cell:
                break
        cells.append(cell)
    return cells


# EXECUTION

try:
    du = Srf.Domain(0)
    dv = Srf.Domain(1)
    domain_u = (du.T0, du.T1)
    domain_v = (dv.T0, dv.T1)

    # Generate Voronoi sites
    sites_uv = generate_curvature_weighted_sites(Srf, N, min_dist=MinDist)
    Sites3d = [Srf.PointAt(u, v) for u, v in sites_uv]

    # Voronoi polygons
    vor_uv = voronoi_on_uv(sites_uv, domain_u, domain_v)
    VoronoiCurves3d = []
    for poly in vor_uv:
        pts3d = []
        for u, v in poly:
            pts3d.append(Srf.PointAt(u, v))
        if len(pts3d) >= 3:
            pts3d.append(pts3d[0])
            VoronoiCurves3d.append(Rhino.Geometry.Polyline(pts3d).ToNurbsCurve())

    # Create pipes
    VoronoiPipeMeshes = []
    tolerance = sc.doc.ModelAbsoluteTolerance
    angle_tol = sc.doc.ModelAngleToleranceRadians
    for crv in VoronoiCurves3d:
        if crv and crv.IsValid and crv.GetLength() > 1e-6:
            breps = Rhino.Geometry.Brep.CreatePipe(
                crv, PipeRadius, False,
                Rhino.Geometry.PipeCapMode.Round,
                True, tolerance, angle_tol
            )
            if breps:
                VoronoiPipeMeshes.extend(breps)

    # Select points based on nearest neighbors
    total_dists = sum_nearest_distances(Sites3d, k=10)
    SelectedPoints3d = select_low_distance_points(Sites3d, total_dists, SelectedCount, MinSpacing)

except Exception as e:
    print("Error:", e)
    Sites3d = []
    SelectedPoints3d = []
    VoronoiCurves3d = []
    VoronoiPipeMeshes = []


# UNIQUE POINTS & RANDOM CLOSEST SELECTION

try:
    tol = 0.01  # duplicate tolerance
    buckets = {}
    UniquePoints = []

    def bucket_key(pt, tol):
        return (int(round(pt.X / tol)),
                int(round(pt.Y / tol)),
                int(round(pt.Z / tol)))

    # Use SelectedPoints3d as GivenPoints
    GivenPoints = SelectedPoints3d
    Curves = VoronoiCurves3d

    for crv in Curves:
        if crv is None:
            continue
        rc, pl = crv.TryGetPolyline()
        if not rc:
            t = Rhino.Geometry.Polyline()
            if crv.ToPolyline(t):
                pl = t
            else:
                continue
        for pt in pl:
            key = bucket_key(pt, tol)
            if key not in buckets:
                buckets[key] = pt
                UniquePoints.append(pt)

    # Select 6 random closest points for each GivenPoint
    SelectedPoints = []
    for gp in GivenPoints:
        dists = [(gp.DistanceTo(up), up) for up in UniquePoints]
        closest30 = heapq.nsmallest(30, dists, key=lambda x: x[0])
        closest_pts = [pt for dist, pt in closest30]
        selected = random.sample(closest_pts, 6)
        SelectedPoints.extend(selected)

except Exception as e:
    print("Error in closest points selection:", e)
    UniquePoints = []
    SelectedPoints = []
    
# Output the original surface
OriginalSurface = Srf






# -------------------------------
# Variation 2) Rectangular Surface Tesselation
# -------------------------------
import Rhino
import scriptcontext as sc
import random
import heapq

random.seed(Seed)


# NEAREST-NEIGHBOR DISTANCE SUM 

def sum_nearest_distances(points, k=10):
    n = len(points)
    sums = []
    for i, pt in enumerate(points):
        distances = [pt.DistanceTo(points[j]) for j in range(n) if j != i]
        nearest_k = heapq.nsmallest(k, distances)
        sums.append(sum(nearest_k))
    return sums

def select_low_distance_points(points, total_distances, count, min_dist):
    pts_with_dist = list(zip(points, total_distances))
    pts_with_dist.sort(key=lambda x: x[1])
    selected = []
    for pt, _ in pts_with_dist:
        if len(selected) >= count:
            break
        if all(pt.DistanceTo(s) >= min_dist for s in selected):
            selected.append(pt)
    return selected


# RECTANGULAR TESSELLATION

try:
    du = Srf.Domain(0)
    dv = Srf.Domain(1)
    umin, umax = du.T0, du.T1
    vmin, vmax = dv.T0, dv.T1

    # Limit grid resolution to avoid crashes
    Nu = max(2, min(int(N), 40))
    Nv = max(2, min(int(N), 40))

    u_vals = [umin + (i / float(Nu - 1)) * (umax - umin) for i in range(Nu)]
    v_vals = [vmin + (j / float(Nv - 1)) * (vmax - vmin) for j in range(Nv)]

    # Build grid of points
    grid = [[Srf.PointAt(u_vals[i], v_vals[j]) for j in range(Nv)] for i in range(Nu)]

    # Sites = all grid vertices
    Sites3d = [pt for row in grid for pt in row]

    # Build only CELL BOUNDARIES (4 edges per cell)
    RectCurves = []
    for i in range(Nu - 1):
        for j in range(Nv - 1):
            p00 = grid[i][j]
            p10 = grid[i+1][j]
            p11 = grid[i+1][j+1]
            p01 = grid[i][j+1]

            # 4 edges
            edges = [
                [p00, p10],
                [p10, p11],
                [p11, p01],
                [p01, p00]
            ]

            for e in edges:
                pl = Rhino.Geometry.Polyline(e)
                RectCurves.append(pl.ToNurbsCurve())

    VoronoiCurves3d = RectCurves

    # Pipe curves (much fewer now)
    VoronoiPipeMeshes = []
    tol = sc.doc.ModelAbsoluteTolerance
    ang = sc.doc.ModelAngleToleranceRadians

    for crv in VoronoiCurves3d:
        if crv and crv.IsValid and crv.GetLength() > 1e-6:
            breps = Rhino.Geometry.Brep.CreatePipe(
                crv, PipeRadius, False,
                Rhino.Geometry.PipeCapMode.Round,
                True, tol, ang
            )
            if breps:
                VoronoiPipeMeshes.extend(breps)

    # Select points based on nearest neighbors
    total_dists = sum_nearest_distances(Sites3d, k=10)
    SelectedPoints3d = select_low_distance_points(Sites3d, total_dists, SelectedCount, MinSpacing)

except Exception as e:
    print("Error:", e)
    Sites3d = []
    SelectedPoints3d = []
    VoronoiCurves3d = []
    VoronoiPipeMeshes = []


# UNIQUE POINTS & RANDOM CLOSEST SELECTION 

try:
    tol = 0.01
    buckets = {}
    UniquePoints = []

    def bucket_key(pt, tol):
        return (int(round(pt.X / tol)),
                int(round(pt.Y / tol)),
                int(round(pt.Z / tol)))

    GivenPoints = SelectedPoints3d
    Curves = VoronoiCurves3d

    for crv in Curves:
        rc, pl = crv.TryGetPolyline()
        if not rc:
            t = Rhino.Geometry.Polyline()
            if crv.ToPolyline(t):
                pl = t
            else:
                continue
        for pt in pl:
            key = bucket_key(pt, tol)
            if key not in buckets:
                buckets[key] = pt
                UniquePoints.append(pt)

    SelectedPoints = []
    for gp in GivenPoints:
        dists = [(gp.DistanceTo(up), up) for up in UniquePoints]
        closest30 = heapq.nsmallest(30, dists, key=lambda x: x[0])
        closest_pts = [pt for dist, pt in closest30]
        if len(closest_pts) >= 6:
            SelectedPoints.extend(random.sample(closest_pts, 6))

except Exception as e:
    print("Error in closest points selection:", e)
    UniquePoints = []
    SelectedPoints = []

OriginalSurface = Srf





# -------------------------------
# Variation 3) Diagrid Surface Tesselation
# -------------------------------
import rhinoscriptsyntax as rs
import Rhino
import random
import math
import heapq
import scriptcontext as sc

EPS = 1e-9
random.seed(Seed)


# NEAREST-NEIGHBOR DISTANCE SUM (REUSED)

def sum_nearest_distances(points, k=10):
    n = len(points)
    sums = []
    for i, pt in enumerate(points):
        distances = [pt.DistanceTo(points[j]) for j in range(n) if j != i]
        nearest_k = heapq.nsmallest(k, distances)
        sums.append(sum(nearest_k))
    return sums

def select_low_distance_points(points, total_distances, count, min_dist):
    pts_with_dist = list(zip(points, total_distances))
    pts_with_dist.sort(key=lambda x: x[1])
    selected = []
    for pt, _ in pts_with_dist:
        if len(selected) >= count:
            break
        too_close = False
        for s in selected:
            if pt.DistanceTo(s) < min_dist:
                too_close = True
                break
        if not too_close:
            selected.append(pt)
    return selected


# DIAGRID TESSELLATION

try:
    du = Srf.Domain(0)
    dv = Srf.Domain(1)
    umin, umax = du.T0, du.T1
    vmin, vmax = dv.T0, dv.T1

    # Use N as number of divisions in each direction
    Nu = max(2, int(N))
    Nv = max(2, int(N))

    u_vals = [umin + (i / float(Nu - 1)) * (umax - umin) for i in range(Nu)]
    v_vals = [vmin + (j / float(Nv - 1)) * (vmax - vmin) for j in range(Nv)]

    # Build grid of points on the surface
    grid_points = []
    for i in range(Nu):
        row = []
        for j in range(Nv):
            pt = Srf.PointAt(u_vals[i], v_vals[j])
            row.append(pt)
        grid_points.append(row)

    # Sites: all grid vertices (you could also use cell centers if you prefer)
    Sites3d = [pt for row in grid_points for pt in row]

    TriCurves3d = []

    # Base grid edges (optional but helpful structurally)
    for i in range(Nu):
        pts = [grid_points[i][j] for j in range(Nv)]
        if len(pts) >= 2:
            pl = Rhino.Geometry.Polyline(pts)
            TriCurves3d.append(pl.ToNurbsCurve())
    for j in range(Nv):
        pts = [grid_points[i][j] for i in range(Nu)]
        if len(pts) >= 2:
            pl = Rhino.Geometry.Polyline(pts)
            TriCurves3d.append(pl.ToNurbsCurve())

    # Diagonals across each quad (forming triangles / diagrid)
    for i in range(Nu - 1):
        for j in range(Nv - 1):
            p00 = grid_points[i][j]
            p10 = grid_points[i + 1][j]
            p11 = grid_points[i + 1][j + 1]
            p01 = grid_points[i][j + 1]

            # Diagonal 1: p00 -> p11
            pl1 = Rhino.Geometry.Polyline([p00, p11])
            TriCurves3d.append(pl1.ToNurbsCurve())

            # Diagonal 2: p10 -> p01
            pl2 = Rhino.Geometry.Polyline([p10, p01])
            TriCurves3d.append(pl2.ToNurbsCurve())

    # Create pipes along all curves
    VoronoiCurves3d = TriCurves3d  # reuse output name
    VoronoiPipeMeshes = []
    tolerance = sc.doc.ModelAbsoluteTolerance
    angle_tol = sc.doc.ModelAngleToleranceRadians

    for crv in VoronoiCurves3d:
        if crv and crv.IsValid and crv.GetLength() > 1e-6:
            breps = Rhino.Geometry.Brep.CreatePipe(
                crv, PipeRadius, False,
                Rhino.Geometry.PipeCapMode.Round,
                True, tolerance, angle_tol
            )
            if breps:
                VoronoiPipeMeshes.extend(breps)

    # Select points based on nearest neighbors
    total_dists = sum_nearest_distances(Sites3d, k=10)
    SelectedPoints3d = select_low_distance_points(Sites3d, total_dists, SelectedCount, MinSpacing)

except Exception as e:
    print("Error:", e)
    Sites3d = []
    SelectedPoints3d = []
    VoronoiCurves3d = []
    VoronoiPipeMeshes = []


# UNIQUE POINTS & RANDOM CLOSEST SELECTION 

try:
    tol = 0.01  # duplicate tolerance
    buckets = {}
    UniquePoints = []

    def bucket_key(pt, tol):
        return (int(round(pt.X / tol)),
                int(round(pt.Y / tol)),
                int(round(pt.Z / tol)))

    GivenPoints = SelectedPoints3d
    Curves = VoronoiCurves3d

    for crv in Curves:
        if crv is None:
            continue
        rc, pl = crv.TryGetPolyline()
        if not rc:
            t = Rhino.Geometry.Polyline()
            if crv.ToPolyline(t):
                pl = t
            else:
                continue
        for pt in pl:
            key = bucket_key(pt, tol)
            if key not in buckets:
                buckets[key] = pt
                UniquePoints.append(pt)

    SelectedPoints = []
    for gp in GivenPoints:
        dists = [(gp.DistanceTo(up), up) for up in UniquePoints]
        closest30 = heapq.nsmallest(30, dists, key=lambda x: x[0])
        closest_pts = [pt for dist, pt in closest30]
        if len(closest_pts) >= 6:
            selected = random.sample(closest_pts, 6)
            SelectedPoints.extend(selected)

except Exception as e:
    print("Error in closest points selection:", e)
    UniquePoints = []
    SelectedPoints = []

# Output the original surface
OriginalSurface = Srf





# -------------------------------
# Variation 1, 2 & 3) Branching Support Structure
# -------------------------------
import Rhino.Geometry as rg
import random

## PARAMETERS ##
random.seed(Seed)

PipeBreps = []         # Output: GH-previewable Breps
up_first = 1.5         # vertical trunk height
noise = 0.5            # maximum noise magnitude
use_noise = True

## SURFACE ##
if not isinstance(srf, rg.Surface):
    raise Exception("Input srf must be a Rhino Surface!")

## DEFINE ROOT POINTS FROM SURFACE ##
u_domain = srf.Domain(0)
v_domain = srf.Domain(1)

uv_points_normalized = [(0.25, 0.25), (0.75, 0.75)]
roots = []
for u_norm, v_norm in uv_points_normalized:
    u = u_domain.T0 + u_norm * (u_domain.T1 - u_domain.T0)
    v = v_domain.T0 + v_norm * (v_domain.T1 - v_domain.T0)
    pt = srf.PointAt(u, v)
    projected_pt = rg.Point3d(pt.X, pt.Y, 0)
    roots.append(projected_pt)

## SURFACE TO MESH (for trimming) ##
Canopy = rg.Mesh.CreateFromSurface(srf, rg.MeshingParameters.Default)

## TRIM FUNCTION ##
def trim_line_to_canopy(start_pt, end_pt, canopy):
    line = rg.Line(start_pt, end_pt)
    hits = rg.Intersect.Intersection.MeshLine(canopy, line)
    if hits:
        pt = hits[0] if isinstance(hits[0], rg.Point3d) else hits[0].Point
        return pt, True
    return end_pt, False

## ATTRACTOR ASSIGNMENT ##
def assign_attractors(attractors, roots):
    mapping = {pt: [] for pt in roots}
    for a in attractors:
        closest = min(roots, key=lambda r: r.DistanceTo(a))
        mapping[closest].append(a)
    return mapping

def split_into_two(attractors):
    if len(attractors) <= 1:
        return attractors, []

    max_dist = -1
    seedA, seedB = attractors[0], attractors[1]
    for i in range(len(attractors)):
        for j in range(i+1, len(attractors)):
            d = attractors[i].DistanceTo(attractors[j])
            if d > max_dist:
                max_dist = d
                seedA, seedB = attractors[i], attractors[j]

    groupA = [seedA]
    groupB = [seedB]

    for a in attractors:
        if a in (seedA, seedB):
            continue
        if a.DistanceTo(seedA) < a.DistanceTo(seedB):
            groupA.append(a)
        else:
            groupB.append(a)

    return groupA, groupB

## CLUSTER DIRECTION WITH HEIGHT-DEPENDENT NOISE ##
def cluster_direction(from_pt, attractors, level, max_level):
    v = rg.Vector3d(0,0,0)
    for a in attractors:
        v += (a - from_pt)
    if v.IsZero:
        v = rg.Vector3d(0,0,1)
    v.Unitize()

    if use_noise:
        noise_scale = noise * (level / max_level)
        v.X += (random.random() - 0.5) * noise_scale
        v.Y += (random.random() - 0.5) * noise_scale
        v.Z += (random.random() - 0.5) * noise_scale
        v.Unitize()

    return v

## PIPE CREATION ##
def branch_radius(level, max_level):
    tip_radius = 0.05
    # Linear interpolation from Thickness to tip_radius
    r = tip_radius + (Thickness - tip_radius) * (level / max_level)
    return r


def line_to_pipe(start_pt, end_pt, level, max_level):
    r = branch_radius(level, max_level)
    line_vec = end_pt - start_pt
    length = line_vec.Length
    if length < 1e-6:
        return None
    direction = rg.Vector3d(line_vec)
    direction.Unitize()
    plane = rg.Plane(start_pt, direction)
    circle = rg.Circle(plane, r)
    cyl = rg.Cylinder(circle, length)
    brep = cyl.ToBrep(True, True)
    return brep

## RECURSIVE BINARY BRANCHING ##
def grow_binary(start_pt, level, attractors, max_level):
    if level == 0 or len(attractors) == 0:
        return

    # SNAP last-level branches directly to attractor
    if level == 1 and attractors:
        intended = attractors[0]
        end_pt, hit = trim_line_to_canopy(start_pt, intended, Canopy)
        pipe = line_to_pipe(start_pt, end_pt, level, max_level)
        if pipe:
            PipeBreps.append(pipe)
        return

    if len(attractors) == 1:
        v = cluster_direction(start_pt, attractors, level, max_level)
        intended = start_pt + v * Step
        end_pt, hit = trim_line_to_canopy(start_pt, intended, Canopy)
        pipe = line_to_pipe(start_pt, end_pt, level, max_level)
        if pipe:
            PipeBreps.append(pipe)
        if not hit:
            grow_binary(end_pt, level-1, attractors, max_level)
        return

    groupA, groupB = split_into_two(attractors)

    vA = cluster_direction(start_pt, groupA, level, max_level)
    intendedA = start_pt + vA * Step
    endA, hitA = trim_line_to_canopy(start_pt, intendedA, Canopy)
    pipeA = line_to_pipe(start_pt, endA, level, max_level)
    if pipeA:
        PipeBreps.append(pipeA)
    if not hitA:
        grow_binary(endA, level-1, groupA, max_level)

    vB = cluster_direction(start_pt, groupB, level, max_level)
    intendedB = start_pt + vB * Step
    endB, hitB = trim_line_to_canopy(start_pt, intendedB, Canopy)
    pipeB = line_to_pipe(start_pt, endB, level, max_level)
    if pipeB:
        PipeBreps.append(pipeB)
    if not hitB:
        grow_binary(endB, level-1, groupB, max_level)

## TREE GROWTH ##
def grow_tree(root_pt, attractors, levels):
    trunk_top = rg.Point3d(root_pt.X, root_pt.Y, root_pt.Z + up_first)
    trunk_end, hit = trim_line_to_canopy(root_pt, trunk_top, Canopy)
    trunk_pipe = line_to_pipe(root_pt, trunk_end, levels, levels)
    if trunk_pipe:
        PipeBreps.append(trunk_pipe)
    if hit:
        return
    grow_binary(trunk_end, levels-1, attractors, levels)

## MAIN EXECUTION ##
attr_map = assign_attractors(Attractors, roots)

for r in roots:
    grow_tree(r, attr_map[r], Levels)

