"""
Assignment 4: Agent-Based Model for Surface Panelization
Author: Magnus Støvring

"""

"""GhPython Script: Agent Simulator
Purpose:
    This component simulates agents moving across a Rhino surface.
    Each agent has a position and velocity, and responds to four forces:
        - Separation: avoid crowding with neighbors
        - Curvature: move toward areas of high absolute Gaussian curvature
        - Repulsion: steer away from obstacle points
        - Boundary: stay inside the surface domain
    The simulation runs step by step when triggered, and agent state
    is stored persistently using Grasshopper's sc.sticky dictionary.

Inputs:
    Surface          : Rhino surface (rg.Surface or rg.BrepFace)
    Agents           : list of Agent objects (from builder)
    RepulsorPts      : list of Point3d obstacles
    MinDistance      : float, separation threshold
    StepSize         : float, movement step size
    SeparationWeight : float, weight for separation force
    CurvatureWeight  : float, weight for curvature force
    RepulsionWeight  : float, weight for repulsion force
    BoundaryWeight   : float, weight for boundary force
    RunStep          : bool, trigger simulation step
    Reset            : bool, reset simulation

Outputs:
    UpdatedAgents : list of Agent objects with updated positions/velocities
    Positions     : list of Point3d positions for visualization
"""

import Rhino.Geometry as rg
import random
import scriptcontext as sc

# -- Agent class --
class Agent:
    """Represents a single agent with position, velocity, and frozen state."""
    def __init__(self, position, velocity):
        self.position = position
        self.velocity = velocity
        self.frozen = False

# -- Utility functions --
def _domains(surface):
    """Return U and V parameter domains of the surface."""
    du = surface.Domain(0)
    dv = surface.Domain(1)
    return du, dv

def _clamp_param(du, dv, u, v):
    """Clamp (u, v) parameters to stay inside the surface domain."""
    if u < du.T0: u = du.T0
    if u > du.T1: u = du.T1
    if v < dv.T0: v = dv.T0
    if v > dv.T1: v = dv.T1
    return u, v

def is_near_edge(surface, pt, threshold=1e-3):
    """Check if a point lies close to the parametric edge of the surface."""
    success, u, v = surface.ClosestPoint(pt)
    if not success: return False
    du, dv = _domains(surface)
    return (abs(u - du.T0) < threshold or abs(u - du.T1) < threshold or
            abs(v - dv.T0) < threshold or abs(v - dv.T1) < threshold)

def clamp_to_surface(surface, pt):
    """Project a point back onto the surface if it drifts away."""
    success, u, v = surface.ClosestPoint(pt)
    return surface.PointAt(u, v) if success else pt

def compute_separation(agent, agents, min_dist):
    """Compute a force vector that pushes agents apart if they are too close."""
    vec = rg.Vector3d(0, 0, 0)
    for other in agents:
        if other is agent or getattr(other, "frozen", False): continue
        d = agent.position.DistanceTo(other.position)
        if 0.0 < d < min_dist:
            away = agent.position - other.position
            if away.IsZero: continue
            away.Unitize()
            vec += away * (min_dist - d)
    return vec

def compute_curvature_gradient(surface, pt, delta=0.01):
    """
    Approximate gradient of ABSOLUTE Gaussian curvature at a point.
    Samples curvature at (u,v), (u+Δu,v), (u,v+Δv) and computes finite differences.
    Maps these differences into 3D using tangent directions derived from PointAt.
    """
    success, u, v = surface.ClosestPoint(pt)
    if not success: return rg.Vector3d(0, 0, 0)

    du, dv = _domains(surface)
    du_len = (du.T1 - du.T0)
    dv_len = (dv.T1 - dv.T0)
    if du_len <= 0 or dv_len <= 0:
        return rg.Vector3d(0, 0, 0)

    du_step = delta * du_len
    dv_step = delta * dv_len

    u1, v1 = _clamp_param(du, dv, u + du_step, v)
    u2, v2 = _clamp_param(du, dv, u, v + dv_step)

    try:
        c0 = surface.CurvatureAt(u, v)
        cu = surface.CurvatureAt(u1, v1)
        cv = surface.CurvatureAt(u2, v2)
    except:
        return rg.Vector3d(0, 0, 0)

    if not c0 or not cu or not cv:
        return rg.Vector3d(0, 0, 0)

    # Use absolute Gaussian curvature values
    k0 = abs(c0.Gaussian)
    k_u = abs(cu.Gaussian)
    k_v = abs(cv.Gaussian)

    dU = k_u - k0
    dV = k_v - k0

    p00 = surface.PointAt(u, v)
    p10 = surface.PointAt(u1, v1)
    p01 = surface.PointAt(u2, v2)

    du_vec = p10 - p00
    dv_vec = p01 - p00
    if du_vec.IsZero or dv_vec.IsZero:
        return rg.Vector3d(0, 0, 0)

    du_vec.Unitize()
    dv_vec.Unitize()

    grad = du_vec * dU + dv_vec * dV
    return grad

def compute_repulsion(agent, repulsorPts):
    """Compute repulsion force away from given obstacle points."""
    vec = rg.Vector3d(0, 0, 0)
    if not repulsorPts: return vec
    for rp in repulsorPts:
        d = agent.position.DistanceTo(rp)
        if d <= 0: continue
        away = agent.position - rp
        if away.IsZero: continue
        away.Unitize()
        vec += away / d
    return vec

def compute_boundary_force(surface, pt, threshold=0.05):
    """
    If agent is close to the parametric boundary, pull it back toward the
    parametric center of the surface.
    """
    success, u, v = surface.ClosestPoint(pt)
    if not success: return rg.Vector3d(0, 0, 0)
    du, dv = _domains(surface)
    du_len = (du.T1 - du.T0)
    dv_len = (dv.T1 - dv.T0)
    if du_len <= 0 or dv_len <= 0:
        return rg.Vector3d(0, 0, 0)

    dist_u = min(abs(u - du.T0), abs(du.T1 - u)) / du_len
    dist_v = min(abs(v - dv.T0), abs(dv.T1 - v)) / dv_len

    if dist_u < threshold or dist_v < threshold:
        center_u = 0.5 * (du.T0 + du.T1)
        center_v = 0.5 * (dv.T0 + dv.T1)
        center = surface.PointAt(center_u, center_v)
        vec = center - pt
        if vec.IsZero: return rg.Vector3d(0, 0, 0)
        vec.Unitize()
        return vec

    return rg.Vector3d(0, 0, 0)

# -- Sticky key --
KEY = "AgentSimulator::sticky_agents"

# -- Reset --
if Reset:
    new_agents = []
    if Agents and len(Agents) > 0:
        # Copy agents from builder; freeze if near edge
        for src in Agents:
            ag = Agent(src.position, src.velocity)
            ag.frozen = getattr(src, "frozen", is_near_edge(Surface, ag.position))
            new_agents.append(ag)
    sc.sticky[KEY] = new_agents

# -- Initialize storage --
if KEY not in sc.sticky:
    sc.sticky[KEY] = []
UpdatedAgents = sc.sticky[KEY]

# -- Step update --
if RunStep and UpdatedAgents:
    for agent in UpdatedAgents:
        if agent.frozen: continue

        # Compute all forces
        sep = compute_separation(agent, UpdatedAgents, MinDistance)
        curv = compute_curvature_gradient(Surface, agent.position)
        rep = compute_repulsion(agent, RepulsorPts)
        bound = compute_boundary_force(Surface, agent.position)

        # Weighted sum of forces
        move = (sep * SeparationWeight +
                curv * CurvatureWeight +
                rep * RepulsionWeight +
                bound * BoundaryWeight)

        # Update agent state
        if move.Length > 0:
            move.Unitize()
            agent.velocity = move * StepSize
            new_pos = agent.position + agent.velocity
            if is_near_edge(Surface, new_pos):
                # Stop at edge; optionally set agent.frozen = True if you want permanent freeze
                agent.velocity = rg.Vector3d(0, 0, 0)
            else:
                agent.position = clamp_to_surface(Surface, new_pos)

# -- Outputs --
sc.sticky[KEY] = UpdatedAgents
Positions = [a.position for a in UpdatedAgents]