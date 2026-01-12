"""
Assignment 4: Agent-Based Model for Surface Panelization
Author: Magnus Støvring

"""

"""GhPython Script: Agent Builder
Purpose:
    This component lays out a grid of agents across a Rhino surface.
    Each agent is initialized at a surface point with zero velocity.
    The resulting list of Agent objects is meant to be fed into a simulator
    component that updates their positions step by step.

Inputs:
    Srf   : Surface (Rhino surface)
    DivU  : int, number of divisions in the U direction
    DivV  : int, number of divisions in the V direction

Outputs:
    Agents : list of Agent objects (each with position, velocity, and id)
    Pts    : list of Point3d positions corresponding to the grid locations
"""

import rhinoscriptsyntax as rs
import Rhino.Geometry as rg

# -- Agent Class --
class Agent:
    """
    Represents a single agent prepared for simulation.
    Stores its initial position on the surface, velocity (starts at zero),
    and placeholder fields that the simulator may use later (e.g., min_distance).
    """
    def __init__(self, position, velocity):
        self.position = position
        self.velocity = velocity
        # Placeholders for downstream simulation use
        self.min_distance = None
        self.curvature_value = None
        self.repulsor_forces = []
        # Optional identifier for tracking agents
        self.id = None

    def initialize_state(self):
        """Reserved for future expansion; not used in the builder."""
        pass

# -- Factory Function --
def build_agents(surface, div_u, div_v):
    """
    Create agents on a uniform UV grid across the given surface.
    - surface : Rhino surface
    - div_u   : number of segments along U (produces div_u+1 points)
    - div_v   : number of segments along V (produces div_v+1 points)

    Returns:
        agents : list of Agent instances with zero velocity
        points : list of Point3d locations where agents are placed
    """
    # Validate inputs
    if surface is None:
        return [], []
    if div_u < 1 or div_v < 1:
        return [], []

    agents, points = [], []
    agent_id = 0

    # Get surface parameter domains
    dom_u = surface.Domain(0)
    dom_v = surface.Domain(1)

    # Loop through a (div_u+1) x (div_v+1) grid in parameter space
    for i in range(div_u + 1):
        # Linearly interpolate U parameter across the domain
        u = dom_u.T0 + (i / float(div_u)) * (dom_u.T1 - dom_u.T0)
        for j in range(div_v + 1):
            # Linearly interpolate V parameter across the domain
            v = dom_v.T0 + (j / float(div_v)) * (dom_v.T1 - dom_v.T0)

            # Evaluate the surface at (u, v) to get a 3D point
            pt = surface.PointAt(u, v)
            points.append(pt)

            # Initialize agent with zero velocity at this point
            vel = rg.Vector3d(0, 0, 0)
            agent = Agent(pt, vel)
            agent.id = agent_id
            agent_id += 1

            agents.append(agent)

    return agents, points

# -- Outputs --
Agents, Pts = build_agents(Srf, DivU, DivV)
