
from manim import *
import numpy as np


# ============================================================================
# Shared
# ============================================================================


def parabola(x, x1, x2, ymin):
    # Vertex is at the midpoint of the intercepts
    xv = (x1 + x2) / 2
    # Calculate 'a' using the fact that V(xv) = ymin
    a = ymin / ((xv - x1) * (xv - x2))
    return a * (x - x1) * (x - x2)

def potential(r, ymin=1.0, sigma=1.0):
    """
    Calculate potential energy using a parabolic profile.
    
    The parabola is configured so that:
        • The potential crosses zero at r = sigma and r = 3 * sigma
        • The minimum value is -epsilon at r = 2 * sigma
    """

    x1 = sigma
    x2 = 3 * sigma
    ymin = -epsilon

    return parabola(r, x1, x2, ymin)

def force(r_vec):
    """
    Calculate the force between two particles from the parabolic potential.

    Returns:
        Force vector on particle i due to particle j
    """
    r = np.linalg.norm(r_vec)
    if r < 1e-6:
        return np.zeros_like(r_vec)

    x1 = sigma
    x2 = 3 * sigma
    ymin = -epsilon
    xv = (x1 + x2) / 2
    a = ymin / ((xv - x1) * (xv - x2))

    dV_dr = a * (2 * r - x1 - x2)
    force_magnitude = -dV_dr

    return force_magnitude * (r_vec / r)


# ============================================================================
# Scenes
# ============================================================================

class MD_T(Scene):
    """Enhanced version with a parabolic interaction potential graph"""
    def construct(self):
        # Simulation parameters
        n_particles = 2
        dt = 0.01
        
        # Potential parameters
        x1 = 0
        x2 = 4
        ymin = -4
        
        # Box boundaries
        box_size = 4.0  # Reduced to make room for graph
        edge_buffer = 0.4  # Keep particles away from edges
        box_center_x = -2.0  # The boundary box is shifted left by 2
        
        # Initialize particles with proper spacing
        np.random.seed(42)
        positions = np.zeros((n_particles, 2))
        min_distance = 2.5  # Allow particles to spawn farther apart for wider range of motion
        
        for i in range(n_particles):
            max_attempts = 100
            for attempt in range(max_attempts):
                # Generate position with buffer from edges, centered at box_center_x
                pos = np.array([
                    np.random.uniform(
                        box_center_x - (box_size/2 - edge_buffer),
                        box_center_x + (box_size/2 - edge_buffer)
                    ),
                    np.random.uniform(
                        -(box_size/2 - edge_buffer),
                        (box_size/2 - edge_buffer)
                    )
                ])
                
                # Check distance from all existing particles
                if i == 0:
                    positions[i] = pos
                    break
                
                distances = np.linalg.norm(positions[:i] - pos, axis=1)
                if np.all(distances > min_distance):
                    positions[i] = pos
                    break
        
        velocities = np.random.uniform(-0.2, 0.2, (n_particles, 2))  # Reduced for gentler dynamics
        masses = np.ones(n_particles)
        
        # Visual scaling factor - compress distances for display (smaller = more compressed)
        visual_scale = 0.5  # Makes particles appear closer together visually
        
        # Calculate visual positions (scaled relative to center of mass)
        def get_visual_positions():
            com = np.mean(positions, axis=0)  # Center of mass
            visual_pos = com + (positions - com) * visual_scale
            return visual_pos
        
        visual_positions = get_visual_positions()
        
        # Create visual elements
        gold_color = "#FFD700"
        particles = VGroup(*[
            Circle(radius=0.2, color=gold_color, fill_opacity=0.9)
            .set_fill(gold_color)
            .set_stroke(GOLD_E, width=3)
            .move_to([visual_positions[i, 0], visual_positions[i, 1], 0])
            for i in range(n_particles)
        ])
        
        # Create line between particles
        particle_line = Line(
            particles[0].get_center(),
            particles[1].get_center(),
            color=BLUE,
            stroke_width=2
        )
        
        boundary = Square(side_length=box_size, color=WHITE, stroke_width=2)
        boundary.shift(LEFT * 2)  # Shift simulation box to the left
        
        # Create potential energy graph
        graph_axes = Axes(
            x_range=[0.8, 3.5, 0.5],
            y_range=[-1.5, 2, 0.5],
            x_length=3.5,
            y_length=4,
            axis_config={"color": WHITE, "include_numbers": True, "font_size": 20},
            tips=False
        )
        graph_axes.shift(RIGHT * 4 + DOWN * 0.5)
        
        # Plot the parabolic potential curve
        potential_curve = graph_axes.plot(
            lambda r: lennard_jones_potential(r, epsilon, sigma),
            x_range=[0.85, 3.5],
            color=YELLOW
        )
        
        # Labels for the graph
        graph_title = Text("Parabolic Potential", font_size=24).next_to(graph_axes, UP, buff=0.2)
        x_label = Text("r", font_size=20).next_to(graph_axes.x_axis, DOWN, buff=0.2)
        y_label = Text("U(r)", font_size=20).next_to(graph_axes.y_axis, LEFT, buff=0.2)
        
        # Create indicator dot on the graph
        initial_distance = np.linalg.norm(positions[1] - positions[0])
        initial_potential = lennard_jones_potential(initial_distance, epsilon, sigma)
        indicator_dot = Dot(
            graph_axes.c2p(initial_distance, initial_potential),
            color=RED,
            radius=0.08
        )
        
        # Distance label - use a simple Text that we'll update manually
        distance_text = Text(
            f"r = {initial_distance:.2f}",
            font_size=20
        ).next_to(graph_axes, DOWN, buff=0.5)
        
        title = Text("MD Simulation with Parabolic Potential", font_size=28)
        title.to_edge(UP)
        
        # Show setup
        self.add(boundary, title, particles, particle_line)
        self.add(graph_axes, potential_curve, graph_title, x_label, y_label, indicator_dot, distance_text)
        
        # Simulation functions
        def update_particles(mob, dt_val):
            nonlocal positions, velocities
            
            # Calculate forces
            forces = np.zeros((n_particles, 2))
            for i in range(n_particles):
                for j in range(i + 1, n_particles):
                    r_vec = positions[j] - positions[i]
                    force = lennard_jones_force(r_vec, epsilon, sigma)
                    forces[i] += force
                    forces[j] -= force
            
            # Velocity Verlet integration
            accelerations = forces / masses[:, np.newaxis]
            velocities += 0.5 * accelerations * dt_val
            positions += velocities * dt_val
            
            # Boundary conditions
            wall_buffer = 0.3  # Distance from wall to reflect particles
            for i in range(n_particles):
                # Check x-dimension (with shifted center)
                if abs(positions[i, 0] - box_center_x) > box_size/2 - wall_buffer:
                    positions[i, 0] = box_center_x + np.sign(positions[i, 0] - box_center_x) * (box_size/2 - wall_buffer)
                    velocities[i, 0] *= -0.85
                # Check y-dimension (centered at 0)
                if abs(positions[i, 1]) > box_size/2 - wall_buffer:
                    positions[i, 1] = np.sign(positions[i, 1]) * (box_size/2 - wall_buffer)
                    velocities[i, 1] *= -0.85
            
            # Update forces
            forces = np.zeros((n_particles, 2))
            for i in range(n_particles):
                for j in range(i + 1, n_particles):
                    r_vec = positions[j] - positions[i]
                    force = lennard_jones_force(r_vec, epsilon, sigma)
                    forces[i] += force
                    forces[j] -= force
            
            accelerations = forces / masses[:, np.newaxis]
            velocities += 0.5 * accelerations * dt_val
            
            # Update visuals with scaled positions
            visual_positions = get_visual_positions()
            for i, particle in enumerate(mob):
                new_pos = [visual_positions[i, 0], visual_positions[i, 1], 0]
                particle.move_to(new_pos)
        
        def update_line(mob):
            mob.put_start_and_end_on(
                particles[0].get_center(),
                particles[1].get_center()
            )
        
        def update_indicator(mob):
            distance = np.linalg.norm(positions[1] - positions[0])
            # Clamp distance to graph range
            # distance = np.clip(distance, 0, 4)
            potential = lennard_jones_potential(distance, epsilon, sigma)
            # Clamp potential to graph range
            # potential = np.clip(potential, -1.5, 2.0)
            mob.move_to(graph_axes.c2p(distance, potential))
        
        def update_distance_text(mob):
            distance = np.linalg.norm(positions[1] - positions[0])
            new_text = Text(
                f"r = {distance:.2f}",
                font_size=20
            ).next_to(graph_axes, DOWN, buff=0.5)
            mob.become(new_text)
        
        # Add updaters
        particles.add_updater(lambda m, dt: update_particles(m, dt))
        particle_line.add_updater(update_line)
        indicator_dot.add_updater(update_indicator)
        distance_text.add_updater(update_distance_text)
        
        # Run simulation
        self.wait(18)  # Extended duration to observe more interaction
        
        particles.remove_updater(update_particles)
        particle_line.remove_updater(update_line)
        indicator_dot.remove_updater(update_indicator)
        distance_text.remove_updater(update_distance_text)
