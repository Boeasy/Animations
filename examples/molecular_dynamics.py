"""
Molecular Dynamics Simulation with Lennard-Jones Potential
6 gold-colored particles interacting via LJ potential
"""
from manim import *
import numpy as np


# ============================================================================
# Shared
# ============================================================================

def lennard_jones_potential(r, epsilon=1.0, sigma=1.0):
    """
    Calculate Lennard-Jones potential energy.
    
    U(r) = 4ε[(σ/r)^12 - (σ/r)^6]
    
    Args:
        r: Distance between particles (scalar or array)
        epsilon: Depth of potential well (default: 1.0)
        sigma: Distance at which potential is zero (default: 1.0)
    
    Returns:
        Potential energy at distance r
    """
    if np.isscalar(r):
        if r < 0.01:
            return 2.0  # Return large positive value to avoid singularity
    else:
        r = np.maximum(r, 0.01)  # Clamp array values
    
    return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)

def lennard_jones_force(r_vec, epsilon=1.0, sigma=1.0):
    """
    Calculate Lennard-Jones force between two particles.
    
    F = 24ε[(2(σ/r)^13 - (σ/r)^7)] * r_hat / r
    
    Args:
        r_vec: Vector from particle i to particle j (2D numpy array)
        epsilon: Depth of potential well (default: 1.0)
        sigma: Distance at which potential is zero (default: 1.0)
    
    Returns:
        Force vector on particle i due to particle j
    """
    r = np.linalg.norm(r_vec)
    if r < 0.01:  # Avoid division by zero
        return np.zeros(2)
    
    # F = 24 * epsilon * (2 * (sigma/r)^13 - (sigma/r)^7) * r_hat / r
    force_magnitude = 24 * epsilon * (2 * (sigma / r)**13 - (sigma / r)**7) / r
    force = force_magnitude * r_vec / r
    return force


# ============================================================================
# Scenes
# ============================================================================

class MolecularDynamics(Scene):
    def construct(self):
        # Simulation parameters
        n_particles = 6
        dt = 0.05  # time step
        n_steps = 200  # number of simulation steps
        
        # Lennard-Jones parameters (in arbitrary units)
        epsilon = 1.0  # depth of potential well
        sigma = 0.5    # distance at which potential is zero
        
        # Box boundaries (scaled to fit the scene)
        box_size = 5.0
        edge_buffer = 0.4  # Keep particles away from edges
        
        # Initialize particles with random positions and velocities
        # Ensure particles don't start too close to walls or each other
        np.random.seed(42)  # for reproducibility
        positions = np.zeros((n_particles, 2))
        min_distance = 0.8  # Minimum distance between particles
        
        for i in range(n_particles):
            max_attempts = 100
            for attempt in range(max_attempts):
                # Generate position with buffer from edges
                pos = np.random.uniform(
                    -(box_size/2 - edge_buffer), 
                    (box_size/2 - edge_buffer), 
                    2
                )
                
                # Check distance from all existing particles
                if i == 0:
                    positions[i] = pos
                    break
                
                distances = np.linalg.norm(positions[:i] - pos, axis=1)
                if np.all(distances > min_distance):
                    positions[i] = pos
                    break
        
        velocities = np.random.uniform(-0.5, 0.5, (n_particles, 2))
        masses = np.ones(n_particles)
        
        # Create visual elements
        # Gold color: #FFD700
        gold_color = "#FFD700"
        particles = VGroup(*[
            Circle(radius=0.2, color=gold_color, fill_opacity=0.8)
            .set_fill(gold_color)
            .set_stroke(GOLD_E, width=2)
            .move_to([positions[i, 0], positions[i, 1], 0])
            for i in range(n_particles)
        ])
        
        # Create boundary box
        boundary = Square(side_length=box_size, color=WHITE, stroke_width=2)
        
        # Add title
        title = Text("Molecular Dynamics: Lennard-Jones Potential", font_size=28)
        title.to_edge(UP)
        
        # Add particle labels
        subtitle = Text("6 Gold Atoms", font_size=20, color=gold_color)
        subtitle.next_to(title, DOWN)
        
        # Show initial setup
        self.play(Create(boundary), Write(title), Write(subtitle))
        self.play(LaggedStart(*[GrowFromCenter(p) for p in particles], lag_ratio=0.1))
        self.wait(0.5)
        
        # Run molecular dynamics simulation
        def update_particles(mob, dt_val):
            """Update particle positions using velocity Verlet integration"""
            nonlocal positions, velocities
            
            # Calculate forces
            forces = np.zeros((n_particles, 2))
            for i in range(n_particles):
                for j in range(i + 1, n_particles):
                    r_vec = positions[j] - positions[i]
                    force = lennard_jones_force(r_vec)
                    forces[i] += force
                    forces[j] -= force  # Newton's third law
            
            # Update velocities (half step)
            accelerations = forces / masses[:, np.newaxis]
            velocities += 0.5 * accelerations * dt_val
            
            # Update positions
            positions += velocities * dt_val
            
            # Apply boundary conditions (elastic collisions with walls)
            wall_buffer = 0.3  # Distance from wall to reflect particles
            for i in range(n_particles):
                for dim in range(2):
                    if abs(positions[i, dim]) > box_size/2 - wall_buffer:
                        positions[i, dim] = np.sign(positions[i, dim]) * (box_size/2 - wall_buffer)
                        velocities[i, dim] *= -0.85  # slightly inelastic collision
            
            # Update velocities (half step)
            forces = np.zeros((n_particles, 2))
            for i in range(n_particles):
                for j in range(i + 1, n_particles):
                    r_vec = positions[j] - positions[i]
                    force = lennard_jones_force(r_vec, epsilon, sigma)
                    forces[i] += force
                    forces[j] -= force
            
            accelerations = forces / masses[:, np.newaxis]
            velocities += 0.5 * accelerations * dt_val
            
            # Update visual positions
            for i, particle in enumerate(mob):
                particle.move_to([positions[i, 0], positions[i, 1], 0])
        
        # Animate the simulation
        particles.add_updater(lambda m, dt: update_particles(m, dt))
        
        # Run simulation
        self.wait(10)  # 10 seconds of animation
        
        particles.remove_updater(update_particles)
        
        # Fade out
        self.play(
            FadeOut(particles),
            FadeOut(boundary),
            FadeOut(title),
            FadeOut(subtitle)
        )

class MolecularDynamicsWithTrails(Scene):
    """Enhanced version with particle trails to visualize motion"""
    def construct(self):
        # Simulation parameters
        n_particles = 6
        dt = 0.05
        
        # Lennard-Jones parameters
        epsilon = 1.0
        sigma = 1.0
        
        # Box boundaries
        box_size = 5.0
        edge_buffer = 0.4  # Keep particles away from edges
        
        # Initialize particles with proper spacing
        np.random.seed(42)
        positions = np.zeros((n_particles, 2))
        min_distance = 0.8  # Minimum distance between particles
        
        for i in range(n_particles):
            max_attempts = 100
            for attempt in range(max_attempts):
                # Generate position with buffer from edges
                pos = np.random.uniform(
                    -(box_size/2 - edge_buffer), 
                    (box_size/2 - edge_buffer), 
                    2
                )
                
                # Check distance from all existing particles
                if i == 0:
                    positions[i] = pos
                    break
                
                distances = np.linalg.norm(positions[:i] - pos, axis=1)
                if np.all(distances > min_distance):
                    positions[i] = pos
                    break
        
        velocities = np.random.uniform(-0.5, 0.5, (n_particles, 2))
        masses = np.ones(n_particles)
        
        # Create visual elements
        gold_color = "#FFD700"
        particles = VGroup(*[
            Circle(radius=0.2, color=gold_color, fill_opacity=0.9)
            .set_fill(gold_color)
            .set_stroke(GOLD_E, width=3)
            .move_to([positions[i, 0], positions[i, 1], 0])
            for i in range(n_particles)
        ])
        
        # Create trails - initialize with starting positions
        trails = VGroup(*[VMobject() for _ in range(n_particles)])
        for i, trail in enumerate(trails):
            trail.set_stroke(gold_color, width=1, opacity=0.3)
            trail.start_new_path([positions[i, 0], positions[i, 1], 0])
        
        boundary = Square(side_length=box_size, color=WHITE, stroke_width=2)
        
        title = Text("MD Simulation with Particle Trails", font_size=28)
        title.to_edge(UP)
        
        # Show setup
        self.add(boundary, title, particles, trails)
        
        # Simulation functions
        def update_particles(mob, dt_val):
            nonlocal positions, velocities
            
            # Calculate forces
            forces = np.zeros((n_particles, 2))
            for i in range(n_particles):
                for j in range(i + 1, n_particles):
                    r_vec = positions[j] - positions[i]
                    force = lennard_jones_force(r_vec)
                    forces[i] += force
                    forces[j] -= force
            
            # Velocity Verlet integration
            accelerations = forces / masses[:, np.newaxis]
            velocities += 0.5 * accelerations * dt_val
            positions += velocities * dt_val
            
            # Boundary conditions
            wall_buffer = 0.3  # Distance from wall to reflect particles
            for i in range(n_particles):
                for dim in range(2):
                    if abs(positions[i, dim]) > box_size/2 - wall_buffer:
                        positions[i, dim] = np.sign(positions[i, dim]) * (box_size/2 - wall_buffer)
                        velocities[i, dim] *= -0.85
            
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
            
            # Update visuals
            for i, particle in enumerate(mob):
                new_pos = [positions[i, 0], positions[i, 1], 0]
                particle.move_to(new_pos)
        
        def update_trails(mob):
            for i, trail in enumerate(mob):
                particle_pos = particles[i].get_center()
                trail.add_line_to(particle_pos)
        
        # Add updaters
        particles.add_updater(lambda m, dt: update_particles(m, dt))
        trails.add_updater(update_trails)
        
        # Run simulation
        self.wait(12)
        
        particles.remove_updater(update_particles)
        trails.remove_updater(update_trails)

class LennardJonesPotential(Scene):
    """Visualize the Lennard-Jones potential curve"""
    def construct(self):
        # Create axes
        axes = Axes(
            x_range=[0.8, 3, 0.5],
            y_range=[-1.5, 2, 0.5],
            x_length=8,
            y_length=6,
            axis_config={"color": BLUE},
            tips=False,
        )
        
        # Labels
        x_label = axes.get_x_axis_label(r"r / \sigma", direction=DOWN)
        y_label = axes.get_y_axis_label(r"U / \epsilon", direction=LEFT)
        
        # Plot the potential using shared function (with default ε=1, σ=1)
        lj_curve = axes.plot(
            lambda r: lennard_jones_potential(r, epsilon=1.0, sigma=1.0),
            x_range=[0.9, 3],
            color=YELLOW,
            stroke_width=3,
        )        # Annotations
        title = Text("Lennard-Jones Potential", font_size=36)
        title.to_edge(UP)
        
        formula = MathTex(
            r"U(r) = 4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right]",
            font_size=32
        )
        formula.next_to(title, DOWN)
        
        # Mark special points
        equilibrium_r = 2**(1/6)  # minimum of potential
        equilibrium_point = Dot(
            axes.c2p(equilibrium_r, lennard_jones_potential(equilibrium_r, epsilon=1.0, sigma=1.0)),
            color=GREEN,
            radius=0.08
        )
        equilibrium_label = Text("Equilibrium", font_size=20, color=GREEN)
        equilibrium_label.next_to(equilibrium_point, RIGHT)
        
        zero_crossing = Dot(axes.c2p(1, 0), color=RED, radius=0.08)
        zero_label = MathTex(r"r = \sigma", font_size=24, color=RED)
        zero_label.next_to(zero_crossing, UP + RIGHT)
        
        # Repulsive and attractive regions
        repulsive_text = Text("Repulsive", font_size=24, color=RED)
        repulsive_text.move_to(axes.c2p(0.95, 1.5))
        
        attractive_text = Text("Attractive", font_size=24, color=BLUE)
        attractive_text.move_to(axes.c2p(2.2, -0.8))
        
        # Animate
        self.play(Write(title), Write(formula))
        self.wait(0.5)
        self.play(Create(axes), Write(x_label), Write(y_label))
        self.wait(0.5)
        self.play(Create(lj_curve), run_time=2)
        self.wait(0.5)
        
        self.play(
            GrowFromCenter(equilibrium_point),
            Write(equilibrium_label),
        )
        self.wait(0.5)
        
        self.play(
            GrowFromCenter(zero_crossing),
            Write(zero_label),
        )
        self.wait(0.5)
        
        self.play(
            Write(repulsive_text),
            Write(attractive_text),
        )
        self.wait(2)
        
        # Fade out
        self.play(
            *[FadeOut(mob) for mob in self.mobjects]
        )

class TwoParticleCollision(Scene):
    """Two particles colliding with live potential energy visualization"""
    def construct(self):
        # Simulation parameters
        # epsilon = 1.0
        # sigma = 1.0
        # dt = 0.03

        # Simulation parameters
        n_particles = 6
        dt = 0.05  # time step
        n_steps = 200  # number of simulation steps
        
        # Lennard-Jones parameters (in arbitrary units)
        epsilon = 1.0  # depth of potential well
        sigma = 0.5    # distance at which potential is zero
        
        # Box boundaries (scaled to fit the scene)
        box_size = 5.0
        edge_buffer = 0.4  # Keep particles away from edges

        # Create two sections: left for particles, right for LJ curve
        v_line = Line(UP * 3.5, DOWN * 3.5, color=WHITE, stroke_width=1)
        
        # Left side: particle collision area
        particle_area_center = LEFT * 3.5
        boundary = Square(side_length=box_size, color=WHITE, stroke_width=2)
        boundary.move_to(particle_area_center)
        
        # Initialize two particles moving toward each other
        # Particle 1: red, starting left moving right
        particle1_start = np.array([-1.5, 0.0])  # Initial separation of 3.0σ
        velocity1 = np.array([1, 0.0])  # Moving right
        
        # Particle 2: blue, starting right moving left
        particle2_start = np.array([1.5, 0.0])  # Initial separation of 3.0σ
        velocity2 = np.array([-1, 0.0])  # Moving left
        
        # Store as simulation state
        positions = np.array([particle1_start, particle2_start])
        velocities = np.array([velocity1, velocity2])
        
        # Create visual particles
        particle1 = Circle(radius=0.15, color=RED, fill_opacity=0.9)
        particle1.set_fill(RED)
        particle1.set_stroke(RED_E, width=3)
        particle1.move_to(particle_area_center + np.array([positions[0, 0], positions[0, 1], 0]))
        
        particle2 = Circle(radius=0.15, color=BLUE, fill_opacity=0.9)
        particle2.set_fill(BLUE)
        particle2.set_stroke(BLUE_E, width=3)
        particle2.move_to(particle_area_center + np.array([positions[1, 0], positions[1, 1], 0]))
        
        # Right side: Lennard-Jones potential curve
        curve_center = RIGHT * 3.5
        axes = Axes(
            x_range=[0.8, 3, 0.5],
            y_range=[-1.5, 2, 0.5],
            x_length=4.5,
            y_length=5,
            axis_config={"color": BLUE_E, "stroke_width": 1},
            tips=False,
        ).move_to(curve_center)
        
        # Plot curve using shared function
        lj_curve = axes.plot(
            lambda r: lennard_jones_potential(r, epsilon, sigma),
            x_range=[0.9, 3],
            color=YELLOW,
            stroke_width=2,
        )
        
        # Axis labels
        x_label = axes.get_x_axis_label(r"r/\sigma", direction=DOWN, buff=0.2).scale(0.7)
        y_label = axes.get_y_axis_label(r"U/\epsilon", direction=LEFT, buff=0.2).scale(0.7)
        
        # Equilibrium line
        equilibrium_r = 2**(1/6)
        eq_line = DashedLine(
            axes.c2p(equilibrium_r, -1.5),
            axes.c2p(equilibrium_r, lennard_jones_potential(equilibrium_r, epsilon, sigma)),
            color=GREEN,
            stroke_width=1,
            dash_length=0.05,
        )
        
        # Calculate initial distance
        initial_distance = np.linalg.norm(positions[1] - positions[0]) / sigma
        
        # Energy tracking dot on the curve - start at initial positions
        initial_r_normalized = initial_distance
        initial_potential = lennard_jones_potential(initial_r_normalized, epsilon, sigma)
        initial_energy_point = axes.c2p(initial_r_normalized, initial_potential)
        
        energy_dot = Dot(color=GREEN, radius=0.06)
        energy_dot.move_to(initial_energy_point + LEFT * 0.08)
        
        # Distance label - show initial distance
        distance_label = MathTex(r"r = ", f"{initial_distance:.2f}", r"\sigma", font_size=28)
        distance_label.next_to(axes, UP, buff=0.2)
        
        # Titles
        particle_title = Text("Particle Collision", font_size=24)
        particle_title.next_to(boundary, UP, buff=0.3)
        
        energy_title = Text("Potential Energy", font_size=24)
        energy_title.next_to(axes, UP, buff=0.5)
        
        # Show setup (all at once)
        self.play(
            Create(boundary),
            Create(v_line),
            Write(particle_title),
        )
        self.play(
            GrowFromCenter(particle1),
            GrowFromCenter(particle2),
        )        
        self.play(
            Create(axes),
            Write(x_label),
            Write(y_label),
            Write(energy_title),
        )
        self.play(Create(lj_curve), Create(eq_line))
        self.play(
            GrowFromCenter(energy_dot),
            Write(distance_label),
        )
        self.wait(0.3)
        
        # Simulation update function
        def update_system(dt_val):
            nonlocal positions, velocities
            
            # Calculate distance and force between particles
            r_vec = positions[1] - positions[0]
            r_distance = np.linalg.norm(r_vec)
            
            # Calculate forces
            force = lennard_jones_force(r_vec, epsilon, sigma)
            forces = np.array([force, -force])  # Newton's third law
            
            # Velocity Verlet integration
            masses = np.array([1.0, 1.0])
            accelerations = forces / masses[:, np.newaxis]
            velocities += 0.5 * accelerations * dt_val
            positions += velocities * dt_val
            
            # Recalculate forces
            r_vec = positions[1] - positions[0]
            r_distance = np.linalg.norm(r_vec)
            force = lennard_jones_force(r_vec, epsilon, sigma)
            forces = np.array([force, -force])
            accelerations = forces / masses[:, np.newaxis]
            velocities += 0.5 * accelerations * dt_val
            
            # Update particle positions
            particle1.move_to(particle_area_center + np.array([positions[0, 0], positions[0, 1], 0]))
            particle2.move_to(particle_area_center + np.array([positions[1, 0], positions[1, 1], 0]))
            
            # Update energy dots on curve
            r_normalized = r_distance / sigma
            r_normalized = max(0.9, min(r_normalized, 2.9))  # Clamp to plot range
            
            potential_energy = lennard_jones_potential(r_normalized, epsilon, sigma)
            potential_energy = max(-1.4, min(potential_energy, 1.9))  # Clamp to plot range
            
            energy_point = axes.c2p(r_normalized, potential_energy)
            energy_dot1.move_to(energy_point + LEFT * 0.08)
            energy_dot2.move_to(energy_point + RIGHT * 0.08)
            
            # Update distance label
            new_distance_label = MathTex(
                r"r = ", 
                f"{r_normalized:.2f}", 
                r"\sigma", 
                font_size=28
            )
            new_distance_label.move_to(distance_label.get_center())
            distance_label.become(new_distance_label)
        
        # Add updater that uses frame dt
        self.add(particle1, particle2, energy_dot1, energy_dot2, distance_label)
        
        # Use always_redraw or add updater to scene
        # Manim Scene updaters receive only the frame dt.
        def scene_updater(dt):
            update_system(dt)

        # Attach updater to the scene so it runs every frame
        self.add_updater(scene_updater)
        
        # Run the collision simulation
        self.wait(15)  # Longer time to see the full collision
        
        self.remove_updater(scene_updater)
        self.wait(1.5)
        
        # Fade out
        self.play(*[FadeOut(mob) for mob in self.mobjects])
