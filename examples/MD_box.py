
from manim import *
import numpy as np

#This one gets the potential right but the MD box wrong, let's mix the two working chunks out of both into one

# ============================================================================
# Shared
# ============================================================================
""" Here we are going to define functions for the potential and force on particles inside of a box, we are only considering pair potentials. 

We will only use a simple Morse potential with general parameters

"""

# Morse Potential Parameters
D_e = 1.0      # Depth of the potential well
alpha = 1.5    # Width of the potential (controls how "steep" it is)
r_e = 1.0      # Equilibrium bond distance

def morse_potential(r):
    """
    Morse potential: V(r) = D_e * (1 - exp(-alpha*(r - r_e)))^2
    
    Args:
        r: Distance between two particles
    
    Returns:
        Potential energy at distance r
    """
    return D_e * (1 - np.exp(-alpha * (r - r_e)))**2

def morse_force(r):
    """
    Force derived from Morse potential: F(r) = -dV/dr
    F(r) = 2 * D_e * alpha * (1 - exp(-alpha*(r - r_e))) * exp(-alpha*(r - r_e))
    
    Args:
        r: Distance between two particles
    
    Returns:
        Magnitude of force at distance r (positive means repulsive, negative means attractive)
    """
    exp_term = np.exp(-alpha * (r - r_e))
    return 2 * D_e * alpha * (1 - exp_term) * exp_term


# ============================================================================
# Scenes
# ============================================================================

""" 
Here we will have 2 scenes to build
 1) The particles floating around in a box, and their equations of motion show on the right hand side. The movement of the particles is governed by their initial conditions and the potential energy defined above.

 2) We remove all but two particles and draw a line in between them that updates with their movement. On the right a graph of the coulomb potential shows the current potential energy difference between the two particles as a function of the distance between them.

"""


class ParticlesInBox(Scene):
    def construct(self):
        # Title
        title = Text("Molecular Dynamics: Particles in a Box", font_size=32)
        title.to_edge(UP)
        self.play(Write(title))
        self.wait(0.5)
        
        # Create divider
        divider = Line(ORIGIN, 4*DOWN, color=WHITE, stroke_width=2)
        divider.shift(0.5*UP)
        
        # Left side: Simulation box
        # Right side: Equations of motion
        
        # Create the box
        box_size = 3.0
        box = Square(side_length=box_size, color=WHITE, stroke_width=3)
        box.shift(3*LEFT + 0.3*DOWN)
        
        self.play(Create(box), Create(divider))
        self.wait(0.3)
        
        # Particle parameters
        num_particles = 6
        particle_radius = 0.15
        particle_mass = 1.0
        buffer = 0.5  # Minimum distance from walls and other particles
        
        # Initialize particles with random positions (avoiding overlaps)
        particles = VGroup()
        positions = []
        velocities = []
        
        np.random.seed(42)
        
        # Generate non-overlapping positions
        for i in range(num_particles):
            while True:
                # Random position within box with buffer from walls
                x = np.random.uniform(-box_size/2 + buffer, box_size/2 - buffer)
                y = np.random.uniform(-box_size/2 + buffer, box_size/2 - buffer)
                pos = np.array([x, y, 0])
                
                # Check if far enough from other particles
                valid = True
                for other_pos in positions:
                    if np.linalg.norm(pos - other_pos) < buffer:
                        valid = False
                        break
                
                if valid:
                    positions.append(pos)
                    break
        
        # Create particle objects and assign random velocities
        colors = [RED, BLUE, GREEN, YELLOW, PURPLE, ORANGE]
        for i in range(num_particles):
            particle = Circle(radius=particle_radius, color=colors[i], fill_opacity=1)
            particle.move_to(positions[i] + box.get_center())
            particles.add(particle)
            
            # Random initial velocity
            vx = np.random.uniform(-0.8, 0.8)
            vy = np.random.uniform(-0.8, 0.8)
            velocities.append(np.array([vx, vy, 0]))
        
        self.play(Create(particles), run_time=1)
        self.wait(0.3)
        
        # Show equations of motion on the right
        eom_title = Text("Equations of Motion", font_size=24, color=YELLOW)
        eom_title.to_edge(RIGHT).shift(2*UP)
        
        # Newton's second law
        newton = MathTex(
            r"m\ddot{\mathbf{r}}_i = \mathbf{F}_i",
            font_size=28
        ).next_to(eom_title, DOWN, buff=0.5)
        
        # Force from potential
        force_eq = MathTex(
            r"\mathbf{F}_i = -\sum_{j \neq i} \nabla V(r_{ij})",
            font_size=26
        ).next_to(newton, DOWN, buff=0.4)
        
        # Morse potential
        morse_eq = MathTex(
            r"V(r) = D_e(1-e^{-\alpha(r-r_e)})^2",
            font_size=24
        ).next_to(force_eq, DOWN, buff=0.4)
        
        # Boundary conditions
        boundary = Text("Elastic Wall Collisions", font_size=20, color=BLUE)
        boundary.next_to(morse_eq, DOWN, buff=0.5)
        
        self.play(
            Write(eom_title),
            Write(newton),
            run_time=1
        )
        self.wait(0.3)
        
        self.play(
            Write(force_eq),
            Write(morse_eq),
            run_time=1.5
        )
        self.wait(0.3)
        
        self.play(Write(boundary), run_time=0.8)
        self.wait(0.5)
        
        # Simulate the dynamics
        dt = 1/60  # Time step for integration (frame rate)
        total_time = 10.0
        
        def update_particles(mob, dt_anim):
            nonlocal positions, velocities
            
            # Use fixed dt for physics simulation
            # Calculate forces on each particle
            forces = [np.array([0.0, 0.0, 0.0]) for _ in range(num_particles)]
            
            for i in range(num_particles):
                for j in range(i + 1, num_particles):
                    # Vector from i to j
                    r_vec = positions[j] - positions[i]
                    r = np.linalg.norm(r_vec)
                    
                    if r > 0.01:  # Avoid division by zero
                        # Force magnitude from Morse potential
                        f_mag = morse_force(r)
                        
                        # Force direction (unit vector from i to j)
                        r_hat = r_vec / r
                        
                        # When f_mag > 0 (repulsive), force on i points away from j (opposite to r_hat)
                        # When f_mag < 0 (attractive), force on i points toward j (along r_hat)
                        force_on_i = f_mag * r_hat
                        
                        forces[i] += force_on_i
                        forces[j] -= force_on_i  # Newton's third law
            
            # Update positions and velocities
            for i in range(num_particles):
                # Velocity Verlet integration
                velocities[i] += (forces[i] / particle_mass) * dt
                positions[i] += velocities[i] * dt
                
                # Wall collisions (elastic)
                box_min = -box_size/2
                box_max = box_size/2
                
                if positions[i][0] - particle_radius < box_min:
                    positions[i][0] = box_min + particle_radius
                    velocities[i][0] = abs(velocities[i][0])
                elif positions[i][0] + particle_radius > box_max:
                    positions[i][0] = box_max - particle_radius
                    velocities[i][0] = -abs(velocities[i][0])
                
                if positions[i][1] - particle_radius < box_min:
                    positions[i][1] = box_min + particle_radius
                    velocities[i][1] = abs(velocities[i][1])
                elif positions[i][1] + particle_radius > box_max:
                    positions[i][1] = box_max - particle_radius
                    velocities[i][1] = -abs(velocities[i][1])
                
                # Update particle position in animation
                mob[i].move_to(positions[i] + box.get_center())
        
        # Add updater
        particles.add_updater(update_particles)
        
        # Run simulation
        self.wait(total_time)
        
        # Remove updater
        particles.remove_updater(update_particles)
        
        self.wait(2)


class TwoParticleMorse(Scene):
    def construct(self):
        # Title
        title = Text("Morse Potential: Two-Particle Interaction", font_size=32)
        title.to_edge(UP)
        self.play(Write(title))
        self.wait(0.5)
        
        # Create divider
        divider = Line(ORIGIN, 4*DOWN, color=WHITE, stroke_width=2)
        divider.shift(0.5*UP)
        
        # Left side: Two particles with distance line
        # Right side: Morse potential graph
        
        # Create the box
        box_size = 3.0
        box = Square(side_length=box_size, color=WHITE, stroke_width=3)
        box.shift(3*LEFT + 0.3*DOWN)
        
        self.play(Create(box), Create(divider))
        self.wait(0.3)
        
        # Create two particles
        particle_radius = 0.15
        particle_mass = 1.0
        
        # Initial positions (not at equilibrium)
        pos1 = np.array([-0.8, 0.5, 0])
        pos2 = np.array([0.6, -0.3, 0])
        
        # Initial velocities
        vel1 = np.array([0.3, -0.2, 0])
        vel2 = np.array([-0.2, 0.3, 0])
        
        particle1 = Circle(radius=particle_radius, color=RED, fill_opacity=1)
        particle1.move_to(pos1 + box.get_center())
        
        particle2 = Circle(radius=particle_radius, color=BLUE, fill_opacity=1)
        particle2.move_to(pos2 + box.get_center())
        
        # Distance line between particles
        distance_line = Line(particle1.get_center(), particle2.get_center(), color=YELLOW, stroke_width=3)
        
        # Distance label
        r_label = MathTex("r", font_size=24, color=YELLOW)
        
        self.play(Create(particle1), Create(particle2), run_time=0.8)
        self.play(Create(distance_line), Write(r_label), run_time=0.8)
        self.wait(0.3)
        
        # Create Morse potential graph on the right
        axes = Axes(
            x_range=[0.5, 3.0, 0.5],
            y_range=[-0.2, 1.5, 0.5],
            x_length=3.5,
            y_length=3.0,
            axis_config={"include_tip": True, "font_size": 20},
            tips=True,
        ).shift(3*RIGHT + 0.3*DOWN)
        
        # Labels
        x_label = axes.get_x_axis_label(MathTex("r", font_size=24))
        y_label = axes.get_y_axis_label(MathTex("V(r)", font_size=24), edge=LEFT, direction=LEFT)
        
        # Plot Morse potential
        morse_curve = axes.plot(
            morse_potential,
            x_range=[0.5, 3.0],
            color=GREEN,
            stroke_width=3
        )
        
        # Equilibrium point marker
        eq_point = axes.coords_to_point(r_e, 0)
        eq_dot = Dot(eq_point, color=WHITE, radius=0.06)
        eq_label = MathTex("r_e", font_size=20).next_to(eq_dot, DOWN, buff=0.1)
        
        self.play(
            Create(axes),
            Write(x_label),
            Write(y_label),
            run_time=1
        )
        
        self.play(Create(morse_curve), run_time=1.5)
        self.play(Create(eq_dot), Write(eq_label), run_time=0.8)
        self.wait(0.5)
        
        # Current energy indicator on the graph
        current_r = np.linalg.norm(pos2 - pos1)
        current_V = morse_potential(current_r)
        indicator_dot = Dot(
            axes.coords_to_point(current_r, current_V),
            color=YELLOW,
            radius=0.08
        )
        
        self.play(Create(indicator_dot), run_time=0.5)
        self.wait(0.3)
        
        # Simulate dynamics
        dt = 1/60  # Time step for integration
        total_time = 12.0
        
        positions = [pos1.copy(), pos2.copy()]
        velocities = [vel1.copy(), vel2.copy()]
        
        def update_system(mob, dt_anim):
            nonlocal positions, velocities
            
            # Use fixed dt for physics
            # Calculate force between particles
            r_vec = positions[1] - positions[0]
            r = np.linalg.norm(r_vec)
            
            if r > 0.01:
                f_mag = morse_force(r)
                r_hat = r_vec / r
                
                # When f_mag > 0 (repulsive), force pushes particles apart
                # When f_mag < 0 (attractive), force pulls particles together
                force_on_0 = f_mag * r_hat
                force_on_1 = -f_mag * r_hat
                
                # Update velocities
                velocities[0] += (force_on_0 / particle_mass) * dt
                velocities[1] += (force_on_1 / particle_mass) * dt
                
                # Update positions
                positions[0] += velocities[0] * dt
                positions[1] += velocities[1] * dt
                
                # Wall collisions for both particles
                box_min = -box_size/2
                box_max = box_size/2
                
                for i in range(2):
                    if positions[i][0] - particle_radius < box_min:
                        positions[i][0] = box_min + particle_radius
                        velocities[i][0] = abs(velocities[i][0])
                    elif positions[i][0] + particle_radius > box_max:
                        positions[i][0] = box_max - particle_radius
                        velocities[i][0] = -abs(velocities[i][0])
                    
                    if positions[i][1] - particle_radius < box_min:
                        positions[i][1] = box_min + particle_radius
                        velocities[i][1] = abs(velocities[i][1])
                    elif positions[i][1] + particle_radius > box_max:
                        positions[i][1] = box_max - particle_radius
                        velocities[i][1] = -abs(velocities[i][1])
        
        def update_particle1(mob):
            mob.move_to(positions[0] + box.get_center())
        
        def update_particle2(mob):
            mob.move_to(positions[1] + box.get_center())
        
        def update_distance_line(mob):
            mob.put_start_and_end_on(particle1.get_center(), particle2.get_center())
        
        def update_r_label(mob):
            mob.move_to(distance_line.get_center() + 0.3*UP)
        
        def update_indicator(mob):
            r = np.linalg.norm(positions[1] - positions[0])
            V = morse_potential(r)
            # Clamp r to graph range
            r_clamped = np.clip(r, 0.5, 3.0)
            mob.move_to(axes.coords_to_point(r_clamped, V))
        
        # Add updaters
        particle1.add_updater(update_particle1)
        particle2.add_updater(update_particle2)
        distance_line.add_updater(update_distance_line)
        r_label.add_updater(update_r_label)
        indicator_dot.add_updater(update_indicator)
        
        # Update system physics
        self.add(particle1, particle2)  # Re-add to ensure updaters work
        particle1.add_updater(lambda m, dt: update_system(m, dt))
        
        # Run simulation
        self.wait(total_time)
        
        # Remove updaters
        particle1.remove_updater(update_particle1)
        particle2.remove_updater(update_particle2)
        distance_line.remove_updater(update_distance_line)
        r_label.remove_updater(update_r_label)
        indicator_dot.remove_updater(update_indicator)
        
        self.wait(2)
