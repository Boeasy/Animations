
from manim import *
import numpy as np

#This one gets the potential right but the MD box wrong, let's mix the two working chunks out of both into one

# ============================================================================
# Shared
# ============================================================================
""" Here we are going to define functions for the potential and force on particles inside of a box, we are only considering pair potentials. 

We will only use a simple Morse potential with general parameters

"""

# Morse Potential Parameters (1.0 , 1.5, 1.0)
D_e = 1.5      # Depth of the potential well
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
    
    #add the safety against division by 0 in here...
    if r < 0.01:
        return 2.0  # Return large positive value to avoid singularity
    elif r < r_e:
        exp_term = np.exp(-alpha * (r - r_e))
        return 2 * D_e * alpha * (1 - exp_term) * exp_term
    else:
        exp_term = np.exp(-alpha * (r - r_e))
        return 5 * D_e * alpha * (1 - exp_term) * exp_term

# ============================================================================
# Scenes
# ============================================================================

""" 
Here we will have 2 scenes to build
 1) The particles floating around in a box, and their equations of motion show on the right hand side. The movement of the particles is governed by their initial conditions and the potential energy defined above.

 2) We remove all but two particles and draw a line in between them that updates with their movement. On the right a graph of the coulomb potential shows the current potential energy difference between the two particles as a function of the distance between them.

"""


class ParticlesInBox(Scene):
    """
    This scene will put a number of particles into a box on the left, with general equations of motion shown on the right. 
    
    Don't put a divider in. 

    At every time step, for each particle in the box, calculate the force on that particle due to all the other particles, and then update the motion of that particle.

    Use vertlet integration.

    forces should be calculated using the morse_force function

    Start with 6 particles, but set up the scene generally to work with a num_particles variable.

    There should be a buffer between the particles and the edge such that the circles don't appear to leave the bounding box.

    """
    def construct(self):
        # Simulation parameters
        num_particles = 6
        dt = 0.05  # time step
        gold_color = GOLD

        # Box boundaries
        box_size = 5.0
        edge_buffer = 0.4  # Keep particles away from edges
        
        # Initialize particles with proper spacing
        np.random.seed(42)
        positions = np.zeros((num_particles, 2))
        min_distance = 0.8  # Minimum distance between particles
        
        for i in range(num_particles):
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
        
        velocities = np.random.uniform(-0.5, 0.5, (num_particles, 2))
        masses = np.ones(num_particles)
        
        # Create boundary box (on the left) FIRST
        boundary = Square(side_length=box_size, color=WHITE, stroke_width=2)
        boundary.shift(3*LEFT + 0.3*DOWN)
        
        # Create visual elements - all gold, positioned INSIDE the box
        particles = VGroup(*[
            Circle(radius=0.2, color=gold_color, fill_opacity=0.9)
            .set_fill(gold_color)
            .set_stroke(GOLD_E, width=3)
            .move_to([positions[i, 0], positions[i, 1], 0] + boundary.get_center())
            for i in range(num_particles)
        ])
        
        # Add title
        title = Text("Molecular Dynamics: Example System", font_size=28)
        title.to_edge(UP)
        
        # Show equations of motion on the right
        eom_title = Text("Equations of Motion", font_size=24, color=YELLOW)
        eom_title.shift(2*UP, 3*RIGHT)
        
        # Newton's second law
        newton = MathTex(
            r"\sum \vec{F}=m\vec{a}",
            font_size=28
        ).next_to(eom_title, DOWN, buff=0.5)
        
        # Acceleration
        acceleration = MathTex(
            r"\vec{a_{i}} = \frac{d^{2}\vec{r}_{i}}{dt^{2}} ",
            font_size=26
        ).next_to(newton, DOWN, buff=0.4)

        # Force Vector
        force_vec = MathTex(
            r"\vec{F_{i}}=-\nabla_{i}U(\vec{r_{1}},\vec{r_{2}},\vec{r_{3}},\vec{r_{4}}, \vec{r_{5}}, \vec{r_{6}})",
            font_size=26
        ).next_to(acceleration, DOWN, buff=0.4)

        # 2nd Order Differential Equations
        order2 = MathTex(
            r"m_{i} \frac{d^{2}\vec{r_{i}}}{dt^{2}}=-\nabla_{i}U(\vec{r_{1}},\dots, \vec{r_{6}})",
            font_size=26
        ).next_to(force_vec, DOWN, buff=0.4)

    
        
        # # Morse potential
        # morse_eq = MathTex(
        #     r"V(r) = D_e(1-e^{-\alpha(r-r_e)})^2",
        #     font_size=24
        # ).next_to(order2, DOWN, buff=0.4)
        
        # Boundary conditions
        # boundary_text = Text("Elastic Wall Collisions", font_size=20, color=BLUE)
        # boundary_text.next_to(morse_eq, DOWN, buff=0.5)
        
        # Show initial setup
        self.add(boundary, title, particles)
        
        self.play(
            Write(eom_title),
            Write(newton),
            run_time=1
        )
        self.wait(0.3)
        
        self.play(
            Write(acceleration),
            Write(force_vec),
            Write(order2),
            # Write(morse_eq),
            run_time=0.5
        )
        
        # Run molecular dynamics simulation
        def update_particles(mob, dt_val):
            """Update particle positions using velocity Verlet integration"""
            nonlocal positions, velocities
            
            # Calculate forces on each particle
            forces = np.zeros((num_particles, 2))
            for i in range(num_particles):
                for j in range(i + 1, num_particles):
                    # Vector from i to j
                    r_vec = positions[j] - positions[i]
                    r = np.linalg.norm(r_vec)
                    f_mag = morse_force(r)
                    r_hat = r_vec / r
                    force = f_mag * r_hat
                    forces[i] += force
                    forces[j] -= force  # Newton's third law
            
            # Velocity Verlet integration: Update velocities (half step)
            accelerations = forces / masses[:, np.newaxis]
            velocities += 0.5 * accelerations * dt_val
            positions += velocities * dt_val
            
            # Apply boundary conditions (elastic collisions with walls)
            wall_buffer = 0.3  # Distance from wall to reflect particles
            for i in range(num_particles):
                for dim in range(2):
                    if abs(positions[i, dim]) > box_size/2 - wall_buffer:
                        positions[i, dim] = np.sign(positions[i, dim]) * (box_size/2 - wall_buffer)
                        velocities[i, dim] *= -0.85  # slightly inelastic collision
            
            # Update forces
            forces = np.zeros((num_particles, 2))
            for i in range(num_particles):
                for j in range(i + 1, num_particles):
                    r_vec = positions[j] - positions[i]
                    r = np.linalg.norm(r_vec)
                    f_mag = morse_force(r)
                    r_hat = r_vec / r
                    force = f_mag * r_hat
                    forces[i] += force
                    forces[j] -= force
            
            accelerations = forces / masses[:, np.newaxis]
            velocities += 0.5 * accelerations * dt_val
            
            # Update visual positions
            # for i, particle in enumerate(mob):
            #     particle.move_to([positions[i, 0], positions[i, 1], 0] + boundary.get_center())


            #.move_to([positions[i, 0], positions[i, 1], 0] + boundary.get_center())

            for i, particle in enumerate(mob):
                particle.move_to([positions[i, 0], positions[i, 1], 0] + boundary.get_center())
        
        # Animate the simulation
        particles.add_updater(lambda m, dt: update_particles(m, dt))
        
        # Run simulation
        self.wait(30)
        
        particles.remove_updater(update_particles)
        
        self.wait(2)


class TwoParticleMorse(Scene):
    """ Need to fix a few things.
    
    Change to r_i,j and make it follow the yellow line.

    start the animation out by putting a label r_i and r_j on the left and right particles, it can fade out and turn into r_i,j afterwards on the line.

    animation needs to run longer, also could raise the starting velocity of the particles.

    Under the title, write the morse potential equation, in terms of U(r_i,j), 

    maybe under the box write, "F = -grad U" and under the graph write the potential eq? 
    
    
    """
    def construct(self):
        # Title
        title = Text("Two-Particle Interaction: Morse Potential", font_size=32)
        title.to_edge(UP)
        self.play(Write(title))
        self.wait(0.5)
        
        # Create the box
        box_size = 3.0
        box = Square(side_length=box_size, color=WHITE, stroke_width=3)
        box.shift(3*LEFT + 0.3*DOWN)
        
        self.play(Create(box))
        self.wait(0.3)
        
        # Create two particles
        particle_radius = 0.15
        particle_mass = 1.0
        
        # Initial positions (at equilibrium)
        pos1 = np.array([-0.2, 0.5, 0])
        pos2 = np.array([0.8, 0.5, 0])
        
        # Initial velocities (raised)
        vel1 = np.array([0.2, -0.5, 0])
        vel2 = np.array([0.0, 2.0, 0])
        
        particle1 = Circle(radius=particle_radius, color=GOLD_A, fill_opacity=1)
        particle1.move_to(pos1 + box.get_center())
        
        particle2 = Circle(radius=particle_radius, color=GOLD_B, fill_opacity=1)
        particle2.move_to(pos2 + box.get_center())

        # position vector labels
        r_1 = MathTex("\vec{r_{i}}", font_size=24, color=YELLOW).move_to(particle1.get_center() + 0.3*UP)
        r_2 = MathTex("\vec{r_{j}}", font_size=24, color=YELLOW).move_to(particle2.get_center() + 0.3*UP)  

        # Distance line between particles
        distance_line = Line(particle1.get_center(), particle2.get_center(), color=RED, stroke_width=3)
        
        # Distance label
        r_label = MathTex("r_{i,j}", font_size=24, color=YELLOW).move_to(distance_line.get_center() + 0.3*UP)


        
        self.play(Create(particle1), Create(particle2), Write(r_1), Write(r_2), run_time=0.8)
        self.wait(1)
        self.play(FadeOut(r_1), FadeOut(r_2))
        self.play(Create(distance_line), Write(r_label), run_time=0.8)
        # Try to put the circles in front of the line
        self.add(particle1, particle2)
        self.wait(0.3)
        
        # Create Morse potential graph on the right
        axes = Axes(
            x_range=[0.5, 3.0, 0.5],
            y_range=[-0.2, 2.0, 0.5],
            x_length=3.5,
            y_length=3.0,
            axis_config={"include_tip": True, "font_size": 20},
            tips=True,
        ).shift(3*RIGHT + 0.3*DOWN)
        
        # Labels
        x_label = axes.get_x_axis_label(MathTex("r_{i,j}", font_size=24))
        y_label = axes.get_y_axis_label(MathTex("U(r_{i,j})", font_size=24), edge=LEFT, direction=LEFT)
        
        # Plot Morse potential
        morse_curve = axes.plot(
            morse_potential,
            x_range=[0.5, 3.0],
            color=GREEN,
            stroke_width=3
        )
        
        # Equilibrium point marker
        eq_point = axes.coords_to_point(r_e, 0)
        eq_dot = Dot(eq_point, color=YELLOW, radius=0.06, fill_opacity=1)

        # Write in the Equations for force and PE

        # Morse potential
        morse_eq = MathTex(
            r"U(r_{i,j}) = D_e(1-e^{-\alpha(r_{i,j}-r_{eq})})^2",
            font_size=36
        ).next_to(axes, UP, buff=0.5)

        # Force due to potential
        force_eq = MathTex(
            r"\vec{F_{i}}=-\nabla_{i}U(r_{i,j})",
            font_size=36
        ).next_to(box, UP, buff=0.5)
        
        self.play(
            Create(axes),
            Write(x_label),
            Write(y_label),
            Write(morse_eq),
            Write(force_eq),
            run_time=0.5
        )
        
        self.play(Create(morse_curve), run_time=1.0)
        self.play(Create(eq_dot), run_time=0.2)
        self.wait(2.0)
        
        # Current energy indicator on the graph
        current_r = np.linalg.norm(pos2 - pos1)
        current_V = morse_potential(current_r)
        indicator_dot = Dot(
            axes.coords_to_point(current_r, current_V),
            color=YELLOW,
            radius=0.08
        )
        
        self.play(
            AnimationGroup(
            FadeOut(eq_dot),
            Create(indicator_dot), run_time=0.5))

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
        self.add(particle1, particle2, r_label)  # Re-add to ensure updaters work
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
