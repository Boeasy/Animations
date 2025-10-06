"""this file will attempt to explain ergodicity in terms of molecular dynamics"""
# Use the conga line example from comp sim of liquids if possible

from manim import *
import numpy as np


class PhaseSpace(Scene):
    def construct(self):
        # Title
        title = Text("Ergodicity: Time Average = Ensemble Average", font_size=32)
        title.to_edge(UP)
        self.play(Write(title))
        self.wait(0.5)
        
        # Parameters for particle system
        n_particles = 40
        box_size = 3.0
        particle_radius = 0.08
        min_distance = 0.3
        
        # Generate random particle positions
        np.random.seed(42)
        positions = []
        velocities = []
        
        max_attempts = 1000
        while len(positions) < n_particles and max_attempts > 0:
            pos = np.random.uniform(-box_size/2, box_size/2, 2)
            
            too_close = False
            for existing_pos in positions:
                if np.linalg.norm(pos - existing_pos) < min_distance:
                    too_close = True
                    break
            
            if not too_close:
                positions.append(pos)
                # Random velocities
                vel = np.random.uniform(-0.5, 0.5, 2)
                velocities.append(vel)
            max_attempts -= 1
        
        # Create box
        box = Square(side_length=box_size, color=WHITE, stroke_width=2)
        box.shift(2.5*LEFT)
        
        # Create particles
        particles = VGroup()
        for pos in positions:
            particle = Circle(
                radius=particle_radius,
                color=BLUE,
                fill_opacity=0.8,
                stroke_width=1
            )
            particle.move_to([pos[0] + 2.5*LEFT[0], pos[1], 0])
            particles.add(particle)
        
        # Show the system
        self.play(Create(box), run_time=0.8)
        self.play(Create(particles), run_time=1.2)
        self.wait(0.5)
        
        # Select one particle to track
        tracked_idx = len(particles) // 2
        tracked_particle = particles[tracked_idx]
        
        # Highlight tracked particle
        highlight = Circle(
            radius=particle_radius * 2,
            color=YELLOW,
            stroke_width=4
        ).move_to(tracked_particle.get_center())
        
        self.play(
            Create(highlight),
            tracked_particle.animate.set_color(YELLOW).set_fill(opacity=1),
            *[particles[i].animate.set_color(GRAY).set_opacity(0.3) 
              for i in range(len(particles)) if i != tracked_idx],
            run_time=1.5
        )
        self.wait(0.5)
        
        # Add velocity vector display on the right
        velocity_box = Rectangle(width=2.5, height=1.5, color=WHITE, stroke_width=2)
        velocity_box.to_edge(RIGHT).shift(UP)
        
        velocity_label = Text("Velocity", font_size=20, color=YELLOW)
        velocity_label.next_to(velocity_box, UP, buff=0.2)
        
        self.play(
            Create(velocity_box),
            Write(velocity_label),
            run_time=1
        )
        
        # Create velocity vector display
        vel_arrow = Arrow(
            ORIGIN,
            RIGHT * 0.5,
            color=RED,
            buff=0,
            stroke_width=4
        ).move_to(velocity_box.get_center())
        
        vel_text = MathTex("v = ", font_size=24, color=RED)
        vel_text.next_to(velocity_box, DOWN, buff=0.2)
        
        self.play(GrowArrow(vel_arrow), Write(vel_text))
        self.wait(0.5)
        
        # Animate particle motion
        time_tracker = ValueTracker(0)
        dt = 0.1
        
        # Store current positions and velocities
        current_positions = [np.array(positions[i]) for i in range(len(positions))]
        current_velocities = [np.array(velocities[i]) for i in range(len(velocities))]
        
        # Track phase space trajectory
        phase_trajectory = []
        
        def update_particles(mob):
            t = time_tracker.get_value()
            
            for i in range(len(mob)):
                # Update position
                current_positions[i] += current_velocities[i] * dt
                
                # Bounce off walls
                for dim in range(2):
                    if abs(current_positions[i][dim]) > box_size/2:
                        current_positions[i][dim] = np.clip(current_positions[i][dim], -box_size/2, box_size/2)
                        current_velocities[i][dim] *= -1
                
                # Update visual position
                mob[i].move_to([current_positions[i][0] - 2.5, current_positions[i][1], 0])
            
            # Track phase space for highlighted particle
            if len(phase_trajectory) < 200:
                phase_trajectory.append({
                    'pos': current_positions[tracked_idx].copy(),
                    'vel': current_velocities[tracked_idx].copy()
                })
        
        def update_highlight(mob):
            mob.move_to(particles[tracked_idx].get_center())
        
        def update_velocity_arrow(mob):
            vel = current_velocities[tracked_idx]
            speed = np.linalg.norm(vel)
            if speed > 0.01:
                direction = vel / speed
                arrow_length = min(speed * 1.5, 1.0)
                mob.put_start_and_end_on(
                    velocity_box.get_center(),
                    velocity_box.get_center() + np.array([direction[0] * arrow_length, direction[1] * arrow_length, 0])
                )
        
        particles.add_updater(update_particles)
        highlight.add_updater(update_highlight)
        vel_arrow.add_updater(update_velocity_arrow)
        
        # Run simulation
        self.play(
            time_tracker.animate.set_value(5),
            run_time=5,
            rate_func=linear
        )
        
        # Remove updaters
        particles.remove_updater(update_particles)
        highlight.remove_updater(update_highlight)
        vel_arrow.remove_updater(update_velocity_arrow)
        
        self.wait(0.5)
        
        # Transition to phase space
        phase_label = Text("Phase Space Trajectory", font_size=24, color=GREEN)
        phase_label.next_to(velocity_box, DOWN, buff=1.5)
        
        self.play(Write(phase_label))
        self.wait(0.5)
        
        # Create phase space axes
        phase_axes = Axes(
            x_range=[-box_size/2, box_size/2, 1],
            y_range=[-0.6, 0.6, 0.2],
            x_length=3,
            y_length=2,
            axis_config={"include_tip": True, "include_numbers": False},
        ).to_edge(RIGHT).shift(DOWN)
        
        x_label = phase_axes.get_x_axis_label(Text("Position x", font_size=16))
        y_label = phase_axes.get_y_axis_label(Text("Velocity vₓ", font_size=16), edge=LEFT, direction=LEFT)
        
        self.play(
            Create(phase_axes),
            Write(x_label),
            Write(y_label),
            run_time=1.5
        )
        
        # Plot phase space trajectory
        if len(phase_trajectory) > 1:
            phase_points = VGroup()
            for data in phase_trajectory[::2]:  # Sample every other point
                point = Dot(
                    phase_axes.c2p(data['pos'][0], data['vel'][0]),
                    radius=0.02,
                    color=GREEN,
                    fill_opacity=0.6
                )
                phase_points.add(point)
            
            self.play(Create(phase_points), run_time=2)
            self.wait(1)
        
        # Add explanation
        explanation = VGroup(
            Text("Over time, the particle", font_size=18),
            Text("explores the accessible", font_size=18),
            Text("phase space uniformly", font_size=18, color=GREEN)
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.15)
        explanation.to_edge(DOWN).shift(0.3*UP)
        
        self.play(Write(explanation), run_time=2)
        self.wait(3)

class CongaLine(Scene):
    # -  -  -  -  -  -  -  -  -  -
    # D  1  2  3  4  5  6  7  8  -
    # -  58 15 14 13 12 11 10 9  -
    # -  57 16 17 18 19 20 21 22 -
    # -  56 29 28 27 26 25 24 23 -
    # -  55 30 31 32 33 34 35 36 -
    # -  54 43 42 41 40 39 38 37 -
    # -  53 44 45 46 47 59 60 61 -
    # -  52 51 50 49 48 64 63 62 D*
    # -  -  -  -  -  -  -  -  -  -

    #  Above diagram shows how the grid should be set up. D is the first detector position, D* is the second.
    # 1-58 are part of the 'large loop' and 59-64 are the 'small loop'
    # there are vertices between each consecutive circle, ie "1 - 2", and 58 connects to 1, and 64 connects to 59. 
    # Each circle should be a random color R G B in the large loop
    # The small loop they will only be R
    # The animation will be first, the detector triangle is at D and each circle moves to the position of the next consecutive circle, until back to its original position. 
    # The detector will keep count of the color of each circle by a count, R G or B and then show a circle that is the average color of all of these circles.
    # It will then do this same process for the small loop with detector at position D*
    
    def construct(self):
        # Title
        title = Text("Ergodicity: The Conga Line Analogy", font_size=36)
        title.to_edge(UP)
        self.play(Write(title))
        self.wait(0.5)
        
        # Parameters
        circle_radius = 0.15
        spacing = 0.6
        
        # Define the grid layout based on the diagram
        # The grid is 10x10, with row 0 and col 0,9 being empty borders
        # Map from circle number to (row, col) position
        grid_layout = {}
        
        # Row 1: D 1 2 3 4 5 6 7 8 -
        positions_row1 = [1, 2, 3, 4, 5, 6, 7, 8]
        for i, num in enumerate(positions_row1):
            grid_layout[num] = (1, i + 1)
        
        # Row 2: - 58 15 14 13 12 11 10 9 -
        positions_row2 = [58, 15, 14, 13, 12, 11, 10, 9]
        for i, num in enumerate(positions_row2):
            grid_layout[num] = (2, i + 1)
        
        # Row 3: - 57 16 17 18 19 20 21 22 -
        positions_row3 = [57, 16, 17, 18, 19, 20, 21, 22]
        for i, num in enumerate(positions_row3):
            grid_layout[num] = (3, i + 1)
        
        # Row 4: - 56 29 28 27 26 25 24 23 -
        positions_row4 = [56, 29, 28, 27, 26, 25, 24, 23]
        for i, num in enumerate(positions_row4):
            grid_layout[num] = (4, i + 1)
        
        # Row 5: - 55 30 31 32 33 34 35 36 -
        positions_row5 = [55, 30, 31, 32, 33, 34, 35, 36]
        for i, num in enumerate(positions_row5):
            grid_layout[num] = (5, i + 1)
        
        # Row 6: - 54 43 42 41 40 39 38 37 -
        positions_row6 = [54, 43, 42, 41, 40, 39, 38, 37]
        for i, num in enumerate(positions_row6):
            grid_layout[num] = (6, i + 1)
        
        # Row 7: - 53 44 45 46 47 59 60 61 -
        positions_row7 = [53, 44, 45, 46, 47, 59, 60, 61]
        for i, num in enumerate(positions_row7):
            grid_layout[num] = (7, i + 1)
        
        # Row 8: - 52 51 50 49 48 64 63 62 -
        positions_row8 = [52, 51, 50, 49, 48, 64, 63, 62]
        for i, num in enumerate(positions_row8):
            grid_layout[num] = (8, i + 1)
        
        # Detector positions
        detector_D = (1, 0)  # Left of position 1
        detector_D_star = (8, 9)  # Right of position 62
        
        # Convert grid positions to coordinates
        def grid_to_coord(row, col):
            x = (col - 4.5) * spacing
            y = (4.5 - row) * spacing
            return np.array([x, y, 0])
        
        # Create position map for all circles
        circle_positions = {}
        for num, (row, col) in grid_layout.items():
            circle_positions[num] = grid_to_coord(row, col)
        
        # Define loop order
        large_loop_order = list(range(1, 59))  # 1 to 58
        small_loop_order = list(range(59, 65))  # 59 to 64
        
        # Assign colors
        np.random.seed(42)
        colors_list = [RED, GREEN, BLUE]
        circle_colors = {}
        
        # Large loop: random RGB
        for num in large_loop_order:
            circle_colors[num] = colors_list[np.random.randint(0, 3)]
        
        # Small loop: all RED
        for num in small_loop_order:
            circle_colors[num] = RED
        
        # Create circles
        circles = {}
        for num in range(1, 65):
            circle = Circle(
                radius=circle_radius,
                color=circle_colors[num],
                fill_opacity=0.8,
                stroke_width=2,
                stroke_color=WHITE
            )
            circle.move_to(circle_positions[num])
            circles[num] = circle
        
        # Create all circle objects
        all_circles = VGroup(*[circles[i] for i in range(1, 65)])
        
        # Show circles
        self.play(Create(all_circles), run_time=2)
        self.wait(0.5)
        
        # Create connections for large loop
        large_loop_lines = VGroup()
        for i in range(len(large_loop_order)):
            num1 = large_loop_order[i]
            num2 = large_loop_order[(i + 1) % len(large_loop_order)]
            line = Line(
                circles[num1].get_center(),
                circles[num2].get_center(),
                color=YELLOW,
                stroke_width=3
            )
            large_loop_lines.add(line)
        
        # Create connections for small loop
        small_loop_lines = VGroup()
        for i in range(len(small_loop_order)):
            num1 = small_loop_order[i]
            num2 = small_loop_order[(i + 1) % len(small_loop_order)]
            line = Line(
                circles[num1].get_center(),
                circles[num2].get_center(),
                color=ORANGE,
                stroke_width=3
            )
            small_loop_lines.add(line)
        
        # Show connections
        self.play(
            Create(large_loop_lines),
            Create(small_loop_lines),
            run_time=2
        )
        self.wait(0.5)
        
        # Add labels
        large_label = Text("Ergodic System (58 circles)", font_size=20, color=YELLOW)
        large_label.to_edge(LEFT).shift(UP * 3)
        
        small_label = Text("Non-ergodic System (6 circles)", font_size=20, color=ORANGE)
        small_label.next_to(large_label, DOWN, aligned_edge=LEFT)
        
        self.play(Write(large_label), Write(small_label), run_time=1.5)
        self.wait(0.5)
        
        # Create detector at position D
        detector_pos = grid_to_coord(*detector_D)
        detector = Triangle(color=WHITE, fill_opacity=1)
        detector.scale(0.2)
        detector.rotate(3*PI/2)  # Point right
        detector.move_to(detector_pos)
        
        detector_label = Text("Detector", font_size=18, color=WHITE)
        detector_label.next_to(detector, LEFT, buff=0.2)
        
        self.play(Create(detector), Write(detector_label), run_time=1)
        self.wait(0.5)
        
        # Animate large loop: circles move through their positions
        
                # Animate large loop: circles move through their positions
        color_counts = {RED: 0, GREEN: 0, BLUE: 0}
        
        # Add updater for lines to follow circles
        def update_large_lines(mob):
            for i in range(len(large_loop_order)):
                num1 = large_loop_order[i]
                num2 = large_loop_order[(i + 1) % len(large_loop_order)]
                mob[i].put_start_and_end_on(
                    circles[num1].get_center(),
                    circles[num2].get_center()
                )
        
        large_loop_lines.add_updater(update_large_lines)
        
        # Create a smooth path for each circle to follow
        # Each circle follows the loop path continuously
        total_time = 12  # Total time for one complete loop
        time_tracker = ValueTracker(0)
        
        # Store original positions of all circles
        original_positions = {num: circle_positions[num].copy() for num in large_loop_order}
        
        # Create list of all position coordinates in order
        path_positions = [original_positions[num] for num in large_loop_order]
        
        # Track which circles we've counted
        circles_counted = set()
        
        # Highlight circle 1 throughout, used 4 on last animation, need it bigger
        circles[1].set_stroke(width=10, color=YELLOW)
        
        def update_circles(mob, dt):
            t = time_tracker.get_value()
            total_positions = len(large_loop_order)
            
            for i, circle_num in enumerate(large_loop_order):
                # Calculate where this circle should be along the path
                # Each circle starts at position i and moves forward
                progress = (t + i) % total_positions
                
                # Get the two positions to interpolate between
                pos_index = int(progress)
                next_pos_index = (pos_index + 1) % total_positions
                alpha = progress - pos_index
                
                # Interpolate between positions
                current_pos = path_positions[pos_index]
                next_pos = path_positions[next_pos_index]
                interpolated_pos = current_pos * (1 - alpha) + next_pos * alpha
                
                circles[circle_num].move_to(interpolated_pos)
                
                # Check if this circle is currently at position 1 (detector)
                # Position 1 is index 0 in large_loop_order
                circle_at_detector_progress = t % total_positions
                if 0 <= circle_at_detector_progress < 1 and i == 0:
                    # Circle that started at position matching current time step
                    shifted_i = int(t) % total_positions
                    counting_circle = large_loop_order[shifted_i]
                    if counting_circle not in circles_counted:
                        circles_counted.add(counting_circle)
                        color = circle_colors[counting_circle]
                        color_counts[color] += 1
                        
                        # Change detector color to match last circle that passed ... this doesn't seem to be working correctly
                        # detector.set_color(color)
        
        # Add updater to all circles
        for circle_num in large_loop_order:
            circles[circle_num].add_updater(update_circles)
        
        # Animate the time tracker for one complete loop
        self.play(
            time_tracker.animate.set_value(len(large_loop_order)),
            run_time=total_time,
            rate_func=linear
        )
        
        # Remove updaters
        for circle_num in large_loop_order:
            circles[circle_num].remove_updater(update_circles)
        large_loop_lines.remove_updater(update_large_lines)
        
        # Reset circle 1 stroke
        circles[1].set_stroke(width=2)
        self.wait(1)
        
        # Calculate average color for large loop
        total = len(large_loop_order)
        avg_color = [
            color_counts[RED] / total,
            color_counts[GREEN] / total,
            color_counts[BLUE] / total
        ]
        
        # Create average color circle
        avg_circle = Circle(radius=0.6, fill_opacity=0.8, stroke_width=3, stroke_color=WHITE)
        avg_circle.set_color(rgb_to_color(avg_color))
        avg_circle.to_edge(RIGHT).shift(UP * 1.5)
        
        avg_label = Text("Average Color\n(Large Loop)", font_size=18)
        avg_label.next_to(avg_circle, DOWN, buff=0.3)
        
        self.play(Create(avg_circle), Write(avg_label), run_time=1.5)
        self.wait(2)
        
        # Transition to small loop
        # Move detector to D* position
        detector_star_pos = grid_to_coord(*detector_D_star)
        
        self.play(
            detector.animate.move_to(detector_star_pos).rotate(PI).set_color(WHITE),  # Face left and reset color
            detector_label.animate.next_to(detector_star_pos + RIGHT * 0.5, RIGHT, buff=0.2),
            run_time=2
        )
        self.wait(0.5)
        
        # Animate small loop
        color_counts_small = {RED: 0, GREEN: 0, BLUE: 0}
        
        # Add updater for small loop lines
        def update_small_lines(mob):
            for i in range(len(small_loop_order)):
                num1 = small_loop_order[i]
                num2 = small_loop_order[(i + 1) % len(small_loop_order)]
                mob[i].put_start_and_end_on(
                    circles[num1].get_center(),
                    circles[num2].get_center()
                )
        
        small_loop_lines.add_updater(update_small_lines)
        
        # Create smooth animation for small loop
        total_time_small = 3  # Faster since fewer circles
        time_tracker_small = ValueTracker(0)
        
        # Store original positions
        original_positions_small = {num: circle_positions[num].copy() for num in small_loop_order}
        
        # Create list of position coordinates in order
        path_positions_small = [original_positions_small[num] for num in small_loop_order]
        
        # Track counted circles
        circles_counted_small = set()
        
        # Find which index in small_loop_order is position 62
        detector_pos_index = small_loop_order.index(62)
        
        # Highlight circle 59 throughout
        circles[59].set_stroke(width=4, color=YELLOW)
        
        def update_small_circles(mob, dt):
            t = time_tracker_small.get_value()
            total_positions = len(small_loop_order)
            
            for i, circle_num in enumerate(small_loop_order):
                # Calculate position along path
                progress = (t + i) % total_positions
                
                # Interpolate between positions
                pos_index = int(progress)
                next_pos_index = (pos_index + 1) % total_positions
                alpha = progress - pos_index
                
                current_pos = path_positions_small[pos_index]
                next_pos = path_positions_small[next_pos_index]
                interpolated_pos = current_pos * (1 - alpha) + next_pos * alpha
                
                circles[circle_num].move_to(interpolated_pos)
                
                # Check if circle is at detector position (62)
                if i == detector_pos_index and 0 <= (t % total_positions) < 1:
                    shifted_i = int(t) % total_positions
                    counting_circle = small_loop_order[shifted_i]
                    if counting_circle not in circles_counted_small:
                        circles_counted_small.add(counting_circle)
                        color = circle_colors[counting_circle]
                        color_counts_small[color] += 1
                        
                        # Change detector color to match last circle that passed
                        detector.set_color(color)
        
        # Add updaters
        for circle_num in small_loop_order:
            circles[circle_num].add_updater(update_small_circles)
        
        # Animate
        self.play(
            time_tracker_small.animate.set_value(len(small_loop_order)),
            run_time=total_time_small,
            rate_func=linear
        )
        
        # Remove updaters
        for circle_num in small_loop_order:
            circles[circle_num].remove_updater(update_small_circles)
        small_loop_lines.remove_updater(update_small_lines)
        
        # Reset circle 59 stroke
        circles[59].set_stroke(width=2)
        self.wait(1)
        
        # Calculate average color for small loop
        total_small = len(small_loop_order)
        avg_color_small = [
            color_counts_small[RED] / total_small,
            color_counts_small[GREEN] / total_small,
            color_counts_small[BLUE] / total_small
        ]
        
        # Create average color circle for small loop
        avg_circle_small = Circle(radius=0.6, fill_opacity=0.8, stroke_width=3, stroke_color=WHITE)
        avg_circle_small.set_color(rgb_to_color(avg_color_small))
        avg_circle_small.to_edge(RIGHT).shift(DOWN * 1.5)
        
        avg_label_small = Text("Average Color\n(Small Loop)", font_size=18)
        avg_label_small.next_to(avg_circle_small, DOWN, buff=0.3)
        
        self.play(Create(avg_circle_small), Write(avg_label_small), run_time=1.5)
        self.wait(2)
        
        # Final explanation
        explanation = VGroup(
            Text("Ergodic: Samples all colors uniformly", font_size=20, color=YELLOW),
            Text("→ Mixed average color", font_size=18),
            Text("", font_size=14),
            Text("Non-ergodic: Trapped in subset", font_size=20, color=ORANGE),
            Text("→ Only sees red", font_size=18)
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.2)
        explanation.to_edge(DOWN).shift(0.3*UP)
        
        self.play(Write(explanation), run_time=3)
        self.wait(4)