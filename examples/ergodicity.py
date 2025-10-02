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
    def construct(self):
        # Title
        title = Text("Ergodicity: The Conga Line Analogy", font_size=36)
        title.to_edge(UP)
        self.play(Write(title))
        self.wait(0.5)
        
        # Parameters
        grid_size = 8
        circle_radius = 0.15
        spacing = 0.8
        
        # Define colors
        colors_list = [RED, BLUE, GREEN]
        
        # Create 8x8 grid positions
        grid_positions = []
        for i in range(grid_size):
            for j in range(grid_size):
                x = (i - grid_size/2 + 0.5) * spacing
                y = (j - grid_size/2 + 0.5) * spacing
                grid_positions.append(np.array([x, y, 0]))
        
        # Assign random colors
        np.random.seed(123)
        circle_colors = [colors_list[np.random.randint(0, 3)] for _ in range(64)]
        
        # Create snake path for large loop (58 circles)
        # Snake goes left-to-right on even rows, right-to-left on odd rows
        # This creates a continuous path with only horizontal/vertical connections
        large_loop_path = []
        for row in range(grid_size):
            if row % 2 == 0:  # Even rows: left to right
                for col in range(grid_size):
                    idx = row * grid_size + col
                    large_loop_path.append(idx)
            else:  # Odd rows: right to left
                for col in range(grid_size - 1, -1, -1):
                    idx = row * grid_size + col
                    large_loop_path.append(idx)
        
        # Remove 6 consecutive positions for the small loop (create a gap)
        # Remove indices 9-14 (part of rows 1-2)
        small_loop_grid_indices = [9, 10, 17, 18, 25, 26]  # 2x3 block in grid
        large_loop_indices = [idx for idx in large_loop_path if idx not in small_loop_grid_indices]
        
        # Create the small loop as a 2x3 grid positioned outside to the right
        # Position it to the right of the main grid
        small_loop_positions = []
        small_loop_colors = []
        for i in range(2):  # 2 rows
            for j in range(3):  # 3 columns
                x = (grid_size/2 + 1 + j) * spacing
                y = (i - 0.5) * spacing
                small_loop_positions.append(np.array([x, y, 0]))
                small_loop_colors.append(colors_list[np.random.randint(0, 3)])
        
        # Create path for small loop (snake through the 2x3 grid)
        small_loop_path = [0, 1, 2, 5, 4, 3]  # Snake pattern in the 6 positions
        
        # Create circles for large loop
        circles = VGroup()
        for i, pos in enumerate(grid_positions):
            if i not in small_loop_grid_indices:  # Skip positions reserved for small loop gap
                circle = Circle(
                    radius=circle_radius,
                    color=circle_colors[i],
                    fill_opacity=0.8,
                    stroke_width=2,
                    stroke_color=WHITE
                )
                circle.move_to(pos)
                circles.add(circle)
            else:
                # Add placeholder (empty) for removed positions
                circle = Circle(radius=circle_radius, fill_opacity=0, stroke_opacity=0)
                circle.move_to(pos)
                circles.add(circle)
        
        # Create circles for small loop (positioned outside)
        small_circles = VGroup()
        for i, pos in enumerate(small_loop_positions):
            circle = Circle(
                radius=circle_radius,
                color=small_loop_colors[i],
                fill_opacity=0.8,
                stroke_width=2,
                stroke_color=WHITE
            )
            circle.move_to(pos)
            small_circles.add(circle)
        
        # Show all circles
        self.play(Create(circles), Create(small_circles), run_time=2)
        self.wait(0.5)
        
        # Create conga line connections for large loop
        large_loop_lines = VGroup()
        for i in range(len(large_loop_indices)):
            idx1 = large_loop_indices[i]
            idx2 = large_loop_indices[(i + 1) % len(large_loop_indices)]
            
            line = Line(
                circles[idx1].get_center(),
                circles[idx2].get_center(),
                color=YELLOW,
                stroke_width=3
            )
            large_loop_lines.add(line)
        
        # Create conga line connections for small loop (snake pattern)
        small_loop_lines = VGroup()
        for i in range(len(small_loop_path)):
            idx1 = small_loop_path[i]
            idx2 = small_loop_path[(i + 1) % len(small_loop_path)]
            
            line = Line(
                small_circles[idx1].get_center(),
                small_circles[idx2].get_center(),
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
        large_label = Text("Ergodic (58 circles)", font_size=20, color=YELLOW)
        large_label.to_edge(LEFT).shift(UP)
        
        small_label = Text("Non-ergodic (6 circles)", font_size=20, color=ORANGE)
        small_label.next_to(large_label, DOWN, aligned_edge=LEFT)
        
        self.play(
            Write(large_label),
            Write(small_label),
            run_time=1.5
        )
        self.wait(0.5)
        
        # Place observer at a fixed position to the left of the grid
        observer_position = np.array([-4, 0, 0])
        observer = Triangle(color=WHITE, fill_opacity=1)
        observer.scale(0.3)
        observer.move_to(observer_position)
        observer.rotate(PI/2)  # Point to the right
        
        observer_label = Text("Observer", font_size=18, color=WHITE)
        observer_label.next_to(observer, LEFT, buff=0.2)
        
        self.play(
            Create(observer),
            Write(observer_label),
            run_time=1
        )
        self.wait(0.5)
        
        # Create color counter display
        counter_box = Rectangle(width=3, height=2, color=WHITE, stroke_width=2)
        counter_box.to_edge(RIGHT).shift(UP)
        
        counter_title = Text("Color Average", font_size=20)
        counter_title.next_to(counter_box, UP, buff=0.2)
        
        red_count = Integer(0, color=RED, font_size=24)
        blue_count = Integer(0, color=BLUE, font_size=24)
        green_count = Integer(0, color=GREEN, font_size=24)
        
        red_label = Text("Red: ", font_size=20, color=RED)
        blue_label = Text("Blue: ", font_size=20, color=BLUE)
        green_label = Text("Green: ", font_size=20, color=GREEN)
        
        red_group = VGroup(red_label, red_count).arrange(RIGHT)
        blue_group = VGroup(blue_label, blue_count).arrange(RIGHT)
        green_group = VGroup(green_label, green_count).arrange(RIGHT)
        
        counter_display = VGroup(red_group, blue_group, green_group).arrange(DOWN, aligned_edge=LEFT, buff=0.2)
        counter_display.move_to(counter_box.get_center())
        
        self.play(
            Create(counter_box),
            Write(counter_title),
            Write(counter_display),
            run_time=1
        )
        self.wait(0.5)
        
        # Animate circles moving through the loop (conga line effect)
        # All circles shift to the next position in the path
        color_counts = {RED: 0, BLUE: 0, GREEN: 0}
        step_time = 0.08
        
        # Store original positions for the large loop
        large_positions = [circles[idx].get_center() for idx in large_loop_indices]
        
        # Add updater for lines to follow circles automatically
        def update_large_lines(mob):
            for i in range(len(large_loop_indices)):
                idx1 = large_loop_indices[i]
                idx2 = large_loop_indices[(i + 1) % len(large_loop_indices)]
                mob[i].put_start_and_end_on(
                    circles[idx1].get_center(),
                    circles[idx2].get_center()
                )
        
        large_loop_lines.add_updater(update_large_lines)
        
        # Only show enough steps to count all colors (not all 58)
        num_steps = len(large_loop_indices)
        
        for step in range(num_steps):
            # Count the color of the circle at observer position
            observer_circle_idx = large_loop_indices[step % len(large_loop_indices)]
            circle_color = circle_colors[observer_circle_idx]
            color_counts[circle_color] += 1
            
            # Create animation list for all circles moving simultaneously
            animations = []
            
            # Each circle moves to the next position in the loop
            for i in range(len(large_loop_indices)):
                current_idx = large_loop_indices[i]
                next_position_idx = (i + 1) % len(large_loop_indices)
                next_pos = large_positions[next_position_idx]
                
                if current_idx == observer_circle_idx:
                    # Highlight the circle passing the observer
                    animations.append(circles[current_idx].animate.move_to(next_pos).set_stroke(width=4, color=WHITE))
                else:
                    animations.append(circles[current_idx].animate.move_to(next_pos))
            
            # Execute the movement
            self.play(*animations, run_time=step_time)
            
            # Reset highlighted circle stroke
            if step < num_steps - 1:  # Don't reset on last step
                circles[observer_circle_idx].set_stroke(width=2)
            
            # Update counter display periodically (every 10 steps)
            if (step + 1) % 10 == 0 or step == num_steps - 1:
                red_count.set_value(color_counts[RED])
                blue_count.set_value(color_counts[BLUE])
                green_count.set_value(color_counts[GREEN])
        
        large_loop_lines.remove_updater(update_large_lines)
        
        self.wait(1)
        
        # Calculate percentages
        total = len(large_loop_indices)
        red_pct = color_counts[RED] / total * 100
        blue_pct = color_counts[BLUE] / total * 100
        green_pct = color_counts[GREEN] / total * 100
        
        result_large = VGroup(
            Text(f"Red: {red_pct:.1f}%", font_size=18, color=RED),
            Text(f"Blue: {blue_pct:.1f}%", font_size=18, color=BLUE),
            Text(f"Green: {green_pct:.1f}%", font_size=18, color=GREEN)
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.1)
        result_large.next_to(counter_box, DOWN, buff=0.3)
        
        self.play(Write(result_large), run_time=1.5)
        self.wait(1)
        
        # Now observer looks at small loop (on the right edge)
        # Move observer to right side to observe the small loop
        small_observer_position = np.array([4, 0, 0])
        
        self.play(
            FadeOut(result_large),
            red_count.animate.set_value(0),
            blue_count.animate.set_value(0),
            green_count.animate.set_value(0),
            observer.animate.move_to(small_observer_position).rotate(PI),  # Face left
            observer_label.animate.next_to(small_observer_position + RIGHT * 0.5, RIGHT, buff=0.2),
            run_time=1.5
        )
        self.wait(0.5)
        
        # Animate small loop circles moving through their conga line
        color_counts_small = {RED: 0, BLUE: 0, GREEN: 0}
        
        # Store original positions for the small loop
        small_positions = [small_circles[idx].get_center() for idx in small_loop_path]
        
        # Add updater for lines to follow circles automatically
        def update_small_lines(mob):
            for i in range(len(small_loop_path)):
                idx1 = small_loop_path[i]
                idx2 = small_loop_path[(i + 1) % len(small_loop_path)]
                mob[i].put_start_and_end_on(
                    small_circles[idx1].get_center(),
                    small_circles[idx2].get_center()
                )
        
        small_loop_lines.add_updater(update_small_lines)
        
        for step in range(len(small_loop_path)):
            # Count the color of the circle at observer position
            observer_circle_idx = small_loop_path[step % len(small_loop_path)]
            circle_color = small_loop_colors[observer_circle_idx]
            color_counts_small[circle_color] += 1
            
            # Create animation list for all circles moving simultaneously
            animations = []
            
            # Each circle moves to the next position in the loop
            for i in range(len(small_loop_path)):
                current_idx = small_loop_path[i]
                next_position_idx = (i + 1) % len(small_loop_path)
                next_pos = small_positions[next_position_idx]
                
                if current_idx == observer_circle_idx:
                    # Highlight the circle passing the observer
                    animations.append(small_circles[current_idx].animate.move_to(next_pos).set_stroke(width=4, color=WHITE))
                else:
                    animations.append(small_circles[current_idx].animate.move_to(next_pos))
            
            # Execute the movement
            self.play(*animations, run_time=0.2)
            
            # Reset highlighted circle stroke
            if step < len(small_loop_path) - 1:
                small_circles[observer_circle_idx].set_stroke(width=2)
        
        small_loop_lines.remove_updater(update_small_lines)
        
        # Update counters
        self.play(
            red_count.animate.set_value(color_counts_small[RED]),
            blue_count.animate.set_value(color_counts_small[BLUE]),
            green_count.animate.set_value(color_counts_small[GREEN]),
            run_time=0.8
        )
        self.wait(1)
        
        # Calculate percentages for small loop
        total_small = len(small_loop_path)
        red_pct_small = color_counts_small[RED] / total_small * 100
        blue_pct_small = color_counts_small[BLUE] / total_small * 100
        green_pct_small = color_counts_small[GREEN] / total_small * 100
        
        result_small = VGroup(
            Text(f"Red: {red_pct_small:.1f}%", font_size=18, color=RED),
            Text(f"Blue: {blue_pct_small:.1f}%", font_size=18, color=BLUE),
            Text(f"Green: {green_pct_small:.1f}%", font_size=18, color=GREEN)
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.1)
        result_small.next_to(counter_box, DOWN, buff=0.3)
        
        self.play(Write(result_small), run_time=1.5)
        self.wait(1)
        
        # Final explanation
        explanation = VGroup(
            Text("Ergodic: Time average = Ensemble average", font_size=22, color=YELLOW),
            Text("(Visits all accessible states uniformly)", font_size=18),
            Text("", font_size=14),
            Text("Non-ergodic: Trapped in subset", font_size=22, color=ORANGE),
            Text("(Cannot access full ensemble)", font_size=18)
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.2)
        explanation.to_edge(DOWN).shift(0.3*UP)
        
        self.play(Write(explanation), run_time=3)
        self.wait(4)