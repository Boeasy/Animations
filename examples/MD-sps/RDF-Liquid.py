""" Visualizing the Radial Distribution Function in a Crystalline material vs a liquid"""

from manim import *
import numpy as np

class LiquidRDF(Scene):
    def construct(self):
        # Title
        title = Text("Radial Distribution Function: Liquid", font_size=36)
        title.to_edge(UP)
        self.play(Write(title))
        self.wait(0.5)
        
        # Parameters
        n_atoms = 80  # More atoms for liquid
        atom_radius = 0.08
        box_size = 6.0
        min_distance = 0.35  # Minimum separation between atoms
        
        # Generate random liquid structure (2D projection)
        np.random.seed(42)  # For reproducibility
        liquid_positions = []
        
        # Generate positions with minimum distance constraint
        max_attempts = 1000
        while len(liquid_positions) < n_atoms and max_attempts > 0:
            pos = np.random.uniform(-box_size/2, box_size/2, 2)
            
            # Check if too close to existing atoms
            too_close = False
            for existing_pos in liquid_positions:
                if np.linalg.norm(pos - existing_pos) < min_distance:
                    too_close = True
                    break
            
            if not too_close:
                liquid_positions.append(pos)
            max_attempts -= 1
        
        # Create atoms
        atoms = VGroup()
        for pos_2d in liquid_positions:
            atom = Circle(radius=atom_radius, color=GOLD, fill_opacity=0.9, stroke_width=1, stroke_color=YELLOW)
            atom.move_to([pos_2d[0], pos_2d[1], 0])
            atoms.add(atom)
        
        # Shift to left
        atoms.shift(3*LEFT)
        liquid_positions = [pos + np.array([-3, 0]) for pos in liquid_positions]
        
        # Show the liquid structure
        self.play(Create(atoms), run_time=2)
        self.wait(0.5)
        
        # Choose reference atom (roughly in the middle)
        reference_idx = len(liquid_positions) // 2
        reference_atom = atoms[reference_idx]
        
        # Highlight reference atom
        reference_highlight = Circle(
            radius=atom_radius * 1.8,
            color=RED,
            stroke_width=4
        ).move_to(reference_atom.get_center())
        
        reference_label = Text("Reference", font_size=18, color=RED).next_to(reference_atom, DOWN, buff=0.2)
        
        self.play(
            Create(reference_highlight),
            Write(reference_label),
            reference_atom.animate.set_color(RED),
            run_time=1
        )
        self.wait(0.5)
        
        # Calculate distances from reference atom
        ref_pos = liquid_positions[reference_idx]
        distances = []
        distance_vectors = VGroup()
        
        for i, pos_2d in enumerate(liquid_positions):
            if i != reference_idx:
                # Calculate 2D distance
                dist = np.linalg.norm(pos_2d - ref_pos)
                distances.append(dist)
                
                # Create visual vector
                vector = Arrow(
                    [ref_pos[0], ref_pos[1], 0],
                    [pos_2d[0], pos_2d[1], 0],
                    buff=atom_radius,
                    color=BLUE,
                    stroke_width=2,
                    max_tip_length_to_length_ratio=0.1
                )
                distance_vectors.add(vector)
        
        # Show distance vectors
        self.play(
            *[GrowArrow(vec) for vec in distance_vectors],
            run_time=2
        )
        self.wait(1)
        
        # Fade out atoms but keep vectors briefly
        self.play(
            FadeOut(atoms),
            FadeOut(reference_highlight),
            FadeOut(reference_label),
            run_time=1
        )
        self.wait(0.3)
        
        # Fade out vectors and transition to graph
        self.play(
            FadeOut(distance_vectors),
            run_time=0.8
        )
        self.wait(0.3)
        
        # Create RDF plot
        axes = Axes(
            x_range=[0, 2, 0.5],
            y_range=[0, 3, 0.5],
            x_length=10,
            y_length=5,
            axis_config={"include_tip": True},
            tips=True,
        ).shift(0.3*DOWN)
        
        x_label = axes.get_x_axis_label(Text("Distance r", font_size=28))
        y_label = axes.get_y_axis_label(Text("g(r)", font_size=28), edge=LEFT, direction=LEFT)
        
        self.play(
            Create(axes),
            Write(x_label),
            Write(y_label),
            run_time=1.5
        )
        self.wait(0.5)
        
        import math
        sqrt2 = math.sqrt(2)
        
        def rdf_curve(r):
            """Generate liquid RDF curve with sharp first peak and rounded second peak"""
            g_r = 0
            
            # Stay near zero until first neighbor distance (around √2)
            if r < .4:
                return 0.05 * r  # Very slight rise from zero
            
            # First peak: SHARP and STEEP at first coordination shell
            first_peak_pos = 0.5  # First shell around √2
            first_peak_height = 2.0
            first_peak_width = 0.2  # Narrow peak = steep
            g_r += first_peak_height * np.exp(-((r - first_peak_pos) ** 2) / (2 * first_peak_width ** 2))
            
            # Second peak: BROAD and ROUNDED at second shell
            second_peak_pos = 1.3  # Second coordination shell
            second_peak_height = 1.4  # Lower than first
            second_peak_width = 0.5  # Much wider = more rounded
            g_r += second_peak_height * np.exp(-((r - second_peak_pos) ** 2) / (2 * second_peak_width ** 2))
            
            # Add asymptotic approach to g(r) = 1 at large distances
            g_r += 1.0 * (1 - np.exp(-r / 3))
            
            return g_r
        
        # Create the RDF curve
        rdf_curve_plot = axes.plot(
            rdf_curve,
            x_range=[0.0, 2],
            color=BLUE,
            stroke_width=4,
        )
        
        self.play(Create(rdf_curve_plot), run_time=3, rate_func=linear)
        self.wait(1)
        
        self.wait(3)
