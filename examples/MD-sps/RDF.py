""" Visualizing the Radial Distribution Function in a Crystalline material vs a liquid"""

from manim import *
import numpy as np


class RDF(Scene):
    def construct(self):
        # Title
        title = Text("Radial Distribution Function: FCC Gold", font_size=36)
        title.to_edge(UP)
        self.play(Write(title))
        self.wait(0.5)
        
        # FCC lattice parameters
        a = 2.0  # Lattice constant (scaled for visualization)
        atom_radius = 0.12
        
        # Define FCC unit cell atom positions
        # Corner atoms (8 corners, but only 1/8 of each belongs to this cell)
        # Face-centered atoms (6 faces, 1/2 of each belongs to this cell)
        fcc_positions = [
            # Corner atoms
            np.array([0, 0, 0]),
            np.array([a, 0, 0]),
            np.array([a, a, 0]),
            np.array([0, a, 0]),
            np.array([0, 0, a]),
            np.array([a, 0, a]),
            np.array([a, a, a]),
            np.array([0, a, a]),
            # Face-centered atoms
            np.array([a/2, a/2, 0]),  # Front face
            np.array([a/2, a/2, a]),  # Back face
            np.array([a/2, 0, a/2]),  # Bottom face
            np.array([a/2, a, a/2]),  # Top face
            np.array([0, a/2, a/2]),  # Left face
            np.array([a, a/2, a/2]),  # Right face
        ]
        
        # Apply isometric projection for 3D appearance
        def apply_perspective(pos):
            x, y, z = pos
            iso_x = x - z * 0.6
            iso_y = (x + z) * 0.3 + y * 0.8
            return np.array([iso_x, iso_y, 0])
        
        # Create unit cell cube edges
        def create_cube_edges(size, color=WHITE, opacity=0.3):
            edges = VGroup()
            vertices_3d = [
                np.array([0, 0, 0]),
                np.array([size, 0, 0]),
                np.array([size, size, 0]),
                np.array([0, size, 0]),
                np.array([0, 0, size]),
                np.array([size, 0, size]),
                np.array([size, size, size]),
                np.array([0, size, size]),
            ]
            
            # Define edges by vertex indices
            edge_pairs = [
                (0, 1), (1, 2), (2, 3), (3, 0),  # Bottom face
                (4, 5), (5, 6), (6, 7), (7, 4),  # Top face
                (0, 4), (1, 5), (2, 6), (3, 7),  # Vertical edges
            ]
            
            for i, j in edge_pairs:
                v1 = apply_perspective(vertices_3d[i])
                v2 = apply_perspective(vertices_3d[j])
                edge = Line(v1, v2, color=color, stroke_width=2, stroke_opacity=opacity)
                edges.add(edge)
            
            return edges
        
        # Center the visualization
        center_offset = apply_perspective(np.array([a/2, a/2, a/2]))
        
        # Create atoms
        atoms = VGroup()
        atom_positions_2d = []
        
        for pos_3d in fcc_positions:
            pos_2d = apply_perspective(pos_3d) - center_offset
            atom = Circle(radius=atom_radius, color=GOLD, fill_opacity=1, stroke_width=2, stroke_color=YELLOW)
            atom.move_to(pos_2d)
            atoms.add(atom)
            atom_positions_2d.append(pos_2d)
        
        # Create cube
        cube = create_cube_edges(a)
        cube.shift(-center_offset)
        
        # Shift everything to the left
        atoms.shift(2.5*LEFT)
        cube.shift(2.5*LEFT)
        atom_positions_2d = [pos + 2.5*LEFT for pos in atom_positions_2d]
        
        # Show the FCC structure
        self.play(Create(cube), run_time=1)
        self.play(Create(atoms), run_time=1.5)
        self.wait(0.5)
        
        # Choose reference atom (corner atom at origin)
        reference_idx = 0
        reference_atom = atoms[reference_idx]
        
        # Highlight reference atom
        reference_highlight = Circle(
            radius=atom_radius * 1.5,
            color=RED,
            stroke_width=4
        ).move_to(reference_atom.get_center())
        
        reference_label = Text("Reference", font_size=20, color=RED).next_to(reference_atom, DOWN, buff=0.3)
        
        self.play(
            Create(reference_highlight),
            Write(reference_label),
            reference_atom.animate.set_color(RED),
            run_time=1
        )
        self.wait(0.5)
        
        # Calculate distances from reference atom
        ref_pos = atom_positions_2d[reference_idx]
        distances_3d = []
        distance_vectors = VGroup()
        
        for i, pos_3d in enumerate(fcc_positions):
            if i != reference_idx:
                # Calculate actual 3D distance
                dist_3d = np.linalg.norm(fcc_positions[i] - fcc_positions[reference_idx])
                distances_3d.append(dist_3d)
                
                # Create visual vector
                vector = Arrow(
                    ref_pos,
                    atom_positions_2d[i],
                    buff=atom_radius,
                    color=BLUE,
                    stroke_width=3,
                    max_tip_length_to_length_ratio=0.15
                )
                distance_vectors.add(vector)
        
        # Show distance vectors
        self.play(
            *[GrowArrow(vec) for vec in distance_vectors],
            run_time=2
        )
        self.wait(1)
        
        # Group distances by length (round to handle floating point)
        distance_groups = {}
        for i, dist in enumerate(distances_3d):
            dist_rounded = round(dist, 2)
            if dist_rounded not in distance_groups:
                distance_groups[dist_rounded] = []
            distance_groups[dist_rounded].append(i)
        
        # Sort distances
        unique_distances = sorted(distance_groups.keys())
        
        # Color code vectors by distance group
        colors = [GREEN, YELLOW, ORANGE, PURPLE, PINK, TEAL]
        
        # Animate grouping by color
        color_anims = []
        for group_idx, dist in enumerate(unique_distances):
            color = colors[group_idx % len(colors)]
            for vec_idx in distance_groups[dist]:
                color_anims.append(distance_vectors[vec_idx].animate.set_color(color))
        
        self.play(*color_anims, run_time=1.5)
        self.wait(1)
        
        # Fade out atoms and cube, keep vectors
        self.play(
            FadeOut(atoms),
            FadeOut(cube),
            FadeOut(reference_highlight),
            FadeOut(reference_label),
            run_time=1
        )
        self.wait(0.5)
        
        # Move vectors to the left and group them
        grouped_vectors = VGroup()
        y_offset = 2.5
        
        for group_idx, dist in enumerate(unique_distances):
            group_vecs = VGroup()
            count = len(distance_groups[dist])
            color = colors[group_idx % len(colors)]
            
            for i, vec_idx in enumerate(distance_groups[dist]):
                vec = distance_vectors[vec_idx]
                group_vecs.add(vec)
            
            grouped_vectors.add(group_vecs)
        
        # Rearrange vectors vertically on the left
        # Sort the vectors from shortest to longest, left to right, where the fcc lattice was. make them straight up and down.
        self.play(
            distance_vectors.animate.arrange(DOWN, center=True, buff=0.2).shift(3.5*LEFT + 0.5*UP),
            run_time=1.5
        )
        self.wait(0.5)
        
        # Create RDF plot on the right
        axes = Axes(
            x_range=[0, 3.5, 0.5],
            y_range=[0, 8, 2],
            x_length=5,
            y_length=4,
            axis_config={"include_tip": True},
            tips=True,
        ).shift(2*RIGHT + 0.3*DOWN)
        
        x_label = axes.get_x_axis_label(Text("Distance r", font_size=24))
        y_label = axes.get_y_axis_label(Text("g(r)", font_size=24), edge=LEFT, direction=LEFT)
        
        self.play(
            Create(axes),
            Write(x_label),
            Write(y_label),
            run_time=1.5
        )
        self.wait(0.5)
        
        # Build RDF bars
        bars = VGroup()
        bar_labels = VGroup()
        
        for group_idx, dist in enumerate(unique_distances):
            count = len(distance_groups[dist])
            color = colors[group_idx % len(colors)]
            
            # Create bar
            x_pos = axes.c2p(dist, 0)
            bar_height = count
            bar_top = axes.c2p(dist, bar_height)
            
            bar = Line(
                x_pos,
                bar_top,
                color=color,
                stroke_width=20
            )
            bars.add(bar)
            
            # Add count label
            count_label = Text(str(count), font_size=18, color=color).next_to(bar_top, UP, buff=0.1)
            bar_labels.add(count_label)
        
        # Animate bars growing from vectors
        bar_anims = []
        for i, bar in enumerate(bars):
            bar_anims.append(Create(bar))
        
        self.play(*bar_anims, run_time=2)
        self.play(*[Write(label) for label in bar_labels], run_time=1)
        self.wait(1)
        
        # Add distance labels on x-axis
        distance_labels = VGroup()
        for dist in unique_distances:
            pos = axes.c2p(dist, 0)
            label = Text(f"{dist:.2f}", font_size=16).next_to(pos, DOWN, buff=0.2)
            distance_labels.add(label)
        
        self.play(Write(distance_labels), run_time=1)
        
        # Add explanation
        explanation = VGroup(
            Text("Radial Distribution Function:", font_size=22, color=YELLOW),
            Text("Counts atoms at each distance", font_size=18),
            Text("from the reference atom", font_size=18),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.2)
        explanation.to_edge(DOWN).shift(0.3*UP)
        
        self.play(Write(explanation), run_time=1.5)
        self.wait(3)
class LiquidRDF(Scene):
    def construct(self):
        # Title
        title = Text("Radial Distribution Function: Liquid Gold", font_size=36)
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
            x_range=[0, 6, 1],
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
        
        # Calculate RDF histogram
        max_dist = 5.5
        n_bins = 50
        bin_edges = np.linspace(0, max_dist, n_bins + 1)
        hist, _ = np.histogram(distances, bins=bin_edges)
        
        # Normalize by shell area (2πr * dr) and density
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        dr = bin_edges[1] - bin_edges[0]
        
        # Calculate g(r) - normalize by ideal gas (random distribution)
        area = box_size * box_size
        number_density = len(liquid_positions) / area
        
        g_r = []
        for i, r in enumerate(bin_centers):
            if r > 0:
                shell_area = 2 * np.pi * r * dr
                expected_count = number_density * shell_area
                if expected_count > 0:
                    g_r.append(hist[i] / expected_count)
                else:
                    g_r.append(0)
            else:
                g_r.append(0)
        
        # Smooth the curve
        from scipy.ndimage import gaussian_filter1d
        g_r_smooth = gaussian_filter1d(g_r, sigma=1.5)
        
        # Create smooth curve
        points = []
        for i, (r, g) in enumerate(zip(bin_centers, g_r_smooth)):
            if r <= max_dist:
                g_clamped = max(0, min(g, 3))  # Clamp to axis range
                points.append(axes.c2p(r, g_clamped))
        
        # Draw the RDF curve
        rdf_curve = VMobject(color=BLUE, stroke_width=4)
        rdf_curve.set_points_as_corners(points)
        
        self.play(Create(rdf_curve), run_time=3, rate_func=linear)
        self.wait(1)
        
        # Add annotations for liquid characteristics
        annotation1 = Text("No atoms closer than", font_size=20, color=YELLOW)
        annotation2 = Text("minimum distance", font_size=20, color=YELLOW)
        annotation_group1 = VGroup(annotation1, annotation2).arrange(DOWN, buff=0.1)
        annotation_group1.next_to(axes.c2p(0.3, 0.2), UP+RIGHT, buff=0.3)
        
        arrow1 = Arrow(
            annotation_group1.get_bottom(),
            axes.c2p(0.3, 0),
            color=YELLOW,
            stroke_width=2,
            max_tip_length_to_length_ratio=0.2
        )
        
        self.play(
            Write(annotation_group1),
            GrowArrow(arrow1),
            run_time=1.5
        )
        self.wait(1)
        
        # Find first peak
        first_peak_idx = np.argmax(g_r_smooth[:15])
        peak_r = bin_centers[first_peak_idx]
        peak_g = g_r_smooth[first_peak_idx]
        
        peak_annotation = Text("First coordination shell", font_size=20, color=GREEN)
        peak_annotation.next_to(axes.c2p(peak_r, peak_g), UP, buff=0.3)
        
        peak_dot = Dot(axes.c2p(peak_r, peak_g), color=GREEN, radius=0.08)
        
        self.play(
            Create(peak_dot),
            Write(peak_annotation),
            run_time=1.5
        )
        self.wait(1)
        
        # Add final explanation
        explanation = VGroup(
            Text("Liquid RDF:", font_size=24, color=YELLOW),
            Text("• Smooth, continuous curve", font_size=20),
            Text("• Oscillates around g(r) = 1", font_size=20),
            Text("• Short-range order only", font_size=20),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.2)
        explanation.to_edge(RIGHT).shift(0.5*UP)
        
        self.play(Write(explanation), run_time=2)
        self.wait(3)
