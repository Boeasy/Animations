""" Visualizing the Radial Distribution Function in a FCC Crystal"""

from manim import *
import numpy as np


class RDF_FCC(Scene):
    def construct(self):
        # Title
        title = Text("Radial Distribution Function: FCC Gold Crystal", font_size=36)
        title.to_edge(UP)
        self.play(Write(title))
        self.wait(0.5)
        
        # FCC lattice parameters
        a = 1.2  # Lattice constant (scaled for visualization)
        atom_radius = 0.08
        
        # Define conventional FCC unit cell (14 atoms - corners + faces for visualization)
        fcc_conventional_positions = [
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
        
        # Define standard 4-atom FCC primitive unit cell
        # This is the minimal repeating unit
        fcc_primitive_positions = [
            np.array([0, 0, 0]),      # Corner
            np.array([a/2, a/2, 0]),  # Front face center
            np.array([a/2, 0, a/2]),  # Bottom face center
            np.array([0, a/2, a/2]),  # Left face center
        ]
        
        # Apply isometric projection for 3D appearance
        def apply_perspective(pos):
            x, y, z = pos
            iso_x = x - z * 0.6
            iso_y = (x + z) * 0.3 + y * 0.8
            return np.array([iso_x, iso_y, 0])
        
        # Create unit cell cube edges
        def create_cube_edges(offset_3d, size, color=WHITE, opacity=0.3):
            edges = VGroup()
            vertices_3d = [
                offset_3d + np.array([0, 0, 0]),
                offset_3d + np.array([size, 0, 0]),
                offset_3d + np.array([size, size, 0]),
                offset_3d + np.array([0, size, 0]),
                offset_3d + np.array([0, 0, size]),
                offset_3d + np.array([size, 0, size]),
                offset_3d + np.array([size, size, size]),
                offset_3d + np.array([0, size, size]),
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
        
        # First, show a single unit cell (using conventional 14-atom representation)
        center_offset = apply_perspective(np.array([a/2, a/2, a/2]))
        
        unit_cell_atoms = VGroup()
        for pos_3d in fcc_conventional_positions:
            pos_2d = apply_perspective(pos_3d) - center_offset
            atom = Circle(radius=atom_radius, color=GOLD, fill_opacity=1, stroke_width=2, stroke_color=YELLOW)
            atom.move_to(pos_2d)
            unit_cell_atoms.add(atom)
        
        unit_cube = create_cube_edges(np.array([0, 0, 0]), a)
        unit_cube.shift(-center_offset)
        
        # Show the single unit cell
        unit_cell_label = Text("FCC Crystal Structure", font_size=24).next_to(unit_cube, DOWN, buff=0.5)
        self.play(Create(unit_cube), run_time=1)
        self.play(Create(unit_cell_atoms), Write(unit_cell_label), run_time=1.5)
        self.wait(1)
        
        # Now create a 4x4x4 supercell using the standard 4-atom primitive cell
        # This avoids duplicates by using the primitive unit cell
        unique_positions_3d = []
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for unit_pos in fcc_primitive_positions:
                        pos = unit_pos + np.array([i*a, j*a, k*a])
                        unique_positions_3d.append(pos)
        
        # Center the lattice around (2a, 2a, 2a) for 4x4x4
        lattice_center = np.array([2*a, 2*a, 2*a])
        large_center_offset = apply_perspective(lattice_center)
        
        # Create all atoms for large lattice
        large_atoms = VGroup()
        atom_positions_2d = []
        atom_positions_3d_centered = []
        
        for pos_3d in unique_positions_3d:
            pos_2d = apply_perspective(pos_3d) - large_center_offset
            atom = Circle(radius=atom_radius, color=GOLD, fill_opacity=1, stroke_width=2, stroke_color=YELLOW)
            atom.move_to(pos_2d)
            large_atoms.add(atom)
            atom_positions_2d.append(pos_2d)
            atom_positions_3d_centered.append(pos_3d)
        
        # Create cube edges for 4x4x4 lattice
        large_cubes = VGroup()
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    cube = create_cube_edges(np.array([i*a, j*a, k*a]), a, WHITE, 0.2)
                    cube.shift(-large_center_offset)
                    large_cubes.add(cube)
        
        # Scale down to fit on screen (larger than before)
        scale_factor = 1.0
        large_atoms.scale(scale_factor)
        large_cubes.scale(scale_factor)
        atom_positions_2d = [pos * scale_factor for pos in atom_positions_2d]
        
        # Center on screen
        large_atoms.shift(0.3*DOWN)
        large_cubes.shift(0.3*DOWN)
        atom_positions_2d = [pos + 0.3*DOWN for pos in atom_positions_2d]
        
        # Transition from single unit cell to large lattice
        self.play(
            FadeOut(unit_cell_atoms),
            FadeOut(unit_cube),
            FadeOut(unit_cell_label),
            run_time=0.5
        )
        
        large_label = Text("4×4×4 FCC Supercell", font_size=24).next_to(large_cubes, DOWN, buff=0.8)
        self.play(
            Create(large_cubes),
            Create(large_atoms),
            Write(large_label),
            run_time=2
        )
        self.wait(1)
        
        # Find reference atom near the center
        center_pos = np.array([2*a, 2*a, 2*a])
        min_dist = float('inf')
        reference_idx = 0
        for i, pos_3d in enumerate(atom_positions_3d_centered):
            dist = np.linalg.norm(pos_3d - center_pos)
            if dist < min_dist:
                min_dist = dist
                reference_idx = i
        
        reference_atom = large_atoms[reference_idx]
        ref_pos_3d = atom_positions_3d_centered[reference_idx]

        # Make all atoms semi-transparent
        self.play(
            *[atom.animate.set_fill(opacity=0.2).set_stroke(opacity=0.2) for atom in large_atoms],
            FadeOut(large_label),
            run_time=1
        )
        
        # Make reference atom opaque and highlighted
        reference_highlight = Circle(
            radius=atom_radius * scale_factor * 1.8,
            color=WHITE,
            stroke_width=5,
            fill_opacity=0
        ).move_to(reference_atom.get_center())
        
        self.play(
            reference_atom.animate.set_fill(opacity=1).set_stroke(opacity=1).set_color(GOLD),
            Create(reference_highlight),
            run_time=1
        )
        self.wait(0.5)
        
        # Calculate distances from reference atom
        distances_3d = []
        atom_indices = []
        
        for i, pos_3d in enumerate(atom_positions_3d_centered):
            if i != reference_idx:
                dist_3d = np.linalg.norm(pos_3d - ref_pos_3d)
                distances_3d.append(dist_3d)
                atom_indices.append(i)
        
        # Group atoms by distance (nearest neighbors, next-nearest, etc.)
        distance_groups = {}
        for i, dist in enumerate(distances_3d):
            dist_normalized = dist / a  # Normalize by lattice constant
            dist_rounded = round(dist_normalized, 2)
            if dist_rounded not in distance_groups:
                distance_groups[dist_rounded] = []
            distance_groups[dist_rounded].append(atom_indices[i])
        
        # Sort distances
        unique_distances = sorted(distance_groups.keys())
        
        # Define the three neighbor shells we want with colors
        # 1st nearest: 12 atoms at distance a/sqrt(2) ≈ 0.707a (face diagonal)
        # 2nd nearest: 6 atoms at distance a (edge)
        # 3rd nearest: 24 atoms at distance a*sqrt(3/2) ≈ 1.225a
        neighbor_shells = [
            {"color": RED, "label": "1st Nearest Neighbors", "expected_count": 12},
            {"color": GREEN, "label": "2nd Nearest Neighbors", "expected_count": 6},
            {"color": BLUE, "label": "3rd Nearest Neighbors", "expected_count": 24},
        ]
        
        # Take the first 3 distance groups
        shell_distances = unique_distances[:3]
        
        # Store count labels to persist
        count_labels = []
        
        # Process each neighbor shell
        for shell_idx, dist_normalized in enumerate(shell_distances):
            shell_info = neighbor_shells[shell_idx]
            color = shell_info["color"]
            label_text = shell_info["label"]
            neighbor_indices = distance_groups[dist_normalized]
            count = len(neighbor_indices)
            
            # Make current shell neighbors opaque
            highlight_anims = []
            for idx in neighbor_indices:
                highlight_anims.append(
                    large_atoms[idx].animate.set_fill(opacity=1).set_stroke(opacity=1).set_color(GOLD)
                )
            
            # Create arrows to all neighbors in this shell
            arrows = VGroup()
            for idx in neighbor_indices:
                arrow = Arrow(
                    reference_atom.get_center(),
                    atom_positions_2d[idx],
                    buff=atom_radius * scale_factor * 1.2,
                    color=color,
                    stroke_width=4,
                    max_tip_length_to_length_ratio=0.12
                )
                arrows.add(arrow)
            
            self.play(
                *highlight_anims,
                *[GrowArrow(arrow) for arrow in arrows],
                run_time=1.5
            )
            self.wait(0.5)
            
            # Create count label at bottom
            count_label = Text(
                f"{label_text}: {count}",
                font_size=24,
                color=color
            )
            
            # Position labels in a row at the bottom
            if shell_idx == 0:
                count_label.to_edge(LEFT)
            elif shell_idx == 1:
                count_label.to_edge(LEFT).shift(1.5*DOWN)
            else:
                count_label.to_edge(LEFT).shift(3.0*DOWN)
            
            count_labels.append(count_label)
            
            self.play(Write(count_label), run_time=0.5)
            self.wait(0.5)
            
            # Fade arrows and return atoms to semi-transparent
            fade_anims = []
            for idx in neighbor_indices:
                fade_anims.append(
                    large_atoms[idx].animate.set_fill(opacity=0.2).set_stroke(opacity=0.2).set_color(GOLD)
                )
            
            self.play(
                *[FadeOut(arrow) for arrow in arrows],
                *fade_anims,
                run_time=1
            )
            self.wait(0.3)
        
        # Keep count labels, fade out the lattice
        self.play(
            FadeOut(large_atoms),
            FadeOut(large_cubes),
            FadeOut(reference_highlight),
            *[FadeOut(lbl) for lbl in count_labels],
            run_time=1
        )
        self.wait(0.5)
        
        # Create RDF graph
        axes = Axes(
            x_range=[0, 2.0, 0.5],
            y_range=[0, 30, 5],
            x_length=10,
            y_length=5,
            axis_config={"include_tip": True, "include_numbers": True},
            tips=True,
        ).shift(0.5*UP)
        
        x_label = axes.get_x_axis_label(
            Text("Distance", font_size=20),
            edge=DOWN,
            direction=DOWN
            )
        y_label = axes.get_y_axis_label(
            Text("Number of atoms", font_size=20),
            edge=LEFT,
            direction=UP)
        y_label.rotate(PI/2)
        y_label.next_to(axes.y_axis, LEFT, buff=0.2)
        
        #axes_title = Text("Radial Distribution Function", font_size=28).next_to(axes, UP, buff=0.4)
        
        self.play(
            Create(axes),
            Write(x_label),
            Write(y_label),
            #Write(axes_title),
            run_time=1.5
        )
        self.wait(0.5)
        
        # Create RDF peaks corresponding to the three neighbor shells
        peaks = []
        for shell_idx, dist_normalized in enumerate(shell_distances):
            shell_info = neighbor_shells[shell_idx]
            color = shell_info["color"]
            count = len(distance_groups[dist_normalized])
            
            # Create peak as a thin vertical line
            bar_width = 8
            x_pos = axes.c2p(dist_normalized, 0)
            bar_top = axes.c2p(dist_normalized, count)
            
            bar = Line(
                x_pos,
                bar_top,
                color=color,
                stroke_width=bar_width
            )
            
            count_text = Text(str(count), font_size=20, color=color).next_to(bar_top, UP, buff=0.15)
            
            peaks.append((bar, count_text))
            
            self.play(
                Create(bar),
                Write(count_text),
                run_time=1
            )
            self.wait(0.3)
        
        # Create smooth RDF curve passing through peak tops
        # create a curve with Gaussian-like peaks at each distance
        height_factors = [1.05, 1.0, 1.17] #fudge factor to make graph look nice :D
        def rdf_curve(r):
            """Generate realistic RDF curve with peaks at shell distances"""
            g_r = 0
            # Add Gaussian peaks at each shell distance
            for shell_idx, dist in enumerate(shell_distances):
                count = len(distance_groups[dist])
                count *= height_factors[shell_idx]
                # Peak height proportional to count, with some width
                peak_width = 0.03  # Width of each peak
                g_r += count * np.exp(-((r - dist) ** 2) / (2 * peak_width ** 2))
            return g_r
        
        # Create the RDF curve
        rdf_graph = axes.plot(
            rdf_curve,
            x_range=[0.4, 1.6],  # Range covering our three peaks
            color=YELLOW,
            stroke_width=3,
        )
        
        self.play(Create(rdf_graph), run_time=2)
        self.wait(2)