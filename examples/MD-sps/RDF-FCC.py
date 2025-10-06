""" Visualizing the Radial Distribution Function in a Crystalline material vs a liquid"""

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
        a = 1.5  # Lattice constant (scaled for visualization)
        atom_radius = 0.1
        
        # Color scheme for distance groups
        # These will be reused for both arrows and RDF bars
        neighbor_colors = [RED, GREEN, BLUE, YELLOW, PURPLE, ORANGE]
        
        # Define FCC unit cell atom positions
        fcc_unit_positions = [
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
        
        # First, show a single unit cell
        center_offset = apply_perspective(np.array([a/2, a/2, a/2]))
        
        unit_cell_atoms = VGroup()
        for pos_3d in fcc_unit_positions:
            pos_2d = apply_perspective(pos_3d) - center_offset
            atom = Circle(radius=atom_radius, color=GOLD, fill_opacity=1, stroke_width=2, stroke_color=YELLOW)
            atom.move_to(pos_2d)
            unit_cell_atoms.add(atom)
        
        unit_cube = create_cube_edges(np.array([0, 0, 0]), a)
        unit_cube.shift(-center_offset)
        
        # Show the single unit cell
        unit_cell_label = Text("Single FCC Unit Cell", font_size=24).next_to(unit_cube, DOWN, buff=0.5)
        self.play(Create(unit_cube), run_time=1)
        self.play(Create(unit_cell_atoms), Write(unit_cell_label), run_time=1.5)
        self.wait(1)
        
        # Now create a 2x2x2 lattice (smaller to fit better)
        # Generate all atom positions in the 2x2x2 supercell
        all_positions_3d = []
        for i in range(3):  # 3 cells to ensure complete 2x2x2 visible
            for j in range(3):
                for k in range(3):
                    for unit_pos in fcc_unit_positions:
                        pos = unit_pos + np.array([i*a, j*a, k*a])
                        all_positions_3d.append(pos)
        
        # Remove duplicates (atoms shared between unit cells)
        unique_positions_3d = []
        tolerance = 0.01
        for pos in all_positions_3d:
            is_duplicate = False
            for existing_pos in unique_positions_3d:
                if np.linalg.norm(pos - existing_pos) < tolerance:
                    is_duplicate = True
                    break
            if not is_duplicate:
                unique_positions_3d.append(pos)
        
        # Center the lattice around (1a, 1a, 1a) for 2x2x2
        lattice_center = np.array([1*a, 1*a, 1*a])
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
        
        # Create cube edges for 2x2x2 lattice
        large_cubes = VGroup()
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    cube = create_cube_edges(np.array([i*a, j*a, k*a]), a, WHITE, 0.2)
                    cube.shift(-large_center_offset)
                    large_cubes.add(cube)
        
        # Scale down to fit
        scale_factor = 0.8
        large_atoms.scale(scale_factor)
        large_cubes.scale(scale_factor)
        atom_positions_2d = [pos * scale_factor for pos in atom_positions_2d]
        
        # Shift to left side
        large_atoms.shift(3.5*LEFT)
        large_cubes.shift(3.5*LEFT)
        atom_positions_2d = [pos + 3.5*LEFT for pos in atom_positions_2d]
        
        # Transition from single unit cell to large lattice
        self.play(
            FadeOut(unit_cell_atoms),
            FadeOut(unit_cube),
            FadeOut(unit_cell_label),
            run_time=0.5
        )
        
        large_label = Text("2×2×2 FCC Lattice", font_size=24).to_edge(LEFT).shift(0.5*DOWN)
        self.play(
            Create(large_cubes),
            Create(large_atoms),
            Write(large_label),
            run_time=2
        )
        self.wait(1)
        
        # Create RDF axes on the right
        axes = Axes(
            x_range=[0, 2.5, 0.5],
            y_range=[0, 25, 5],
            x_length=5,
            y_length=4,
            axis_config={"include_tip": True, "include_numbers": True},
            tips=True,
        ).shift(3*RIGHT + 0.5*DOWN)
        
        x_label = axes.get_x_axis_label(Text("r/a (lattice constant)", font_size=18))
        y_label = axes.get_y_axis_label(Text("# of atoms", font_size=20), edge=LEFT, direction=LEFT)
        
        axes_title = Text("Radial Distribution", font_size=22).next_to(axes, UP, buff=0.3)
        
        self.play(
            Create(axes),
            Write(x_label),
            Write(y_label),
            Write(axes_title),
            run_time=1.5
        )
        self.wait(0.5)
        
        # Find reference atom near the center
        # Look for atom closest to center position
        center_pos = np.array([1*a, 1*a, 1*a])
        min_dist = float('inf')
        reference_idx = 0
        for i, pos_3d in enumerate(atom_positions_3d_centered):
            dist = np.linalg.norm(pos_3d - center_pos)
            if dist < min_dist:
                min_dist = dist
                reference_idx = i
        
        reference_atom = large_atoms[reference_idx]
        
        # Highlight reference atom
        reference_highlight = Circle(
            radius=atom_radius * scale_factor * 1.5,
            color=WHITE,
            stroke_width=5
        ).move_to(reference_atom.get_center())
        
        reference_label = Text("Reference\nAtom", font_size=18, color=WHITE).next_to(reference_atom, DOWN, buff=0.3)
        
        self.play(
            Create(reference_highlight),
            Write(reference_label),
            reference_atom.animate.set_color(WHITE),
            run_time=1
        )
        self.wait(0.5)
        
        # Calculate distances from reference atom
        ref_pos_3d = atom_positions_3d_centered[reference_idx]
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
        
        # Limit to first few neighbor shells
        max_shells = min(3, len(unique_distances))
        unique_distances = unique_distances[:max_shells]
        
        # Process each neighbor shell
        for shell_idx, dist_normalized in enumerate(unique_distances):
            color = neighbor_colors[shell_idx % len(neighbor_colors)]
            neighbor_indices = distance_groups[dist_normalized]
            count = len(neighbor_indices)
            
            # Fade out all non-relevant atoms (make them transparent)
            # Keep reference atom and current shell neighbors opaque
            fade_anims = []
            for i, atom in enumerate(large_atoms):
                if i != reference_idx and i not in neighbor_indices:
                    fade_anims.append(atom.animate.set_fill_opacity(0.15).set_stroke_opacity(0.15))
            
            # Make current shell neighbors opaque and colored
            highlight_anims = []
            for idx in neighbor_indices:
                highlight_anims.append(large_atoms[idx].animate.set_fill_opacity(1).set_stroke_opacity(1).set_color(color))
            
            # Create arrows to all neighbors in this shell
            arrows = VGroup()
            for idx in neighbor_indices:
                arrow = Arrow(
                    reference_atom.get_center(),
                    atom_positions_2d[idx],
                    buff=atom_radius * scale_factor,
                    color=color,
                    stroke_width=3,
                    max_tip_length_to_length_ratio=0.15
                )
                arrows.add(arrow)
            
            # Show all arrows and fade other atoms
            shell_label = Text(
                f"Shell {shell_idx + 1}: {count} neighbors at r = {dist_normalized:.2f}a",
                font_size=20,
                color=color
            ).to_edge(DOWN)
            
            self.play(
                *fade_anims,
                *highlight_anims,
                *[GrowArrow(arrow) for arrow in arrows],
                Write(shell_label),
                run_time=1.5
            )
            self.wait(0.5)
            
            # Create RDF bar with better spacing
            bar_width = 15  # Thinner bars for better spacing
            x_pos = axes.c2p(dist_normalized, 0)
            bar_top = axes.c2p(dist_normalized, count)
            
            bar = Line(
                x_pos,
                bar_top,
                color=color,
                stroke_width=bar_width
            )
            
            count_label = Text(str(count), font_size=18, color=color).next_to(bar_top, UP, buff=0.1)
            
            # Fade arrows sequentially while bar grows
            arrow_fade_anims = [FadeOut(arrow, run_time=0.1) for arrow in arrows]
            
            # Play fade outs with bar growing incrementally
            self.play(
                Succession(*arrow_fade_anims),
                Create(bar),
                run_time=2
            )
            self.play(Write(count_label), run_time=0.5)
            self.wait(0.5)
            
            # Fade out shell label
            self.play(FadeOut(shell_label), run_time=0.3)
        
        self.wait(1)
        
        # Add final explanation
        explanation = VGroup(
            Text("FCC Crystal RDF:", font_size=22, color=YELLOW),
            Text("• Discrete peaks", font_size=18),
            Text("• Exact distances", font_size=18),
            Text("• Long-range order", font_size=18),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.2)
        explanation.to_edge(DOWN).shift(0.5*UP)
        
        self.play(Write(explanation), run_time=1.5)
        self.wait(3)