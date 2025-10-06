from manim import *
import numpy as np


class MLPotential(Scene):
    """Demonstrate the difference between EAM and ML potentials."""
    #This is basically where it needs to be, but let's try exaggerating the shapes and movement for easy viewing
    #and then add the equations (general form) along with it.
    
    def construct(self):
        # Title
        title = Text("ML Potential vs. EAM", font_size=40)
        title.to_edge(UP)
        self.play(FadeIn(title))
        
        # Parameters
        atom_radius = 0.2
        initial_distance = 3.5
        gold_color = GOLD
        cloud_color = YELLOW
        pair_color = GREEN
        repulsion_color = RED
        
        # ============================================================
        # PART 1: Traditional EAM with symmetric electron clouds
        # ============================================================
        
        eam_label = Text("Traditional EAM", font_size=30, color=BLUE)
        eam_label.next_to(title, DOWN, buff=0.3)
        self.play(Write(eam_label))
        
        # Create two atoms
        atom1 = Circle(radius=atom_radius, stroke_width=3)
        atom1.set_fill(gold_color, opacity=1.0)
        atom1.set_stroke(GOLD_E, width=3)
        atom1.move_to(LEFT * initial_distance / 2)
        
        atom2 = Circle(radius=atom_radius, stroke_width=3)
        atom2.set_fill(gold_color, opacity=1.0)
        atom2.set_stroke(GOLD_E, width=3)
        atom2.move_to(RIGHT * initial_distance / 2)
        
        self.play(FadeIn(atom1), FadeIn(atom2))
        self.wait(0.5)
        
        # Show pair potential arrows (attractive)
        pair_label = Text("1. Pair Potential (Attractive)", font_size=24, color=pair_color)
        pair_label.to_corner(UL, buff=0.5).shift(DOWN * 2)
        
        arrow1_to_2 = Arrow(
            start=atom1.get_right() + RIGHT * 0.1,
            end=atom2.get_left() - RIGHT * 0.1,
            color=pair_color,
            stroke_width=4,
            max_tip_length_to_length_ratio=0.15
        )
        
        arrow2_to_1 = Arrow(
            start=atom2.get_left() - RIGHT * 0.1,
            end=atom1.get_right() + RIGHT * 0.1,
            color=pair_color,
            stroke_width=4,
            max_tip_length_to_length_ratio=0.15
        )
        
        self.play(Write(pair_label))
        self.play(GrowArrow(arrow1_to_2), GrowArrow(arrow2_to_1))
        self.wait(1)
        
        # Atoms move closer due to attraction
        new_distance = 2.8
        self.play(
            atom1.animate.move_to(LEFT * new_distance / 2),
            atom2.animate.move_to(RIGHT * new_distance / 2),
            FadeOut(arrow1_to_2),
            FadeOut(arrow2_to_1),
            run_time=1.5
        )
        self.wait(0.5)
        
        # Show symmetric electron clouds
        cloud_label = Text("2. Symmetric Electron Clouds", font_size=24, color=cloud_color)
        cloud_label.move_to(pair_label.get_center())
        
        cloud_radius = 0.6
        cloud1 = Circle(radius=cloud_radius, stroke_width=2)
        cloud1.set_fill(cloud_color, opacity=0.2)
        cloud1.set_stroke(cloud_color, width=2, opacity=0.4)
        cloud1.move_to(atom1.get_center())
        
        cloud2 = Circle(radius=cloud_radius, stroke_width=2)
        cloud2.set_fill(cloud_color, opacity=0.2)
        cloud2.set_stroke(cloud_color, width=2, opacity=0.4)
        cloud2.move_to(atom2.get_center())
        
        self.play(FadeOut(pair_label), Write(cloud_label))
        self.play(FadeIn(cloud1, scale=0.5), FadeIn(cloud2, scale=0.5))
        self.wait(0.5)
        
        # Show repulsion from overlapping clouds
        repulsion_label = Text("3. Electron Repulsion", font_size=24, color=repulsion_color)
        repulsion_label.move_to(cloud_label.get_center())
        
        # Repulsion arrows pointing outward
        repulsion1 = Arrow(
            start=atom1.get_center() + RIGHT * 0.1,
            end=atom1.get_center() + LEFT * 0.6,
            color=repulsion_color,
            stroke_width=4,
            max_tip_length_to_length_ratio=0.2
        )
        
        repulsion2 = Arrow(
            start=atom2.get_center() + LEFT * 0.1,
            end=atom2.get_center() + RIGHT * 0.6,
            color=repulsion_color,
            stroke_width=4,
            max_tip_length_to_length_ratio=0.2
        )
        
        self.play(FadeOut(cloud_label), Write(repulsion_label))
        self.play(GrowArrow(repulsion1), GrowArrow(repulsion2))
        self.wait(1)
        
        # Atoms move slightly apart due to repulsion
        equilibrium_distance = 3.0
        self.play(
            atom1.animate.move_to(LEFT * equilibrium_distance / 2),
            atom2.animate.move_to(RIGHT * equilibrium_distance / 2),
            cloud1.animate.move_to(LEFT * equilibrium_distance / 2),
            cloud2.animate.move_to(RIGHT * equilibrium_distance / 2),
            FadeOut(repulsion1),
            FadeOut(repulsion2),
            run_time=1.5
        )
        self.wait(0.5)
        
        # Show equilibrium
        equilibrium_text = Text("Equilibrium Distance", font_size=24, color=GREEN)
        equilibrium_text.move_to(repulsion_label.get_center())
        
        distance_line = Line(
            atom1.get_right(),
            atom2.get_left(),
            color=WHITE,
            stroke_width=2
        )
        distance_value = MathTex(r"d_{EAM}", font_size=30, color=WHITE)
        distance_value.next_to(distance_line, DOWN, buff=0.2)
        
        self.play(
            FadeOut(repulsion_label),
            Write(equilibrium_text),
            Create(distance_line),
            Write(distance_value)
        )
        self.wait(1.5)
        
        # Clear the scene for ML potential demonstration
        self.play(
            FadeOut(atom1),
            FadeOut(atom2),
            FadeOut(cloud1),
            FadeOut(cloud2),
            FadeOut(equilibrium_text),
            FadeOut(distance_line),
            FadeOut(distance_value),
            FadeOut(eam_label)
        )
        self.wait(0.5)
        
        # ============================================================
        # PART 2: ML Potential with non-symmetric electron clouds
        # ============================================================
        
        ml_label = Text("ML Potential", font_size=30, color=PURPLE)
        ml_label.next_to(title, DOWN, buff=0.3)
        self.play(Write(ml_label))
        
        # Create two fresh atoms at initial distance
        atom1_ml = Circle(radius=atom_radius, stroke_width=3)
        atom1_ml.set_fill(gold_color, opacity=1.0)
        atom1_ml.set_stroke(GOLD_E, width=3)
        atom1_ml.move_to(LEFT * initial_distance / 2)
        
        atom2_ml = Circle(radius=atom_radius, stroke_width=3)
        atom2_ml.set_fill(gold_color, opacity=1.0)
        atom2_ml.set_stroke(GOLD_E, width=3)
        atom2_ml.move_to(RIGHT * initial_distance / 2)
        
        self.play(FadeIn(atom1_ml), FadeIn(atom2_ml))
        self.wait(0.5)
        
        # Show non-symmetric electron clouds
        asymmetric_label = Text("Non-Symmetric Electron Density", font_size=24, color=PURPLE_B)
        asymmetric_label.to_corner(UL, buff=0.5).shift(DOWN * 2)
        self.play(Write(asymmetric_label))
        
        # Create asymmetric electron clouds using ellipses
        # Atom 1: cloud stretched toward atom 2
        cloud1_ml = Ellipse(
            width=1.6,  # Wider horizontally
            height=1.0,
            stroke_width=2
        )
        cloud1_ml.set_fill(PURPLE_B, opacity=0.25)
        cloud1_ml.set_stroke(PURPLE_C, width=2, opacity=0.5)
        cloud1_ml.move_to(atom1_ml.get_center())
        cloud1_ml.shift(RIGHT * 0.15)  # Shift toward atom2
        
        # Atom 2: cloud stretched toward atom 1
        cloud2_ml = Ellipse(
            width=1.6,
            height=1.0,
            stroke_width=2
        )
        cloud2_ml.set_fill(PURPLE_B, opacity=0.25)
        cloud2_ml.set_stroke(PURPLE_C, width=2, opacity=0.5)
        cloud2_ml.move_to(atom2_ml.get_center())
        cloud2_ml.shift(LEFT * 0.15)  # Shift toward atom1
        
        self.play(
            FadeIn(cloud1_ml, scale=0.7),
            FadeIn(cloud2_ml, scale=0.7)
        )
        self.wait(1)
        
        # Show complex force vectors from ML model
        ml_force_label = Text("ML-Predicted Forces", font_size=24, color=PURPLE_C)
        ml_force_label.move_to(asymmetric_label.get_center())
        
        # Multiple force vectors at different angles (representing complex ML interactions)
        ml_arrow1 = Arrow(
            start=atom1_ml.get_center() + UP * 0.15,
            end=atom1_ml.get_center() + UP * 0.15 + RIGHT * 0.5 + UP * 0.2,
            color=PURPLE_C,
            stroke_width=3,
            max_tip_length_to_length_ratio=0.2
        )
        
        ml_arrow2 = Arrow(
            start=atom1_ml.get_center() + DOWN * 0.15,
            end=atom1_ml.get_center() + DOWN * 0.15 + RIGHT * 0.4 + DOWN * 0.15,
            color=PURPLE_D,
            stroke_width=3,
            max_tip_length_to_length_ratio=0.2
        )
        
        ml_arrow3 = Arrow(
            start=atom2_ml.get_center() + UP * 0.15,
            end=atom2_ml.get_center() + UP * 0.15 + LEFT * 0.5 + UP * 0.25,
            color=PURPLE_C,
            stroke_width=3,
            max_tip_length_to_length_ratio=0.2
        )
        
        ml_arrow4 = Arrow(
            start=atom2_ml.get_center() + DOWN * 0.15,
            end=atom2_ml.get_center() + DOWN * 0.15 + LEFT * 0.45 + DOWN * 0.1,
            color=PURPLE_D,
            stroke_width=3,
            max_tip_length_to_length_ratio=0.2
        )
        
        self.play(
            FadeOut(asymmetric_label),
            Write(ml_force_label)
        )
        self.play(
            GrowArrow(ml_arrow1),
            GrowArrow(ml_arrow2),
            GrowArrow(ml_arrow3),
            GrowArrow(ml_arrow4)
        )
        self.wait(1.5)
        
        # Atoms move to a DIFFERENT equilibrium distance
        ml_equilibrium_distance = 2.6  # Different from EAM!
        self.play(
            atom1_ml.animate.move_to(LEFT * ml_equilibrium_distance / 2 + UP * 0.15),
            atom2_ml.animate.move_to(RIGHT * ml_equilibrium_distance / 2 + UP * 0.1),
            cloud1_ml.animate.move_to(LEFT * ml_equilibrium_distance / 2 + UP * 0.15 + RIGHT * 0.1),
            cloud2_ml.animate.move_to(RIGHT * ml_equilibrium_distance / 2 + UP * 0.1 + LEFT * 0.1),
            FadeOut(ml_arrow1),
            FadeOut(ml_arrow2),
            FadeOut(ml_arrow3),
            FadeOut(ml_arrow4),
            run_time=2
        )
        self.wait(0.5)
        
        # Show ML equilibrium
        ml_equilibrium_text = Text("ML Equilibrium", font_size=24, color=PURPLE_A)
        ml_equilibrium_text.move_to(ml_force_label.get_center())
        
        ml_distance_line = Line(
            atom1_ml.get_right(),
            atom2_ml.get_left(),
            color=PURPLE_B,
            stroke_width=2
        )
        ml_distance_value = MathTex(r"d_{ML}", font_size=30, color=PURPLE_B)
        ml_distance_value.next_to(ml_distance_line, DOWN, buff=0.2)
        
        self.play(
            FadeOut(ml_force_label),
            Write(ml_equilibrium_text),
            Create(ml_distance_line),
            Write(ml_distance_value)
        )
        self.wait(1)
        
        # ============================================================
        # PART 3: Comparison
        # ============================================================
        
        comparison_label = Text("Key Difference", font_size=32, color=YELLOW)
        comparison_label.move_to(ml_label.get_center())
        
        self.play(
            FadeOut(ml_label),
            FadeOut(ml_equilibrium_text),
            Transform(title, comparison_label)
        )
        
        # Show both distances for comparison
        eam_comparison = VGroup(
            Text("EAM:", font_size=24, color=BLUE).shift(LEFT * 3 + DOWN * 2),
            Text("• Symmetric electron clouds", font_size=20, color=GRAY).shift(LEFT * 2.5 + DOWN * 2.5),
            Text("• Fixed functional form", font_size=20, color=GRAY).shift(LEFT * 2.5 + DOWN * 3)
        )
        
        ml_comparison = VGroup(
            Text("ML:", font_size=24, color=PURPLE).shift(RIGHT * 2.5 + DOWN * 2),
            Text("• Non-symmetric densities", font_size=20, color=GRAY).shift(RIGHT * 3 + DOWN * 2.5),
            Text("• Learned from data", font_size=20, color=GRAY).shift(RIGHT * 3 + DOWN * 3),
            Text("• More accurate predictions", font_size=20, color=GRAY).shift(RIGHT * 3 + DOWN * 3.5)
        )
        
        self.play(
            FadeIn(eam_comparison, shift=RIGHT * 0.5),
            FadeIn(ml_comparison, shift=LEFT * 0.5)
        )
        
        self.wait(3)
