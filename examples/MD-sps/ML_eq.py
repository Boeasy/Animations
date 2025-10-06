"""Machine Learning Potential Visualization"""

from manim import *
import numpy as np


class MLPotentialIntro(Scene):
    def construct(self):
        # Title
        title = Text("Machine Learning Interatomic Potentials", font_size=36)
        title.to_edge(UP)
        self.play(Write(title))
        self.wait(0.5)
        
        # Create nucleus (center of atom)
        nucleus = Dot(radius=0.15, color=GOLD_A)
        nucleus.shift(2*LEFT)
        
        # Create symmetrical electron cloud (circular)
        # Using multiple concentric circles with varying opacity
        electron_cloud_symmetrical = VGroup()
        for i in range(8):
            radius = 0.3 + i * 0.15
            opacity = 0.8 - i * 0.08
            circle = Circle(
                radius=radius,
                color=GOLD_B,
                stroke_width=2,
                fill_opacity=opacity * 0.3,
                stroke_opacity=opacity
            )
            circle.move_to(nucleus.get_center())
            electron_cloud_symmetrical.add(circle)
        
        # Add some electron dots orbiting
        electrons = VGroup()
        for angle in [0, PI/2, PI, 3*PI/2]:
            electron = Dot(radius=0.06, color=GOLD_C)
            electron.move_to(nucleus.get_center() + 1.2 * np.array([np.cos(angle), np.sin(angle), 0]))
            electrons.add(electron)
        
        atom_label = Text("Atomic Orbitals\n(Spherical Symmetry)", font_size=20)
        atom_label.next_to(nucleus, DOWN, buff=1.5)
        
        # Machine Learning Potential equation - show with first orbital
        ml_equation = MathTex(r"\sum_{i}U_{i}=\sum_{i}\sum_{\alpha}c_{\alpha}B_{\alpha}^{(i)}", font_size=44)
        ml_equation.shift(3*RIGHT + 0.5*UP)
        
        equation_label = Text("ML Potential Energy", font_size=24, color=GOLD)
        equation_label.next_to(ml_equation, UP, buff=0.5)
        
        # Term explanations
        term_explanations = VGroup(
            MathTex(r"i", font_size=28, color=BLUE).next_to(MathTex(r"\text{: atom index}", font_size=24)),
            MathTex(r"\alpha", font_size=28, color=GREEN).next_to(MathTex(r"\text{: descriptor index}", font_size=24)),
            MathTex(r"B_{\alpha}^{(i)}", font_size=28, color=YELLOW).next_to(MathTex(r"\text{: descriptor value for atom } i", font_size=24)),
            MathTex(r"c_{\alpha}", font_size=28, color=PURPLE).next_to(MathTex(r"\text{: learned coefficient}", font_size=24)),
        )
        
        # Better formatting for explanations
        explanations = VGroup(
            VGroup(
                MathTex(r"i", font_size=26, color=BLUE),
                Text(": atom index", font_size=20)
            ).arrange(RIGHT, buff=0.2),
            VGroup(
                MathTex(r"\alpha", font_size=26, color=GREEN),
                Text(": descriptor index", font_size=20)
            ).arrange(RIGHT, buff=0.2),
            VGroup(
                MathTex(r"B_{\alpha}^{(i)}", font_size=26, color=YELLOW),
                Text(": descriptor value for atom ", font_size=20),
                MathTex(r"i", font_size=20, color=BLUE)
            ).arrange(RIGHT, buff=0.1),
            VGroup(
                MathTex(r"c_{\alpha}", font_size=26, color=PURPLE),
                Text(": learned coefficient (weight)", font_size=20)
            ).arrange(RIGHT, buff=0.2),
        ).arrange(DOWN, buff=0.3, aligned_edge=LEFT)
        explanations.next_to(ml_equation, DOWN, buff=0.8)
        explanations.shift(0.3*LEFT)
        
        # Show the symmetrical atom with equation
        self.play(
            Create(nucleus),
            Create(electron_cloud_symmetrical),
            Create(electrons),
            run_time=2
        )
        self.play(Write(atom_label), run_time=1)
        self.wait(0.5)
        
        # Bring in equation and label
        self.play(
            FadeIn(ml_equation, shift=LEFT),
            Write(equation_label),
            run_time=2
        )
        self.wait(0.5)
        
        # Show explanations
        self.play(Write(explanations), run_time=2.5)
        self.wait(2)
        
        # Create asymmetrical electron cloud (elliptical/distorted)
        electron_cloud_asymmetrical = VGroup()
        for i in range(8):
            radius_x = 0.3 + i * 0.18  # Stretched in x
            radius_y = 0.3 + i * 0.12  # Compressed in y
            opacity = 0.8 - i * 0.08
            
            # Create ellipse
            ellipse = Ellipse(
                width=radius_x * 2,
                height=radius_y * 2,
                color=GOLD_E,
                stroke_width=2,
                fill_opacity=opacity * 0.3,
                stroke_opacity=opacity
            )
            ellipse.move_to(nucleus.get_center())
            ellipse.rotate(PI/6)  # Rotate slightly for asymmetry
            electron_cloud_asymmetrical.add(ellipse)
        
        # Asymmetric electron positions (all clustered on the back/left side)
        electrons_asymmetric = VGroup()
        asymmetric_angles = [PI/2, 3*PI/4, 5*PI/6, PI]
        asymmetric_radii = [0.95, 0.9, 0.95, 1.1]
        asymmetric_positions = [
            r * np.array([np.cos(angle), np.sin(angle), 0])
            for r, angle in zip(asymmetric_radii, asymmetric_angles)
        ]
        for pos in asymmetric_positions:
            electron = Dot(radius=0.06, color=GOLD_D)
            electron.move_to(nucleus.get_center() + pos)
            electrons_asymmetric.add(electron)
        
        new_label = Text("Atomic Orbitals\n(Asymmetrical)", font_size=20, color=GOLD_D)
        new_label.next_to(nucleus, DOWN, buff=1.5)
        
        # Morph electron cloud
        self.play(
            Transform(electron_cloud_symmetrical, electron_cloud_asymmetrical),
            Transform(electrons, electrons_asymmetric),
            nucleus.animate.set_color(GOLD_E),
            Transform(atom_label, new_label),
            run_time=2.5
        )
        self.wait(2)
        
        # Fade out explanations
        self.play(FadeOut(explanations), run_time=1)
        self.wait(0.5)
        
        # Fade out
        self.play(
            *[FadeOut(mob) for mob in [
                title, nucleus, electron_cloud_symmetrical, electrons,
                atom_label, ml_equation, equation_label
            ]],
            run_time=1.5
        )
