from manim import *

class IntroTitle(Scene):
    def construct(self):
        # Title - largest
        title = Text(
            "Progress Towards a Quantum-Accurate\n" 
            "Classical SNAP Machine Learning\n"
            "Interaction Potential for Gold",
            font_size=50,
            weight=BOLD,
            color=GOLD_A
        )
        title.to_edge(UP, buff=1.5)
        
        # Authors - medium size with superscript numbers
        # We'll build this with multiple Text objects positioned together
        author_parts = VGroup(
            Text("Tyrel Boese", font_size=32, color=WHITE),
            Text("[1]", font_size=20, color=WHITE).shift(0.3*UP),
            Text(", Ian Anderson", font_size=32, color=WHITE),
            Text("[2]", font_size=20, color=WHITE).shift(0.3*UP),
            Text(", James Goff", font_size=32, color=WHITE),
            Text("[3]", font_size=20, color=WHITE).shift(0.3*UP),
            Text(",", font_size=32, color=WHITE),
        ).arrange(RIGHT, buff=0.05, aligned_edge=DOWN)
        
        # Second line of authors
        author_parts2 = VGroup(
            Text("Jarrod Schiffbauer", font_size=32, color=WHITE),
            Text("[1] [4]", font_size=20, color=WHITE).shift(0.3*UP),
        ).arrange(RIGHT, buff=0.05, aligned_edge=DOWN)
        
        # Combine author lines
        authors = VGroup(author_parts, author_parts2).arrange(DOWN, buff=0.3, center=True)
        authors.next_to(title, DOWN, buff=0.5)
        
        # Affiliations - smallest with superscript numbers
        affil1_parts = VGroup(
            Text("[1]", font_size=16, color=GRAY).shift(0.2*UP),
            Text(" Colorado Mesa University Department of Physical and Environmental Sciences", font_size=24, color=GRAY),
        ).arrange(RIGHT, buff=0.05, aligned_edge=DOWN)
        
        affil2_parts = VGroup(
            Text("[2]", font_size=16, color=GRAY).shift(0.2*UP),
            Text(" Colorado State University School of Materials Science and Engineering", font_size=24, color=GRAY),
        ).arrange(RIGHT, buff=0.05, aligned_edge=DOWN)
        
        affil3_parts = VGroup(
            Text("[3]", font_size=16, color=GRAY).shift(0.2*UP),
            Text(" Sandia National Laboratories", font_size=24, color=GRAY),
        ).arrange(RIGHT, buff=0.05, aligned_edge=DOWN)
        
        affil4_parts = VGroup(
            Text("[4]", font_size=16, color=GRAY).shift(0.2*UP),
            Text(" University of Colorado Boulder", font_size=24, color=GRAY),
        ).arrange(RIGHT, buff=0.05, aligned_edge=DOWN)
        
        affiliations = VGroup(affil1_parts, affil2_parts, affil3_parts, affil4_parts).arrange(
            DOWN, buff=0.2, aligned_edge=LEFT
        )
        affiliations.next_to(authors, DOWN, buff=0.8)
        
        # Animate everything in
        self.play(Write(title), run_time=2)
        self.wait(0.5)
        self.play(Write(authors), run_time=2)
        self.wait(0.5)
        self.play(Write(affiliations), run_time=2)
        self.wait(3)
        
        # Fade out
        self.play(
            FadeOut(title),
            FadeOut(authors),
            FadeOut(affiliations),
            run_time=1
        )


# Original text data for reference
"""
"Progress Towards a Quantum-Accurate Classical SNAP Machine Learning Interaction Potential for Gold"

"Tyrel Boese[1], Ian Anderson[2], James Goff [3], Jarrod Schiffbauer [1] [4]"

"Colorado Mesa University Department of Physical and Environmental Sciences[1], Colorado State University School of Materials Science and Engineering[2], Sandia National Laboratories[3], University of Colorado Boulder[4]"
"""
