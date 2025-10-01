from manim import *


class HelloWorld(Scene):
    def construct(self):
        text = Text("Hello, Manim!", font_size=72)
        self.play(Write(text))
        self.wait(1)
