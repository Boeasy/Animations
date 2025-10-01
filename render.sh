#!/usr/bin/env bash
# Simple wrapper to call manim with sane defaults
SCENE=${1:-HelloWorld}
FILE=${2:-examples/simple_scene.py}
# Default to preview, low quality
conda run -n math360 python -m manim -pql "$FILE" "$SCENE"
