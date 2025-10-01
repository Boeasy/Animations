# Manim setup for this project

This repository contains a minimal scaffold to create and render Manim Community animations on macOS.

## What's included

- `requirements.txt` - Python deps for Manim Community
- `.gitignore` - common ignores for Python/manim
- `examples/simple_scene.py` - a tiny example Manim scene you can render
- `examples/molecular_dynamics.py` - **6-particle molecular dynamics with Lennard-Jones potential**
- `render.sh` - small helper script to run the Manim CLI
- `.vscode/tasks.json` - a VS Code Task to run the renderer

## Quick Examples

**Render the simple HelloWorld scene:**
```bash
conda run -n math360 python -m manim -pql examples/simple_scene.py HelloWorld
```

**Render the molecular dynamics simulation:**
```bash
# Basic version - 6 particles
conda run -n math360 python -m manim -pql examples/molecular_dynamics.py MolecularDynamics

# With particle trails
conda run -n math360 python -m manim -pql examples/molecular_dynamics.py MolecularDynamicsWithTrails

# Two-particle collision with live energy visualization ⭐
conda run -n math360 python -m manim -pql examples/molecular_dynamics.py TwoParticleCollision
```

Quick start (recommended using a virtual environment)

1) Install system dependencies (Homebrew)

   # If you don't have Homebrew, install it from https://brew.sh/
   brew update
   brew install python ffmpeg pkg-config cairo pango libffi

   # On Apple Silicon (M1/M2) you may need to ensure brew's PATH is set correctly. Example:
   # echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
   # eval "$(/opt/homebrew/bin/brew shellenv)"

2) Create and activate a virtual environment

   python3 -m venv .venv
   source .venv/bin/activate

3) Install Python dependencies

   pip install --upgrade pip
   pip install -r requirements.txt

4) Render the example

   # make the script executable once
   chmod +x ./render.sh
   # run the default example
   ./render.sh

The renderer script is a thin wrapper around `python -m manim` and uses the `-pql` (preview, low quality) flags by default. To render a different scene or change quality, pass arguments to `render.sh`.

Troubleshooting
- If manim reports missing system libraries, ensure Homebrew packages installed above are present and on your PATH.
- If you prefer conda: create a conda env and `pip install -r requirements.txt` inside it.

Extra troubleshooting
- If importing manim fails after pip install, ensure your virtualenv is activated and `python -m pip show manim` lists version 0.19.0 (or the version you installed).
- If rendering fails with Cairo/Pango errors, confirm `pkg-config --cflags cairo` returns values and that Homebrew's pkg-config is on PATH.

Notes on versions
- This scaffold targets Manim Community (manim) 0.17+ API. If you need a different version, update `requirements.txt`.
