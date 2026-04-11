# PANAIR Geometry Kernel (MATLAB NURBS)

A modern MATLAB-based NURBS geometry kernel for aerodynamic surface generation, structured meshing, and export to panel methods such as PANAIR.

## Features

### Curves (NURBSCurve)
- Evaluation, derivatives, curvature
- Arc length and inversion
- Knot insertion/removal
- Degree elevation/reduction
- Global interpolation and least-squares fitting

### Surfaces (NURBSSurface)
- Evaluation and derivatives
- Normals, curvature
- Lofting, Coons, Gordon, sweep, revolution
- Knot refinement and splitting

### Meshing
- Structured grids via isoMesh(nu, nv)
- Supports trimmed surfaces and degenerate quads

### Mesh Export
- WGS (PANAIR)
- VTK
- STL
- CSV

## Workflow

1. Fit airfoil curve
2. Split at leading edge
3. Loft upper/lower surfaces
4. Mesh surfaces
5. Export meshes

## Usage

```matlab
addpath(genpath('path_to_repo'));
demo_wing_mesh_export
```

## Notes
- Uses exact NURBS geometry
- Structured meshes for aerodynamic solvers
- Degenerate quads supported

## Author
Matthew Rombold