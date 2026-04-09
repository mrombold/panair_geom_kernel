# Aircraft Geometry Kernel — MATLAB

A NURBS-based geometry kernel for aircraft surface modeling, isoparametric meshing,
and preparation of aerodynamic analysis grids (PanAir, mesh generators, etc.).

## Package Structure

```
geom_kernel/
  +geom/
    BasisFunctions.m      % Core NURBS algorithms (Piegl & Tiller)
    NURBSCurve.m          % Degree-p NURBS curve
    NURBSSurface.m        % Bi-degree NURBS surface
    LoftedSurface.m       % Skinning / lofting through section curves
    Aircraft.m            % Factory: wing, fuselage, conic arcs, SOR
    Demo.m                % NACA 4-digit airfoil generator + helpers
  demo_geom_kernel.m      % Full demonstration script
```

## Quick Start

```matlab
addpath('path/to/geom_kernel')   % adds +geom package to path

% NURBS curve
P = [0 0 0; 1 2 0; 2 1 0; 3 0 0];
C = geom.NURBSCurve(P, 3);
pts = C.evaluate(linspace(0,1,100));

% NURBS surface
P_net = zeros(4,4,3);  % fill in control net...
S = geom.NURBSSurface(P_net, 3, 3);
mesh = S.isoMesh(30, 30);

% Conic arc (Roy Liming method)
Ca = geom.Aircraft.conicArc([0 0 0], [2 0 0], [0 1 0], [0 -1 0], 0.5);

% Lofted wing surface
S_wing = geom.Aircraft.wingSurface(airfoil_cell, spans, chords, sweeps, twists, dihedrals);

% Isoparametric mesh
mesh = S_wing.isoMesh(40, 12);
% mesh.X, .Y, .Z          — [nu x nv] coordinate grids
% mesh.normals             — [nu x nv x 3] surface normals
% mesh.connectivity        — [(nu-1)*(nv-1) x 4] quad connectivity
```

## Algorithms (Piegl & Tiller Reference)

| Algorithm | Method | P&T Reference |
|-----------|--------|---------------|
| A2.1 | `BasisFunctions.FindSpan` | Knot span binary search |
| A2.2 | `BasisFunctions.BasisFuns` | Non-zero B-spline basis functions |
| A2.3 | `BasisFunctions.DersBasisFuns` | Derivatives of basis functions |
| A3.1 | `NURBSCurve.evaluate` | NURBS curve point |
| A3.2 | `NURBSCurve.derivative` | NURBS curve derivatives |
| A4.3 | `NURBSSurface.evaluate` | NURBS surface point |
| A4.4 | `NURBSSurface.partialDerivatives` | Rational surface derivatives |
| A5.4 | `NURBSCurve.refine` | Knot refinement |
| §10.3 | `LoftedSurface` | Surface skinning |

## Road Map / Next Steps

### Near Term
- [ ] `NURBSCurve.elevate(t)` — degree elevation
- [ ] `NURBSSurface.refine(Xu, Xv)` — surface knot refinement
- [ ] Closed/periodic NURBS curves (for fuselage sections)
- [ ] `geom.Intersect` — curve-surface and surface-surface intersection
- [ ] Global surface interpolation (not just skinning)

### Medium Term
- [ ] `geom.IsoMesh` — structured quadrilateral mesher with cosine spacing
- [ ] PanAir WGS format writer (`mesh_to_wgs.m`)
- [ ] SU2 `.su2` mesh writer
- [ ] STL export
- [ ] Curvature-adaptive mesh refinement

### Aircraft-Specific
- [ ] Wing-body junction intersection curve
- [ ] Control surface (flap/aileron) splitting
- [ ] Nacelle/pylon geometry
- [ ] Full closed fuselage from conic sections (periodic NURBS)

## Design Notes

### Coordinate Convention
- X: streamwise (aft positive)
- Y: spanwise (right wing positive)
- Z: up positive

### Parameter Convention
- All NURBS parameters normalized to [0, 1]
- `u` = chordwise / profile direction
- `v` = spanwise / axial direction

### Compatibility with PanAir
The `isoMesh()` output is structured to feed directly into a WGS writer:
- Each surface becomes one or more networks
- Normal vectors point outward (into flow)
- Quad connectivity follows PanAir network row/column convention

## Dependencies
- MATLAB R2018b or later (for `inputParser`, `classdef`, handle classes)
- No toolboxes required
- Optional: Parallel Computing Toolbox for large mesh generation

## References
- Piegl, L. & Tiller, W. (1997). *The NURBS Book* (2nd ed.). Springer.
- Liming, R.A. (1944). *Practical Analytic Geometry with Applications to Aircraft*. Macmillan.
- Farin, G. (2002). *Curves and Surfaces for CAGD* (5th ed.). Morgan Kaufmann.
