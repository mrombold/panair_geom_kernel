%% DEMO_GEOM_KERNEL.m
%  Demonstration of the Aircraft Geometry Kernel.
%
%  Exercises:
%    1. NURBS curve basics - construction, evaluation, derivatives
%    2. Conic arc (Roy Liming method)
%    3. Ruled surface between two airfoils
%    4. Lofted wing surface from multiple sections
%    5. Isoparametric mesh extraction
%    6. Surface of revolution (axisymmetric fuselage)

clear; close all; clc;

% ---- PATH SETUP -------------------------------------------------------
%  The +geom package directory must have its PARENT on the MATLAB path.
%  This block handles it automatically regardless of where you run from.
%
%  Layout expected:
%    geom_kernel/              <-- this directory on path
%      +geom/
%        NURBSCurve.m  ...
%      demo_geom_kernel.m      <-- this file
%
%  If you get "Unable to resolve the name 'geom.XXX'" you can also fix it
%  manually in the Command Window:
%    >> addpath('C:/path/to/geom_kernel')
% -----------------------------------------------------------------------
script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir)
    addpath(script_dir);
end
% Verify the package is visible
if ~exist('+geom', 'dir') && isempty(which('geom.NURBSCurve'))
    error(['Cannot find +geom package.\n' ...
           'Expected it at: %s\n' ...
           'Add the geom_kernel folder to your MATLAB path:\n' ...
           '  addpath(''%s'')'], ...
           fullfile(script_dir, '+geom'), script_dir);
end
fprintf('geom package found at: %s\n', script_dir);

%% ======================================================================
%  1.  NURBS CURVE BASICS
%% ======================================================================
fprintf('--- 1. NURBS Curve basics ---\n');

% Cubic B-spline through some 3-D control points
P = [0  0  0;
     1  2  0.5;
     2  1 -0.3;
     3  2  0.2;
     4  0  0];
p = 3;
C = geom.NURBSCurve(P, p);

% Evaluate at parameter values
u_test = linspace(0, 1, 5);
pts    = C.evaluate(u_test);
fprintf('  Evaluated at u=%s\n', mat2str(u_test,3));
fprintf('  Points:\n'); disp(pts);

% Derivatives
D1 = C.derivative(0.5, 1);
D2 = C.derivative(0.5, 2);
T  = C.tangent(0.5);
K  = C.curvature(linspace(0,1,50));

fprintf('  Tangent at u=0.5: [%.4f  %.4f  %.4f]\n', T);
fprintf('  Max curvature: %.4f\n', max(K));

% Arc length
L = C.arcLength();
fprintf('  Total arc length: %.4f\n', L);

% Plot
fig1 = figure('Name','NURBS Curve Demo');
C.plot(200, 'ShowCP', true, 'ShowKnots', true);
title('Cubic NURBS Curve');
view(3);

%% ======================================================================
%  2.  CONIC ARC (Roy Liming)
%% ======================================================================
fprintf('\n--- 2. Roy Liming Conic Arc ---\n');

% T0 and T2 must be CONVERGENT - their tangent lines must intersect
% to define the shoulder point P1.
% T0=[1,1,0] and T2=[-1,1,0] meet at P1=[1,1,0]: a clean symmetric arch.
% (Anti-parallel tangents like [0,1,0] / [0,-1,0] produce parallel lines
%  that never meet - rank-deficient system.)
P0 = [0  0  0];
P2 = [2  0  0];
T0 = [1  1  0];   % northeast at P0
T2 = [-1 1  0];   % northwest at P2  --> shoulder at [1,1,0]

fig2 = figure('Name','Conic Arcs');
hold on; grid on; axis equal;
colors = lines(4);
rho_vals = [0.3, 0.5, 0.7, 0.9];
labels   = {};
for k = 1:4
    Ca = geom.Aircraft.conicArc(P0, P2, T0, T2, rho_vals(k));
    Ca.plot(100, 'ShowCP', false, 'Color', colors(k,:));
    labels{k} = sprintf('\\rho = %.1f', rho_vals(k));
end
legend(labels, 'Location','north');
title('Roy Liming Conic Arcs (varying \rho)');
xlabel('X'); ylabel('Y');

fprintf('  Conic arcs plotted for rho = %s\n', mat2str(rho_vals));

%% ======================================================================
%  3.  NURBS SURFACE - Simple bicubic patch
%% ======================================================================
fprintf('\n--- 3. Bicubic NURBS Surface Patch ---\n');

% 4x4 control net: simple saddle-like patch
[ug, vg] = meshgrid(linspace(-1,1,4), linspace(-1,1,4));
Px = ug';
Py = vg';
Pz = Px.^2 - Py.^2;   % saddle

P3d = zeros(4, 4, 3);
P3d(:,:,1) = Px;
P3d(:,:,2) = Py;
P3d(:,:,3) = Pz * 0.3;

S = geom.NURBSSurface(P3d, 3, 3);

% Evaluate normal at center
N = S.normal(0.5, 0.5);
fprintf('  Normal at (0.5,0.5): [%.4f  %.4f  %.4f]\n', N);

% Extract isoparametric curve at v=0.5
iso_u = S.isoCurveU(0.5);
fprintf('  isoCurveU(0.5): %d control pts, degree %d\n', ...
        size(iso_u.P,1), iso_u.p);

fig3 = figure('Name','NURBS Surface Patch');
S.plot(30, 30, 'ShowCP', true, 'ShowIso', true, 'Alpha', 0.7);
title('Bicubic NURBS Surface (saddle patch)');
hold on;
iso_u.plot(50, 'ShowCP', false, 'Color', 'r', 'LineWidth', 2);
text(0, 0, max(Pz(:))*0.3+0.05, 'v=0.5 iso', 'Color','r');

%% ======================================================================
%  4.  WING SURFACE via LoftedSurface
%% ======================================================================
fprintf('\n--- 4. Lofted Wing Surface ---\n');

% NACA 4-series airfoil generator (normalized, unit chord)
naca4 = @(m_pc, p_pc, t_pct, N) geom.Demo.naca4Coords(m_pc, p_pc, t_pct, N);

% Wing definition
spans     = [0, 0.5, 1.0, 1.5, 2.0];   % m
chords    = [0.4, 0.38, 0.33, 0.26, 0.18];  % m (taper)
sweeps    = [0, 0.04, 0.09, 0.15, 0.22]; % m LE sweep
twists    = [0, -0.5, -1.0, -2.0, -3.0];  % deg washout
dihedrals = [0, 0.03, 0.07, 0.13, 0.20]; % m

% NACA 2412 root -> NACA 0009 tip (linear taper in camber)
n_af_pts = 80;   % more points = smoother LE; exact NURBS interp through each
airfoils = cell(5,1);
camber_m  = linspace(0.02, 0,    5);
camber_p  = linspace(0.40, 0.40, 5);
thickness = linspace(0.12, 0.09, 5);
for k = 1:5
    [x,y] = geom.Demo.naca4Coords(camber_m(k), camber_p(k), thickness(k), n_af_pts);
    airfoils{k} = [x(:), y(:)];
end

S_wing = geom.Aircraft.wingSurface(airfoils, spans, chords, sweeps, twists, dihedrals);

% Extract isoparametric mesh
mesh_w = S_wing.isoMesh(60, 20);
fprintf('  Wing mesh: %d x %d = %d quad cells\n', ...
        mesh_w.nu-1, mesh_w.nv-1, (mesh_w.nu-1)*(mesh_w.nv-1));

fig4 = figure('Name','Lofted Wing Surface');
S_wing.plot(60, 20, 'ShowIso', true, 'FaceColor', [0.8 0.85 0.95], 'Alpha', 0.9);
title('Lofted Wing: NACA 2412 -> 0009');
view(30, 25);

%% ======================================================================
%  5.  ISOPARAMETRIC MESH EXPORT (for PanAir-style use)
%% ======================================================================
fprintf('\n--- 5. Isoparametric Mesh ---\n');

nu = 20; nv = 8;
mesh_fine = S_wing.isoMesh(nu, nv);

fprintf('  Nodes:     %d\n', nu*nv);
fprintf('  Quad cells: %d\n', (nu-1)*(nv-1));
fprintf('  Connectivity array: %dx%d\n', size(mesh_fine.connectivity));

fig5 = figure('Name','Wing Quad Mesh');
surf(mesh_fine.X, mesh_fine.Y, mesh_fine.Z, ...
     'FaceColor','none', 'EdgeColor', [0.2 0.2 0.8]);
hold on;
quiver3(mesh_fine.X, mesh_fine.Y, mesh_fine.Z, ...
        0.015*mesh_fine.normals(:,:,1), ...
        0.015*mesh_fine.normals(:,:,2), ...
        0.015*mesh_fine.normals(:,:,3), 0, 'r');
axis equal; grid on;
title(sprintf('Wing Quad Mesh (%dx%d) with Surface Normals', nu, nv));
xlabel('X (chordwise)'); ylabel('Y (span)'); zlabel('Z');
view(30, 35);

%% ======================================================================
%  6.  SURFACE OF REVOLUTION (axisymmetric fuselage)
%% ======================================================================
fprintf('\n--- 6. Surface of Revolution ---\n');

% Profile: half-body of revolution (Sears-Haack-inspired)
x_fus = linspace(0, 3, 15)';
r_fus = 0.2 * sin(pi * x_fus/3) .* (1 - 0.15*x_fus/3);
profile_3d = [x_fus, zeros(size(x_fus)), r_fus];  % profile in XZ plane

S_fus = geom.Aircraft.revolutionSurface(profile_3d, [0 0 0], [1 0 0]);

mesh_fus = S_fus.isoMesh(20, 24);

fig6 = figure('Name','Fuselage Surface of Revolution');
S_fus.plot(20, 24, 'FaceColor', [0.9 0.85 0.7], 'Alpha', 0.85);
title('Axisymmetric Fuselage (Surface of Revolution)');
view(20, 25);

fprintf('  Fuselage mesh: %dx%d = %d quads\n', ...
        mesh_fus.nu-1, mesh_fus.nv-1, (mesh_fus.nu-1)*(mesh_fus.nv-1));

%% ======================================================================
%  Summary
%% ======================================================================
fprintf('\n=== Geometry Kernel Demo Complete ===\n');
fprintf('Figures: (1) NURBS Curve  (2) Conic Arcs  (3) NURBS Surface\n');
fprintf('         (4) Lofted Wing   (5) Quad Mesh   (6) Fuselage\n');