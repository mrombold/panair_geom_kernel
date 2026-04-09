%% DEMO_NEW_FEATURES.m
%  Demonstrates features added in v0.2 of the geometry kernel:
%
%   1. CompositeCurve - closed fuselage section from conic arcs
%   2. Cosine-spaced isoMesh (better airfoil LE resolution)
%   3. Surface refine() - knot insertion without geometry change
%   4. Projection.toCurve / toSurface - closest-point Newton solver
%   5. NACA 5-series airfoil
%   6. Symmetric wing mirror
%   7. MeshWriter - WGS (PanAir), VTK, STL export

clear; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir), addpath(script_dir); end
fprintf('geom package: %s\n\n', script_dir);

%% ======================================================================
%  1.  CompositeCurve -- closed fuselage cross-section
%% ======================================================================
fprintf('--- 1. Closed Fuselage Section (CompositeCurve) ---\n');

% Oval fuselage cross-section at X=1.0
% Wider than tall, rounded top, flatter bottom
CC = geom.Aircraft.closedFuselageSection(1.0, 0.3, [0.18, -0.12], ...
                                          0.65, 0.55, 0.60);

% Convert to single NURBS for lofting
C_nurbs = CC.toNURBS(3, 60);

fprintf('  Arc length: %.4f m\n', CC.arcLength());
fprintf('  Segments:   %d conic arcs\n', CC.nseg);

fig1 = figure('Name','Closed Fuselage Section');
CC.plot(300, 'ShowBreaks', true, 'Color', [0.1 0.4 0.9]);
C_nurbs.plot(300, 'ShowCP', false, 'Color', [0.9 0.3 0.1]);
legend('Composite (conic arcs)', 'Single NURBS fit', 'Location','best');
title('Closed Conic Fuselage Section');
view(2); axis equal;

%% ======================================================================
%  2.  Exact NURBS circle via CompositeCurve.circle
%% ======================================================================
fprintf('\n--- 2. Exact NURBS Circle ---\n');

CC_circ = geom.CompositeCurve.circle([0 0 0], 1.0, [0 0 1], 4);
% Verify: sample points should all be at radius 1
t_samp = linspace(0, 1, 200);
pts    = CC_circ.evaluate(t_samp);
radii  = sqrt(sum(pts(:,1:2).^2, 2));
fprintf('  Circle radius error: max = %.2e  (exact NURBS representation)\n', ...
        max(abs(radii - 1)));

fig2 = figure('Name','Exact NURBS Circle');
CC_circ.plot(200, 'Color', [0.2 0.6 0.2], 'ShowBreaks', true);
axis equal; title('Exact NURBS Circle (4 rational quadratic arcs)');
view(2);

%% ======================================================================
%  3.  Cosine-spaced mesh vs uniform
%% ======================================================================
fprintf('\n--- 3. Cosine vs Uniform Mesh Spacing ---\n');

% Build a simple swept wing
[x4, y4] = geom.Demo.naca4Coords(0.02, 0.4, 0.12, 80);
af = {[x4(:), y4(:)], [x4(:), y4(:)]};
S_test = geom.Aircraft.wingSurface(af, [0, 1.5], [0.4, 0.2], [0, 0.1], [0, -2], [0, 0.05]);

mesh_unif   = S_test.isoMesh(30, 10);
mesh_cosine = S_test.isoMesh(30, 10, 'SpacingU', 'cosine');
mesh_wrap   = S_test.isoMesh(30, 10, 'SpacingU', 'cosine_wrap');

% Compare chordwise u-parameter distribution at midspan
% The airfoil wraps TE(u=0) -> upper -> LE(u=0.5) -> lower -> TE(u=1)
% so we check spacing near u=0 (TE) and u=0.5 (LE)
u_unif   = mesh_unif.u(:, 5);
u_cosine = mesh_cosine.u(:, 5);
u_wrap   = mesh_wrap.u(:, 5);

du_te_unif   = min(abs(diff(u_unif(1:5))));       % near TE
du_le_unif   = min(abs(diff(u_unif(13:17))));      % near LE (midpoint)
du_te_cosine = min(abs(diff(u_cosine(1:5))));
du_le_cosine = min(abs(diff(u_cosine(13:17))));
du_te_wrap   = min(abs(diff(u_wrap(1:5))));
du_le_wrap   = min(abs(diff(u_wrap(13:17))));

fprintf('  Spacing type    min du @ TE   min du @ LE\n');
fprintf('  uniform         %.4f        %.4f\n', du_te_unif,   du_le_unif);
fprintf('  cosine          %.4f        %.4f  (only TE refined)\n', du_te_cosine, du_le_cosine);
fprintf('  cosine_wrap     %.4f        %.4f  (TE + LE refined)\n', du_te_wrap,   du_le_wrap);

fig3 = figure('Name','Mesh Spacing Comparison');
subplot(1,3,1);
surf(mesh_unif.X,   mesh_unif.Y,   mesh_unif.Z,   'EdgeColor','k','FaceColor','none');
title('Uniform'); axis equal; view(5,30); xlabel('X'); zlabel('Z');
subplot(1,3,2);
surf(mesh_cosine.X, mesh_cosine.Y, mesh_cosine.Z, 'EdgeColor','b','FaceColor','none');
title('cosine (TE only)'); axis equal; view(5,30); xlabel('X'); zlabel('Z');
subplot(1,3,3);
surf(mesh_wrap.X,   mesh_wrap.Y,   mesh_wrap.Z,   'EdgeColor','r','FaceColor','none');
title('cosine\_wrap (TE + LE)'); axis equal; view(5,30); xlabel('X'); zlabel('Z');

%% ======================================================================
%  4.  Surface refine() -- verify geometry invariance
%% ======================================================================
fprintf('\n--- 4. Surface Knot Refinement ---\n');

S_base = S_test.surface;

% Insert midpoint knots in both directions
Xu = [0.25, 0.5, 0.75];
Xv = [0.33, 0.67];
S_ref = S_base.refine(Xu, Xv);

% Evaluate same point on both - should be identical
test_pts = [0.3, 0.6; 0.7, 0.2; 0.5, 0.5];
max_err = 0;
for k = 1:size(test_pts,1)
    p1 = S_base.evaluate(test_pts(k,1), test_pts(k,2));
    p2 = S_ref.evaluate( test_pts(k,1), test_pts(k,2));
    max_err = max(max_err, norm(p1-p2));
end
fprintf('  Max geometry error after refinement: %.2e  (should be ~eps)\n', max_err);
fprintf('  Control pts before: %dx%d\n', S_base.n+1, S_base.m+1);
fprintf('  Control pts after:  %dx%d\n', S_ref.n+1,  S_ref.m+1);

%% ======================================================================
%  5.  Projection -- closest point
%% ======================================================================
fprintf('\n--- 5. Closest-Point Projection ---\n');

% Project a cloud of random points onto the wing surface
rng(42);
query_pts = [0.2 + 0.2*rand(5,1), ...
             0.3 + 0.9*rand(5,1), ...
             0.05 + 0.15*rand(5,1)];

[u_proj, v_proj, pts_proj, dists] = geom.Projection.toSurfaceBatch(S_base, query_pts, 12, 8);

fprintf('  Projected %d points onto wing surface:\n', size(query_pts,1));
for k = 1:size(query_pts,1)
    fprintf('    Query [%.3f %.3f %.3f] -> Surface [%.3f %.3f %.3f]  d=%.4f\n', ...
            query_pts(k,:), pts_proj(k,:), dists(k));
end

fig5 = figure('Name','Surface Projection');
S_base.plot(30, 20, 'Alpha', 0.6, 'EdgeAlpha', 0.1);
hold on;
plot3(query_pts(:,1), query_pts(:,2), query_pts(:,3), 'r*', 'MarkerSize',10);
plot3(pts_proj(:,1),  pts_proj(:,2),  pts_proj(:,3),  'bo', 'MarkerSize',8,'MarkerFaceColor','b');
for k = 1:size(query_pts,1)
    plot3([query_pts(k,1),pts_proj(k,1)],[query_pts(k,2),pts_proj(k,2)],[query_pts(k,3),pts_proj(k,3)],'k--');
end
legend('Surface','Query pts','Projected pts','Location','best');
title('Closest-Point Projection onto Wing Surface');
view(20,30);

%% ======================================================================
%  6.  NACA 5-series airfoil
%% ======================================================================
fprintf('\n--- 6. NACA 5-Series Airfoil ---\n');

des5 = 23012;
[x5, y5] = geom.Aircraft.naca5Coords(des5, 80);
C5 = geom.Aircraft.airfoilFromCoords([x5(:), y5(:)]);

fprintf('  NACA %d: %d data points, %d control pts\n', des5, numel(x5), size(C5.P,1));

% Compare with NACA 2412 from 4-series
[x4, y4] = geom.Demo.naca4Coords(0.02, 0.4, 0.12, 80);
C4 = geom.Aircraft.airfoilFromCoords([x4(:), y4(:)]);

fig6 = figure('Name','NACA 5-Series vs 4-Series');
hold on; grid on; axis equal;
u5 = linspace(0,1,300);
pts5 = C5.evaluate(u5);
pts4 = C4.evaluate(u5);
plot(pts4(:,1), pts4(:,3), 'b-', 'LineWidth', 1.5);
plot(pts5(:,1), pts5(:,3), 'r-', 'LineWidth', 1.5);
legend('NACA 2412 (4-series)', sprintf('NACA %d (5-series)', des5), 'Location','best');
title('NACA Airfoil Comparison');
xlabel('x/c'); ylabel('z/c');

%% ======================================================================
%  7.  Symmetric wing + MeshWriter export
%% ======================================================================
fprintf('\n--- 7. Symmetric Wing + Export ---\n');

% Build half-wing then mirror
nsec = 4;
spans_h   = linspace(0, 2.0, nsec);
chords_h  = linspace(0.5, 0.15, nsec);
sweeps_h  = linspace(0, 0.35, nsec);
twists_h  = linspace(0, -3, nsec);
dihedral_h= linspace(0, 0.12, nsec);

afs = cell(nsec,1);
for k = 1:nsec
    t_k = 0.12 - 0.03*(k-1)/(nsec-1);
    [xk,yk] = geom.Demo.naca4Coords(0.02, 0.4, t_k, 60);
    afs{k}  = [xk(:), yk(:)];
end

S_right = geom.Aircraft.wingSurface(afs, spans_h, chords_h, sweeps_h, twists_h, dihedral_h);
[~, S_left] = geom.Aircraft.symmetricWing(S_right);

fig7 = figure('Name','Full Span Wing');
hold on;
S_right.surface.plot(30, 12, 'FaceColor',[0.7 0.8 1.0], 'Alpha',0.9, 'ShowIso',false);
S_left.plot(30, 12, 'FaceColor',[0.7 0.8 1.0], 'Alpha',0.9, 'ShowIso',false);
title('Full-Span Wing (half + mirror)');
view(0, 60); axis equal; grid on;

% Export to WGS (PanAir)
mesh_r = S_right.isoMesh(20, 8, 'SpacingU','cosine');
mesh_l = S_left.isoMesh(20, 8, 'SpacingU','cosine');
net_r  = geom.MeshWriter.meshToNetwork(mesh_r, 'Wing_Stbd', 0, 0);
net_l  = geom.MeshWriter.meshToNetwork(mesh_l, 'Wing_Port', 0, 0);
geom.MeshWriter.toWGS({net_r, net_l}, fullfile(script_dir, 'wing_full.wgs'));

% Export to VTK (ParaView)
geom.MeshWriter.toVTK({mesh_r, mesh_l}, fullfile(script_dir, 'wing_full.vtk'));

% Export to STL
geom.MeshWriter.toSTL({mesh_r, mesh_l}, fullfile(script_dir, 'wing_full.stl'));

fprintf('\n=== New Features Demo Complete ===\n');
fprintf('Export files written to: %s\n', script_dir);
fprintf('  wing_full.wgs  (PanAir LaWGS)\n');
fprintf('  wing_full.vtk  (ParaView)\n');
fprintf('  wing_full.stl  (STL mesh)\n');