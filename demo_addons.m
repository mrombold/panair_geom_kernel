%% DEMO_ADDONS_VALIDATION.m
% Demonstrates the addon utilities for the panair geometry kernel:
%
% 1. geom.SurfaceEval  - analytic first/second surface derivatives
% 2. geom.PatchOps     - exact patch splitting/extraction
% 3. geom.Intersect    - curve/surface intersection
%
% The script is written in the same style as the existing repo demos and is
% intended to serve both as a usage example and a quick validation check.
%
% Expected layout:
%   repo_root/
%     +geom/
%     demo_geom_kernel.m
%     demo_new_features.m
%     demo_addons_validation.m   <-- this file
%
% If the addon classes were copied into +geom/, this script should run from
% the repo root or directly from this file's folder.

clear; close all; clc;

%% ---- PATH SETUP -------------------------------------------------------
script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir)
    addpath(script_dir);
end

required_classes = {
    'geom.NURBSCurve'
    'geom.NURBSSurface'
    'geom.Aircraft'
    'geom.SurfaceEval'
    'geom.PatchOps'
    'geom.Intersect'
    };

for k = 1:numel(required_classes)
    if isempty(which(required_classes{k}))
        error(['Could not find required class: %s\n' ...
               'Make sure the repo root is on the MATLAB path and the addon\n' ...
               'files were copied into +geom/.'], required_classes{k});
    end
end

fprintf('geom package found at: %s\n', script_dir);
fprintf('Addon classes detected: SurfaceEval, PatchOps, Intersect\n\n');

%% ======================================================================
% 1. BUILD A SIMPLE TEST WING SURFACE
%% ======================================================================
fprintf('--- 1. Build test wing surface ---\n');

% Use the same airfoil family across three span stations.
[xaf, yaf] = geom.Demo.naca4Coords(0.02, 0.4, 0.12, 80);
af = [xaf(:), yaf(:)];
airfoils = {af, af, af};
spans     = [0.0, 0.8, 1.6];
chords    = [0.60, 0.42, 0.26];
sweeps    = [0.00, 0.08, 0.18];
twists    = [0.0, -1.5, -3.0];
dihedrals = [0.00, 0.03, 0.08];

S_loft = geom.Aircraft.wingSurface(airfoils, spans, chords, sweeps, twists, dihedrals);
if isprop(S_loft, 'surface')
    S = S_loft.surface;
else
    S = S_loft;
end

fprintf(' Surface control net: %d x %d\n', S.n+1, S.m+1);
fprintf(' Surface degrees: p=%d, q=%d\n', S.p, S.q);

mesh0 = S.isoMesh(32, 10);
fig1 = figure('Name','Addon Demo - Base Wing Surface');
surf(mesh0.X, mesh0.Y, mesh0.Z, 'EdgeColor', [0.15 0.15 0.15], 'FaceColor', 'none');
axis equal; grid on; view(25, 30);
title('Base Wing Surface'); xlabel('X'); ylabel('Y'); zlabel('Z');

%% ======================================================================
% 2. ANALYTIC SURFACE DERIVATIVES AND CURVATURE
%% ======================================================================
fprintf('\n--- 2. SurfaceEval: analytic derivatives and curvature ---\n');

u0 = 0.37;
v0 = 0.58;
P0 = S.evaluate(u0, v0);
[Su_fd, Sv_fd] = S.partialDerivatives(u0, v0);
[Su, Sv, Suu, Suv, Svv] = geom.SurfaceEval.firstSecond(S, u0, v0);
[K, H, Nhat] = geom.SurfaceEval.curvatures(S, u0, v0);

fprintf(' Evaluation point: [%.6f %.6f %.6f]\n', P0);
fprintf(' ||Su(analytic)-Su(existing)|| = %.3e\n', norm(Su - Su_fd));
fprintf(' ||Sv(analytic)-Sv(existing)|| = %.3e\n', norm(Sv - Sv_fd));
fprintf(' ||Suu|| = %.3e, ||Suv|| = %.3e, ||Svv|| = %.3e\n', ...
    norm(Suu), norm(Suv), norm(Svv));
fprintf(' Gaussian curvature K = %.6e\n', K);
fprintf(' Mean curvature     H = %.6e\n', H);

fig2 = figure('Name','Addon Demo - SurfaceEval');
surf(mesh0.X, mesh0.Y, mesh0.Z, 'EdgeColor', [0.70 0.70 0.70], 'FaceAlpha', 0.15);
hold on;
plot3(P0(1), P0(2), P0(3), 'ro', 'MarkerFaceColor','r', 'MarkerSize',8);
quiver3(P0(1), P0(2), P0(3), Su(1),  Su(2),  Su(3),  0.12, 'b', 'LineWidth', 1.5);
quiver3(P0(1), P0(2), P0(3), Sv(1),  Sv(2),  Sv(3),  0.12, 'g', 'LineWidth', 1.5);
quiver3(P0(1), P0(2), P0(3), Nhat(1),Nhat(2),Nhat(3),0.12, 'm', 'LineWidth', 1.5);
axis equal; grid on; view(35, 28);
legend('Surface','Eval point','S_u','S_v','Normal','Location','best');
title('Analytic Surface Derivatives at One Point');
xlabel('X'); ylabel('Y'); zlabel('Z');

%% ======================================================================
% 3. PATCH SPLITTING AND SUBPATCH EXTRACTION
%% ======================================================================
fprintf('\n--- 3. PatchOps: exact split and extraction ---\n');

usplit = 0.45;
vsplit = 0.55;
[S_uA, S_uB] = geom.PatchOps.splitU(S, usplit);
[S_vA, S_vB] = geom.PatchOps.splitV(S, vsplit);

fprintf(' splitU at u=%.2f -> [%dx%d] + [%dx%d] control nets\n', usplit, ...
    S_uA.n+1, S_uA.m+1, S_uB.n+1, S_uB.m+1);
fprintf(' splitV at v=%.2f -> [%dx%d] + [%dx%d] control nets\n', vsplit, ...
    S_vA.n+1, S_vA.m+1, S_vB.n+1, S_vB.m+1);

ur = [0.20, 0.78];
vr = [0.18, 0.82];
S_sub = geom.PatchOps.extractUV(S, ur, vr);

% Verify exact geometry preservation by mapping original (u,v) to local [0,1]^2.
uv_test = [0.32 0.30; 0.55 0.60; 0.72 0.78];
max_sub_err = 0.0;
for k = 1:size(uv_test,1)
    u = uv_test(k,1);
    v = uv_test(k,2);
    if u < ur(1) || u > ur(2) || v < vr(1) || v > vr(2)
        continue;
    end
    p_base = S.evaluate(u, v);
    uloc = (u - ur(1)) / (ur(2) - ur(1));
    vloc = (v - vr(1)) / (vr(2) - vr(1));
    p_sub  = S_sub.evaluate(uloc, vloc);
    max_sub_err = max(max_sub_err, norm(p_base - p_sub));
end
fprintf(' Max geometry error after extractUV: %.3e\n', max_sub_err);

loops = geom.PatchOps.boundaryCurves(S_sub);
t = linspace(0, 1, 150);
pts_u0 = loops.u0.evaluate(t);
pts_u1 = loops.u1.evaluate(t);
pts_v0 = loops.v0.evaluate(t);
pts_v1 = loops.v1.evaluate(t);

mesh_sub = S_sub.isoMesh(26, 10);
fig3 = figure('Name','Addon Demo - PatchOps');
surf(mesh0.X, mesh0.Y, mesh0.Z, 'EdgeColor', [0.85 0.85 0.85], 'FaceAlpha', 0.10); hold on;
surf(mesh_sub.X, mesh_sub.Y, mesh_sub.Z, 'EdgeColor', 'k', 'FaceColor', [0.8 0.9 1.0], 'FaceAlpha', 0.65);
plot3(pts_u0(:,1), pts_u0(:,2), pts_u0(:,3), 'r-', 'LineWidth', 1.7);
plot3(pts_u1(:,1), pts_u1(:,2), pts_u1(:,3), 'r-', 'LineWidth', 1.7);
plot3(pts_v0(:,1), pts_v0(:,2), pts_v0(:,3), 'b-', 'LineWidth', 1.7);
plot3(pts_v1(:,1), pts_v1(:,2), pts_v1(:,3), 'b-', 'LineWidth', 1.7);
axis equal; grid on; view(30, 28);
title('Extracted Subpatch and Exact Boundary Curves');
legend('Base surface','Extracted subpatch','u-boundaries','', 'v-boundaries','', 'Location','best');
xlabel('X'); ylabel('Y'); zlabel('Z');

%% ======================================================================
% 4. CURVE / SURFACE INTERSECTION
%% ======================================================================
fprintf('\n--- 4. Intersect: curve-surface intersection ---\n');

% Build a line that passes through a known surface point.
u_hit = 0.62;
v_hit = 0.44;
P_hit = S.evaluate(u_hit, v_hit);
dz = 0.25;
C_line = geom.NURBSCurve([P_hit + [0 0 dz]; P_hit - [0 0 dz]], 1);

hits = geom.Intersect.curveSurface(C_line, S, ...
    'nCurve', 41, 'nU', 21, 'nV', 21, 'tolXYZ', 1e-9, 'tolParam', 1e-9, 'maxIter', 50);

fprintf(' Intersections found: %d\n', numel(hits));
if isempty(hits)
    warning('No curve/surface intersection found.');
else
    % Choose the best residual if multiple seeds collapsed to the same hit.
    [~, ibest] = min([hits.residual]);
    h = hits(ibest);
    fprintf(' Best hit residual: %.3e\n', h.residual);
    fprintf(' Expected params:  uc~%.4f (line), us~%.4f, vs~%.4f\n', 0.5, u_hit, v_hit);
    fprintf(' Returned params:  uc=%.4f, us=%.4f, vs=%.4f\n', h.uCurve, h.uSurface, h.vSurface);
    fprintf(' Point error vs planted point: %.3e\n', norm(h.point - P_hit));
end

line_pts = C_line.evaluate(linspace(0,1,100));
fig4 = figure('Name','Addon Demo - Intersect');
surf(mesh0.X, mesh0.Y, mesh0.Z, 'EdgeColor', [0.65 0.65 0.65], 'FaceAlpha', 0.15); hold on;
plot3(line_pts(:,1), line_pts(:,2), line_pts(:,3), 'k-', 'LineWidth', 2.0);
plot3(P_hit(1), P_hit(2), P_hit(3), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
if ~isempty(hits)
    for i = 1:numel(hits)
        plot3(hits(i).point(1), hits(i).point(2), hits(i).point(3), 'ro', ...
            'MarkerFaceColor', 'r', 'MarkerSize', 7);
    end
    legend('Surface','Cut line','Planted point','Intersection(s)','Location','best');
else
    legend('Surface','Cut line','Planted point','Location','best');
end
axis equal; grid on; view(28, 26);
title('Curve / Surface Intersection Demo');
xlabel('X'); ylabel('Y'); zlabel('Z');

%% ======================================================================
% 5. SUMMARY / PASS-FAIL HINTS
%% ======================================================================
fprintf('\n=== Addon Demo Summary ===\n');
fprintf(' SurfaceEval first-derivative agreement should be near machine precision.\n');
fprintf(' extractUV geometry error should be ~eps to 1e-12-ish.\n');
fprintf(' curveSurface should return at least one hit with small residual.\n');
fprintf('\nSuggested next manual checks:\n');
fprintf(' 1) Change the planted line to a swept cut curve and test multiple intersections.\n');
fprintf(' 2) Use PatchOps.extractUV() to isolate a body/wing closeout region.\n');
fprintf(' 3) Replace any finite-difference surface curvature calls with SurfaceEval.\n');