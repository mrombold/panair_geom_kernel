%% demo_patch_tangent_continuity
clear; close all; clc;

if ~exist('+geom', 'dir') && isempty(which('geom.NURBSSurface'))
    error('Run this from the repo root or add the geometry kernel to the MATLAB path.');
end

fprintf('=== Patch tangent continuity demo ===\n');

%% 1) Build a simple lofted half-body patch with a symmetry edge at v0
stations = [0 5 10 15 20];
curvesA = cell(1, numel(stations));
curvesB = cell(1, numel(stations));

for k = 1:numel(stations)
    x = stations(k);
    rA = 1.0 + 0.15*cos(pi*x/20);
    rB = 0.65 + 0.10*sin(pi*x/20);

    % Patch A: crown to shoulder on right side
    PA = [x 0.0 1.2; x 0.4*rA 1.1; x 0.8*rA 0.8; x 1.0*rA 0.45];
    curvesA{k} = geom.NURBSCurve.globalInterp(PA, 3, 'centripetal');

    % Patch B: shoulder to lower side on right side
    PB = [x 1.0*rA 0.45; x 1.15*rA 0.1; x 1.10*rA -0.2; x 0.75*rA -0.7];
    curvesB{k} = geom.NURBSCurve.globalInterp(PB, 3, 'centripetal');
end

SA = geom.NURBSSurface.loft(curvesA, 3, 'centripetal');
SB = geom.NURBSSurface.loft(curvesB, 3, 'centripetal');

patchA = geom.Patch('Upper', SA);
patchB = geom.Patch('Lower', SB);

%% 2) Enforce symmetry tangency at the centerline edge of patchA
% Here the symmetry plane is y = 0 and the crown edge is v0.
repSym = patchA.enforceSymmetryPlane('v0', [0 0 0], [0 1 0]);
fprintf('Symmetry edge max plane distance      = %.3e\n', repSym.maxEdgePlaneDistance);
fprintf('Symmetry edge cross-derivative error  = %.3e\n', repSym.maxCrossDerivativeOffNormal);

%% 3) Enforce G1 between patchA (v1 edge) and patchB (v0 edge)
repG1_before = patchA.g1ReportWithPatch('v1', patchB, 'v0');
fprintf('Before G1: max position err = %.3e, cross-plane mismatch = %.3e\n', ...
    repG1_before.maxPositionError, repG1_before.maxCrossPlaneMismatch);

repG1_after = patchB.enforceG1WithPatch('v0', patchA, 'v1', 'Lambda', 1.0, 'AverageEdge', true);
fprintf('After  G1: max position err = %.3e, cross-plane mismatch = %.3e\n', ...
    repG1_after.maxPositionError, repG1_after.maxCrossPlaneMismatch);

%% 4) Plot a sampled view
figure('Name','Patch tangent continuity demo'); hold on; axis equal; grid on;
title('Upper/Lower patches after symmetry and G1 enforcement');

plotSurfaceWire(patchA.S, 21, 21);
plotSurfaceWire(patchB.S, 21, 21);

% Plot shared boundary and symmetry edge
plotEdge(patchA.S, 'v0', 'k-', 2.0);
plotEdge(patchA.S, 'v1', 'r-', 2.0);
plotEdge(patchB.S, 'v0', 'b--', 1.5);
view(3);

fprintf('Done.\n');

function plotSurfaceWire(S, nu, nv)
    u = linspace(S.domainU(1), S.domainU(2), nu);
    v = linspace(S.domainV(1), S.domainV(2), nv);
    [U,V] = meshgrid(u,v);
    P = zeros(numel(U),3);
    for ii = 1:numel(U)
        P(ii,:) = S.evaluate(U(ii), V(ii));
    end
    X = reshape(P(:,1), size(U));
    Y = reshape(P(:,2), size(U));
    Z = reshape(P(:,3), size(U));
    mesh(X,Y,Z,'FaceColor','none');
end

function plotEdge(S, edgeId, style, lw)
    uv = geom.Patch.edgeSampleUV(S, edgeId, 81);
    P = zeros(size(uv,1),3);
    for ii = 1:size(uv,1)
        P(ii,:) = S.evaluate(uv(ii,1), uv(ii,2));
    end
    plot3(P(:,1), P(:,2), P(:,3), style, 'LineWidth', lw);
end
