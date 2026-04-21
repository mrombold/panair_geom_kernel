%% demo_patch_tangent_continuity_clean
% Cleaner continuity demo for geom.Patch:
%   1) symmetry-plane tangency on a single strip patch
%   2) G1 tangency between two adjacent strip patches
%
% This demo uses simple rectangular interpolation nets so the effect of the
% continuity enforcement is easy to see.

clear; close all; clc;

if ~exist('+geom', 'dir') && isempty(which('geom.NURBSSurface'))
    error('Run this from the repo root or add the geometry kernel to the MATLAB path.');
end
if isempty(which('geom.Patch'))
    error('geom.Patch not found. Add Patch.m to +geom first.');
end

fprintf('=== Clean Patch Tangent Continuity Demo ===\n');

%% ------------------------------------------------------------------------
% 1) Symmetry-plane tangency on a single crown patch
% -------------------------------------------------------------------------
fprintf('\n--- 1) Symmetry plane tangency ---\n');

xStations = [0 4 8 12];
yRowsSym  = [0.0 0.35 0.75 1.15];

Qsym = zeros(numel(yRowsSym), numel(xStations), 3);
for j = 1:numel(yRowsSym)
    y = yRowsSym(j);
    for i = 1:numel(xStations)
        x = xStations(i);
        % Mild crown shape in z. Along y=0 this should become symmetry edge.
        z = 1.15 - 0.22*y^2 + 0.10*cos(pi*x/xStations(end));
        Qsym(j,i,:) = [x, y, z];
    end
end

Ssym = geom.NURBSSurface.globalInterpNet(Qsym, 3, 3, 'centripetal', 'centripetal');
patchSym = geom.Patch('SymmetryStrip', Ssym);

repSymBefore = patchSym.symmetryReport('u0', [0 0 0], [0 1 0]);
repSymAfter  = patchSym.enforceSymmetryPlane('u0', [0 0 0], [0 1 0]);

fprintf('Before symmetry: max edge-plane distance     = %.3e\n', repSymBefore.maxEdgePlaneDistance);
fprintf('After  symmetry: max edge-plane distance     = %.3e\n', repSymAfter.maxEdgePlaneDistance);
fprintf('After  symmetry: max cross-derivative error  = %.3e\n', repSymAfter.maxCrossDerivativeOffNormal);

%% ------------------------------------------------------------------------
% 2) G1 tangency between two adjacent patches
% -------------------------------------------------------------------------
fprintf('\n--- 2) G1 continuity between adjacent patches ---\n');

yRowsA = [0.0 0.45 0.90 1.30];
yRowsB = [1.30 1.80 2.25 2.70];
sharedY = yRowsA(end);

QA = zeros(numel(yRowsA), numel(xStations), 3);
QB = zeros(numel(yRowsB), numel(xStations), 3);

for j = 1:numel(yRowsA)
    y = yRowsA(j);
    for i = 1:numel(xStations)
        x = xStations(i);
        z = 1.00 - 0.18*(y-0.25)^2 + 0.08*cos(pi*x/xStations(end));
        QA(j,i,:) = [x, y, z];
    end
end

for j = 1:numel(yRowsB)
    y = yRowsB(j);
    for i = 1:numel(xStations)
        x = xStations(i);

        if abs(y - sharedY) < 1e-12
            % Match the shared boundary exactly.
            z = 1.00 - 0.18*(y-0.25)^2 + 0.08*cos(pi*x/xStations(end));
        else
            % Deliberately choose a different shoulder-side shape so the
            % initial join is only G0, not G1.
            dy = y - sharedY;
            z = 0.92 - 0.08*dy - 0.05*dy^2 + 0.05*cos(pi*x/xStations(end) + 0.25);
        end

        QB(j,i,:) = [x, y, z];
    end
end

SA = geom.NURBSSurface.globalInterpNet(QA, 3, 3, 'centripetal', 'centripetal');
SB = geom.NURBSSurface.globalInterpNet(QB, 3, 3, 'centripetal', 'centripetal');

patchA = geom.Patch('PatchA', SA);
patchB = geom.Patch('PatchB', SB);

repG1Before = patchB.g1ReportWithPatch('u0', patchA, 'u1');
repG1After  = patchB.enforceG1WithPatch('u0', patchA, 'u1', ...
    'Lambda', 1.0, 'AverageEdge', true);

fprintf('Before G1: max position err          = %.3e\n', repG1Before.maxPositionError);
fprintf('Before G1: max cross-plane mismatch  = %.3e\n', repG1Before.maxCrossPlaneMismatch);
fprintf('After  G1: max position err          = %.3e\n', repG1After.maxPositionError);
fprintf('After  G1: max cross-plane mismatch  = %.3e\n', repG1After.maxCrossPlaneMismatch);

%% ------------------------------------------------------------------------
% 3) Plotting
% -------------------------------------------------------------------------
figure('Name','Clean continuity demo','Color','w');

tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% Symmetry plot
nexttile; hold on; axis equal; grid on; view(3);
title('Symmetry-plane tangency');
plotSurfaceWire(Ssym, 25, 25, [0.75 0.75 0.75]);
plotSurfaceWire(patchSym.S, 25, 25, [0.00 0.45 0.85]);
plotEdge(patchSym.S, 'u0', 'k-', 2.0);

% Draw y = 0 symmetry plane as a light rectangle
xp = [xStations(1) xStations(end) xStations(end) xStations(1)];
yp = [0 0 0 0];
zp = [0.3 0.3 1.4 1.4];
patch(xp, yp, zp, [0.9 0.9 0.9], 'FaceAlpha', 0.25, 'EdgeColor', 'none');

xlabel('x'); ylabel('y'); zlabel('z');
legend({'before','after','symmetry edge','y=0 plane'}, 'Location','best');

% Patch-to-patch G1 plot
nexttile; hold on; axis equal; grid on; view(3);
title('Adjacent patches: G1 enforcement');
plotSurfaceWire(SA, 25, 25, [0.10 0.55 0.10]);
plotSurfaceWire(SB, 25, 25, [0.70 0.25 0.25]);
plotSurfaceWire(patchB.S, 25, 25, [0.00 0.45 0.85]);

plotEdge(SA,       'u1', 'k-',  2.0);
plotEdge(SB,       'u0', 'r--', 1.5);
plotEdge(patchB.S, 'u0', 'b-',  2.0);

xlabel('x'); ylabel('y'); zlabel('z');
legend({'Patch A','Patch B before','Patch B after','shared edge A','shared edge B before','shared edge B after'}, ...
    'Location','best');

fprintf('\nDone.\n');

%% ------------------------------------------------------------------------
function plotSurfaceWire(S, nu, nv, colorRGB)
    u = linspace(S.domainU(1), S.domainU(2), nu);
    v = linspace(S.domainV(1), S.domainV(2), nv);

    % Constant-v lines
    for j = 1:numel(v)
        P = zeros(numel(u),3);
        for i = 1:numel(u)
            P(i,:) = S.evaluate(u(i), v(j));
        end
        plot3(P(:,1), P(:,2), P(:,3), '-', 'Color', colorRGB, 'LineWidth', 0.75);
    end

    % Constant-u lines
    for i = 1:numel(u)
        P = zeros(numel(v),3);
        for j = 1:numel(v)
            P(j,:) = S.evaluate(u(i), v(j));
        end
        plot3(P(:,1), P(:,2), P(:,3), '-', 'Color', colorRGB, 'LineWidth', 0.75);
    end
end

function plotEdge(S, edgeId, style, lw)
    uv = geom.Patch.edgeSampleUV(S, edgeId, 101);
    P = zeros(size(uv,1),3);
    for ii = 1:size(uv,1)
        P(ii,:) = S.evaluate(uv(ii,1), uv(ii,2));
    end
    plot3(P(:,1), P(:,2), P(:,3), style, 'LineWidth', lw);
end
