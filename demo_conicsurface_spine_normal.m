%% demo_conicsurface_spine_normal.m
clear; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir)
    addpath(script_dir);
end
if ~exist('+geom', 'dir') && isempty(which('geom.NURBSCurve'))
    error(['Cannot find +geom package.\n' ...
           'Run this from the repo root or add the geometry kernel to path.']);
end
fprintf('geom package found at: %s\n', script_dir);
fprintf('\n=== ConicSurface spine-normal demo ===\n');

Psp = [ ...
     0.0   0.0   0.0;
    20.0   5.0   2.0;
    45.0  12.0   6.0;
    75.0  18.0  11.0;
   100.0  25.0  15.0];
Spine = geom.NURBSCurve.globalInterp(Psp, 3, 'centripetal');

Pup = [ ...
     0.0   0.0  18.0;
    20.0   5.0  19.0;
    45.0  12.0  21.0;
    75.0  18.0  24.0;
   100.0  25.0  26.0];
Plo = [ ...
     0.0   8.0   7.0;
    20.0  12.0   7.5;
    45.0  19.0   8.5;
    75.0  25.0   9.5;
   100.0  31.0  10.5];
Ptan = [ ...
     0.0   8.0  18.0;
    20.0  12.0  19.0;
    45.0  19.0  20.5;
    75.0  25.0  22.0;
   100.0  31.0  23.0];
Psh = [ ...
     0.0   6.0  15.0;
    20.0   9.5  16.0;
    45.0  15.5  17.8;
    75.0  21.0  19.8;
   100.0  26.5  21.0];

UpperGuide    = geom.NURBSCurve.globalInterp(Pup, 3, 'centripetal');
LowerGuide    = geom.NURBSCurve.globalInterp(Plo, 3, 'centripetal');
TangencyGuide = geom.NURBSCurve.globalInterp(Ptan, 3, 'centripetal');
ShoulderGuide = geom.NURBSCurve.globalInterp(Psh, 3, 'centripetal');

S = geom.ConicSurface.fromSpineNormal( ...
    'Spine', Spine, ...
    'SpineDomain', Spine.domain, ...
    'UpperGuide', UpperGuide, ...
    'LowerGuide', LowerGuide, ...
    'TangencyGuide', TangencyGuide, ...
    'ShoulderGuide', ShoulderGuide, ...
    'PlaneXHint', [0 0 1], ...
    'EnableSectionCache', true);

M = S.isoMesh(41, 31);

figure('Name','ConicSurface spine-normal','Color','w'); hold on;
surf(M.X, M.Y, M.Z, ...
    'FaceAlpha', 0.8, ...
    'EdgeColor', [0.25 0.25 0.25], ...
    'FaceColor', [0.60 0.78 0.96]);

plotSampled(Spine,         'k-', 'LineWidth', 2.0);
plotSampled(UpperGuide,    'b-');
plotSampled(LowerGuide,    'r-');
plotSampled(TangencyGuide, '-', 'Color', [0.95 0.65 0.10]);
plotSampled(ShoulderGuide, 'm-');

svals = linspace(Spine.domain(1), Spine.domain(2), 6);
for k = 1:numel(svals)
    [Ck, meta] = S.sectionAtStation(svals(k));
    tk = linspace(Ck.domain(1), Ck.domain(2), 100);
    Pk = Ck.evaluate(tk);
    plot3(Pk(:,1), Pk(:,2), Pk(:,3), 'k-', 'LineWidth', 1.0);

    O = meta.guideInfo.frame.origin;
    n = meta.guideInfo.frame.normal;
    q = O + 5*n;
    plot3([O(1) q(1)], [O(2) q(2)], [O(3) q(3)], 'g-');
end

axis equal; grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Guide-driven ConicSurface with spine-normal station planes');
view(3);

function h = plotSampled(C, style, varargin)
    t = linspace(C.domain(1), C.domain(2), 200);
    P = C.evaluate(t);
    h = plot3(P(:,1), P(:,2), P(:,3), style, varargin{:});
end
