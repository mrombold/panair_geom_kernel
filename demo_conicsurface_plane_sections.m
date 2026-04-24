
function demo_conicsurface_plane_sections()
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
fprintf('\n=== ConicSurface plane-section demo ===\n');

[UpperGuide,~]    = geom.Loft.limingConic([0 0 40], [100 0 55], [25 0 55], [25 0 50]);
[MaxBSide,~]      = geom.Loft.limingConic([0 0 35], [100 0 40], [25 0 40], [25 0 36.5]);
[MaxBTop,~]       = geom.Loft.limingConic([0 5 0], [100 10 0], [25 10 0], [25 8 0]);
[UprShldrSide,~]  = geom.Loft.limingConic([0 0 38], [100 0 52], [25 0 52], [25 0 47]);
[UprShldrTop,~]   = geom.Loft.limingConic([0 4 0], [100 8 0], [25 8 0], [25 7 0]);

[LowerGuide,~]    = geom.Loft.combinePlanarGuidesTo3D(MaxBTop, MaxBSide);
[ShoulderGuide,~] = geom.Loft.combinePlanarGuidesTo3D(UprShldrTop, UprShldrSide);
[TangencyGuide,~] = geom.Loft.combinePlanarGuidesTo3D(MaxBTop, UpperGuide);

S = geom.ConicSurface( ...
    'UpperGuide', UpperGuide, ...
    'LowerGuide', LowerGuide, ...
    'TangencyGuide', TangencyGuide, ...
    'ShoulderGuide', ShoulderGuide, ...
    'SweepOrigin', [0 0 0], ...
    'SweepVector', [1 0 0], ...
    'StationRange', [0 100], ...
    'EnableSectionCache', true);

M = S.isoMesh(101, 41);

zWater = 42.0;
yButt  = 4.0;
x0     = 50.0;
planePointCant = [x0 4.0 44.0];
nCant  = [1.0 0.12 -0.04];

[Cwater, infoWater] = geom.conicSurfaceSection(S, 'z', zWater, ...
    'SamplesU', 121, ...
    'RefinePasses', 1, ...
    'ExtendToGuides', true, ...
    'MultiBranchAction', 'warn_primary', ...
    'BranchCountAction', 'split');

[Cbutt, infoButt] = geom.conicSurfaceSection(S, 'y', yButt, ...
    'SamplesU', 121, ...
    'RefinePasses', 1, ...
    'ExtendToGuides', true, ...
    'MultiBranchAction', 'warn_primary', ...
    'BranchCountAction', 'split');

[Ccant, infoCant] = geom.conicSurfaceSection(S, planePointCant, nCant, ...
    'SamplesU', 241, ...
    'RefinePasses', 2, ...
    'ExtendToGuides', true, ...
    'MultiBranchAction', 'warn_primary', ...
    'BranchCountAction', 'split');

fprintf('Waterline: branches = %d, fit RMS = %.3e, guide attachments = %d\n', ...
    infoWater.numBranches, infoWater.fitRMS, infoWater.guideAttachment.numAttached);
fprintf('Buttline : branches = %d, fit RMS = %.3e, guide attachments = %d\n', ...
    infoButt.numBranches, infoButt.fitRMS, infoButt.guideAttachment.numAttached);
fprintf('Canted   : branches = %d, fit RMS = %.3e, guide attachments = %d\n', ...
    infoCant.numBranches, infoCant.fitRMS, infoCant.guideAttachment.numAttached);

figure('Name', 'ConicSurface plane sections', 'Color', 'w');
hold on;

UpperGuide.plot();
LowerGuide.plot();
TangencyGuide.plot();
ShoulderGuide.plot();

surf(M.X, M.Y, M.Z, ...
    'FaceAlpha', 0.75, ...
    'EdgeColor', [0.25 0.25 0.25], ...
    'FaceColor', [0.60 0.78 0.96]);

Cwater.plot(250, 'Color', [0.05 0.20 0.95], 'LineWidth', 2.5);
Cbutt.plot(250,  'Color', [0.00 0.60 0.15], 'LineWidth', 2.5);
Ccant.plot(250,  'Color', [0.90 0.35 0.05], 'LineWidth', 2.5);

localDrawPlanePatch(M, [0 0 zWater], [0 0 1], [0.75 0.85 1.00], 0.18);
localDrawPlanePatch(M, [0 yButt 0], [0 1 0], [0.82 1.00 0.82], 0.18);
localDrawPlanePatch(M, planePointCant, nCant, [1.00 0.88 0.78], 0.18);

axis vis3d;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Guide-driven ConicSurface with guide-extended section cuts');
view(3);
%grid on;
ax = gca;
ax.Clipping = 'off';

legend({'UpperGuide','LowerGuide','TangencyGuide','ShoulderGuide', ...
        'ConicSurface', ...
        sprintf('Waterline  Z = %.1f', zWater), ...
        sprintf('Buttline   Y = %.1f', yButt), ...
        sprintf('Canted station plane through x = %.1f', x0)}, ...
        'Location', 'bestoutside');

fprintf('\nDone.\n');

end

function localDrawPlanePatch(M, p0, n, faceColor, alphaVal)
n = n(:);
n = n / norm(n);

X = M.X(:); Y = M.Y(:); Z = M.Z(:);
mins = [min(X) min(Y) min(Z)];
maxs = [max(X) max(Y) max(Z)];
span = norm(maxs - mins);

tmp = [1;0;0];
if abs(dot(tmp,n)) > 0.9
    tmp = [0;1;0];
end
e1 = tmp - dot(tmp,n) * n;
e1 = e1 / norm(e1);
e2 = cross(n, e1);

p0 = p0(:);
s = 0.35 * span;
corn = [ -1 -1
          1 -1
          1  1
         -1  1 ] * s;

P = zeros(4,3);
for k = 1:4
    P(k,:) = (p0 + corn(k,1)*e1 + corn(k,2)*e2).';
end

patch('XData', P(:,1), 'YData', P(:,2), 'ZData', P(:,3), ...
      'FaceColor', faceColor, 'FaceAlpha', alphaVal, ...
      'EdgeColor', 'none');
end
