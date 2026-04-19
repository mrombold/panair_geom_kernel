
%% demo_loft_liming_fuselage_sections.m
clear; close all; clc;

%% Path / package check
script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir)
    addpath(script_dir);
end
if ~exist('+geom', 'dir') && isempty(which('geom.NURBSCurve'))
    error(['Cannot find +geom package.\n' ...
           'Run this from the repo root or add the geometry kernel to path.']);
end
fprintf('geom package found at: %s\n', script_dir);

%% User settings
pGuide = 3;
fitMethod = 'centripetal';
loftDegree = 3;
shoulderParameter = [];   % [] = solve shoulder parameter from geometry
stationX = [0.40 1.00 2.00 4.00 6.50 8.50 10.00 11.20];

% Editable guide points in XYZ
P_top = [
    0.00  0.00  0.00
    0.40  0.00  0.33
    1.00  0.00  0.75
    2.00  0.00  1.03
    4.00  0.00  1.10
    6.50  0.00  1.10
    8.50  0.00  1.08
    10.00 0.00  0.90
    11.20 0.00  0.46
    12.00 0.00  0.00];

P_bot = [
    0.00  0.00  0.00
    0.40  0.00 -0.27
    1.00  0.00 -0.53
    2.00  0.00 -0.80
    4.00  0.00 -1.10
    6.50  0.00 -1.10
    8.50  0.00 -1.00
    10.00 0.00 -0.72
    11.20 0.00 -0.30
    12.00 0.00  0.00];

P_maxB = [
    0.00  0.00  0.00
    0.40  0.24  0.02
    1.00  0.64  0.10
    2.00  0.99  0.16
    4.00  1.10  0.00
    6.50  1.10  0.00
    8.50  1.08  0.04
    10.00 0.84  0.02
    11.20 0.37  0.00
    12.00 0.00  0.00];

P_shU = [
    0.00  0.00  0.00
    0.40  0.06  0.16
    1.00  0.16  0.42
    2.00  0.34  0.72
    4.00  0.52  0.98
    6.50  0.52  0.98
    8.50  0.52  0.94
    10.00 0.38  0.74
    11.20 0.15  0.28
    12.00 0.00  0.00];

P_shL = [
    0.00  0.00  0.00
    0.40  0.08 -0.10
    1.00  0.20 -0.26
    2.00  0.42 -0.56
    4.00  0.60 -0.94
    6.50  0.60 -0.94
    8.50  0.58 -0.82
    10.00 0.42 -0.58
    11.20 0.16 -0.20
    12.00 0.00  0.00];

%% 1) Build longitudinal guide curves
fprintf('\n--- 1. Build longitudinal guide curves ---\n');
C_top  = geom.NURBSCurve.globalInterp(P_top,  pGuide, fitMethod);
C_bot  = geom.NURBSCurve.globalInterp(P_bot,  pGuide, fitMethod);
C_maxB = geom.NURBSCurve.globalInterp(P_maxB, pGuide, fitMethod);
C_shU  = geom.NURBSCurve.globalInterp(P_shU,  pGuide, fitMethod);
C_shL  = geom.NURBSCurve.globalInterp(P_shL,  pGuide, fitMethod);
fprintf(' built top / bottom / max-breadth / shoulder guide curves\n');

%% 2) Build Liming-style fuselage sections at stations
fprintf('\n--- 2. Build fuselage sections from guide curves ---\n');
sections = geom.Loft.buildSectionFamily(stationX, ...
    'UpperProfile', C_top, ...
    'LowerProfile', C_bot, ...
    'MaxBreadth', C_maxB, ...
    'UpperShoulder', C_shU, ...
    'LowerShoulder', C_shL, ...
    'Axis', 'x', ...
    'ShoulderParameter', shoulderParameter);

Cup = cell(size(sections));
Clo = cell(size(sections));
for k = 1:numel(sections)
    Cup{k} = sections{k}.upperCurve;
    Clo{k} = sections{k}.lowerCurve;
    ptsk = sections{k}.points;
    fprintf([' section %2d: x = %7.3f | top z = %7.4f | max y = %7.4f | ' ...
             'bot z = %7.4f\n'], ...
        k, sections{k}.stationValue, ptsk.top(3), ptsk.maxBreadth(2), ptsk.bottom(3));
    fprintf('            upper: u = %.4f, w = %.4f | lower: u = %.4f, w = %.4f\n', ...
        sections{k}.upperMeta.parameter, sections{k}.upperMeta.weight, ...
        sections{k}.lowerMeta.parameter, sections{k}.lowerMeta.weight);
end

%% 3) Loft upper / lower half-fuselage surfaces
fprintf('\n--- 3. Loft upper / lower surfaces ---\n');
S_upper = geom.Loft.loftSections(Cup, loftDegree, 'centripetal');
S_lower = geom.Loft.loftSections(Clo, loftDegree, 'centripetal');
fprintf(' lofted upper and lower half-fuselage surfaces\n');

%% 4) Plot guide curves, sections, and surfaces
figure('Name','Loft demo: Liming fuselage sections');
hold on; grid on; axis equal; view(3);

C_top.plot([],  'Color', [0.85 0.15 0.15], 'LineWidth', 2.0, 'ShowCP', false);
C_bot.plot([],  'Color', [0.15 0.25 0.85], 'LineWidth', 2.0, 'ShowCP', false);
C_maxB.plot([], 'Color', [0.10 0.60 0.10], 'LineWidth', 2.0, 'ShowCP', false);
C_shU.plot([],  'Color', [0.90 0.55 0.10], 'LineWidth', 1.5, 'ShowCP', false);
C_shL.plot([],  'Color', [0.55 0.25 0.85], 'LineWidth', 1.5, 'ShowCP', false);

plot3(P_top(:,1),  P_top(:,2),  P_top(:,3),  'o', 'MarkerSize', 5, 'LineWidth', 1.0);
plot3(P_bot(:,1),  P_bot(:,2),  P_bot(:,3),  'o', 'MarkerSize', 5, 'LineWidth', 1.0);
plot3(P_maxB(:,1), P_maxB(:,2), P_maxB(:,3), 'o', 'MarkerSize', 5, 'LineWidth', 1.0);
plot3(P_shU(:,1),  P_shU(:,2),  P_shU(:,3),  'o', 'MarkerSize', 5, 'LineWidth', 1.0);
plot3(P_shL(:,1),  P_shL(:,2),  P_shL(:,3),  'o', 'MarkerSize', 5, 'LineWidth', 1.0);

for k = 1:numel(Cup)
    Cup{k}.plot(120, 'Color', [0.75 0.35 0.10], 'LineWidth', 1.5, 'ShowCP', false);
    Clo{k}.plot(120, 'Color', [0.10 0.35 0.75], 'LineWidth', 1.5, 'ShowCP', false);

    p = sections{k}.points;
    plot3(p.top(1),          p.top(2),          p.top(3),          'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
    plot3(p.maxBreadth(1),   p.maxBreadth(2),   p.maxBreadth(3),   'kd', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
    plot3(p.bottom(1),       p.bottom(2),       p.bottom(3),       'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
    plot3(p.upperShoulder(1),p.upperShoulder(2),p.upperShoulder(3),'k^', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
    plot3(p.lowerShoulder(1),p.lowerShoulder(2),p.lowerShoulder(3),'kv', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
end

nuPlot = 41;
nvPlot = 31;
meshU = S_upper.isoMesh(nuPlot, nvPlot, 'SpacingU', 'cosine', 'SpacingV', 'linear');
meshL = S_lower.isoMesh(nuPlot, nvPlot, 'SpacingU', 'cosine', 'SpacingV', 'linear');

surf(meshU.X, meshU.Y, meshU.Z, ...
    'FaceColor', [0.95 0.70 0.55], 'FaceAlpha', 0.55, ...
    'EdgeColor', 'k', 'EdgeAlpha', 0.15);

surf(meshL.X, meshL.Y, meshL.Z, ...
    'FaceColor', [0.60 0.70 0.95], 'FaceAlpha', 0.55, ...
    'EdgeColor', 'k', 'EdgeAlpha', 0.15);

plot3([0 12], [0 0], [0 0], 'k--', 'LineWidth', 1.0);

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Liming-style fuselage sections and lofted half-body surfaces');

%% 5) Plot one station in local section coordinates
kShow = min(4, numel(sections));
sec = sections{kShow};

uPlot = linspace(sec.upperCurve.domain(1), sec.upperCurve.domain(2), 150);
vPlot = linspace(sec.lowerCurve.domain(1), sec.lowerCurve.domain(2), 150);
Qup2 = geom.Loft.toPlane2D(sec.upperCurve.evaluate(uPlot), sec.frame);
Qlo2 = geom.Loft.toPlane2D(sec.lowerCurve.evaluate(vPlot), sec.frame);

figure('Name','Loft demo: section plane view');
hold on; grid on; axis equal;
plot(Qup2(:,1), Qup2(:,2), 'LineWidth', 2.0);
plot(Qlo2(:,1), Qlo2(:,2), 'LineWidth', 2.0);

p2 = sec.points2D;
plot(p2.top(1),           p2.top(2),           'ks', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(p2.maxBreadth(1),    p2.maxBreadth(2),    'kd', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(p2.bottom(1),        p2.bottom(2),        'ks', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(p2.upperShoulder(1), p2.upperShoulder(2), 'k^', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(p2.lowerShoulder(1), p2.lowerShoulder(2), 'kv', 'MarkerSize', 7, 'MarkerFaceColor', 'k');

tUp2 = geom.Loft.toPlane2D(sec.tangentIntersectionUpper, sec.frame);
tLo2 = geom.Loft.toPlane2D(sec.tangentIntersectionLower, sec.frame);
plot(tUp2(1), tUp2(2), 'ro', 'MarkerSize', 7, 'LineWidth', 1.5);
plot(tLo2(1), tLo2(2), 'bo', 'MarkerSize', 7, 'LineWidth', 1.5);

xlabel('Local breadth coordinate');
ylabel('Local vertical coordinate');
title(sprintf('Section plane view at station x = %.3f', sec.stationValue));

fprintf('\nDone.\n');
fprintf('Try editing the guide-point arrays near the top of this script.\n');
