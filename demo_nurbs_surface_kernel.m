%% DEMO_NURBS_SURFACE_KERNEL.m
% Demonstration / validation script for the NURBS surface kernel,
% including construction algorithms:
%   - loft / skinning
%   - Coons patch
%   - Gordon surface
%   - surface of revolution
%   - sweep surface

clear; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir)
    addpath(script_dir);
end

if ~exist('+geom', 'dir') && isempty(which('geom.NURBSSurface'))
    error(['Cannot find +geom package.\n' ...
           'Expected it at: %s\n' ...
           'Add the geom kernel folder to your MATLAB path:\n' ...
           '  addpath(''%s'')'], ...
           fullfile(script_dir, '+geom'), script_dir);
end
fprintf('geom package found at: %s\n', script_dir);

%% ======================================================================
%  1. BASE SURFACE
%% ======================================================================
fprintf('\n--- 1. Base surface kernel checks ---\n');

P = zeros(4,4,3);
xu = [0.0 1.0 2.2 3.2];
yv = [0.0 0.8 1.8 3.0];

for i = 1:4
    for j = 1:4
        x = xu(i);
        y = yv(j);
        z = 0.28 * sin(0.9*x) * cos(0.8*y) + 0.08 * x - 0.05 * y;
        P(i,j,:) = [x, y, z];
    end
end

S = geom.NURBSSurface(P, 3, 3);
fprintf('  validate() = %d\n', S.validate());

u0 = 0.37;
v0 = 0.61;
Spt = S.evaluate(u0, v0);
[Su, Sv] = S.partialDerivatives(u0, v0);
[K, H] = S.curvatures(u0, v0);

fprintf('  S(%.3f, %.3f) = [%.6f %.6f %.6f]\n', u0, v0, Spt);
fprintf('  |Su| = %.6f, |Sv| = %.6f\n', norm(Su), norm(Sv));
fprintf('  K = %.6e, H = %.6e\n', K, H);

fig1 = figure('Name','1 - Base surface');
S.plot(24, 24, 'ShowCP', true, 'ShowIso', true);
title('Base surface');
view(3);

%% ======================================================================
%  2. LOFT / SKINNING
%% ======================================================================
fprintf('\n--- 2. Loft / skinning ---\n');

try
    c1 = geom.NURBSCurve([0 0 0; 1 0.4 0.1; 2 0.1 0.0; 3 0 0], 3);
    c2 = geom.NURBSCurve([0 0.5 0.6; 1 0.9 0.8; 2 0.8 0.5; 3 0.6 0.3], 3);
    c3 = geom.NURBSCurve([0 1.1 1.0; 1 1.4 1.2; 2 1.5 0.7; 3 1.4 0.4], 3);
    c4 = geom.NURBSCurve([0 1.8 1.2; 1 2.0 1.0; 2 2.2 0.5; 3 2.3 0.2], 3);

    Sloft = geom.NURBSSurface.loft({c1,c2,c3,c4}, 3, 'centripetal');
    fprintf('  loft net size = %d x %d\n', size(Sloft.P,1), size(Sloft.P,2));

    fig2 = figure('Name','2 - Lofted surface');
    Sloft.plot(24, 24, 'ShowCP', true, 'ShowIso', true);
    hold on;
    c1.plot(120, 'ShowCP', false, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
    c2.plot(120, 'ShowCP', false, 'Color', [0.1 0.5 0.9], 'LineWidth', 2.0);
    c3.plot(120, 'ShowCP', false, 'Color', [0.1 0.7 0.2], 'LineWidth', 2.0);
    c4.plot(120, 'ShowCP', false, 'Color', [0.8 0.4 0.1], 'LineWidth', 2.0);
    title('Loft / skinning');
    view(3);
catch ME
    fprintf('  Loft section failed: %s\n', ME.message);
end

%% ======================================================================
%  3. COONS PATCH
%% ======================================================================
fprintf('\n--- 3. Coons patch ---\n');

try
    Cu0 = geom.NURBSCurve([0 0 0; 1 0.0 0.3; 2 0.0 0.2; 3 0 0], 3);
    Cu1 = geom.NURBSCurve([0 2 0.2; 1 2.0 0.8; 2 2.0 0.7; 3 2 0.1], 3);
    Cv0 = geom.NURBSCurve([0 0 0; 0 0.8 0.4; 0 1.4 0.5; 0 2 0.2], 3);
    Cv1 = geom.NURBSCurve([3 0 0; 3 0.7 0.2; 3 1.4 0.4; 3 2 0.1], 3);

    Scoons = geom.NURBSSurface.coons(Cu0, Cu1, Cv0, Cv1, 3, 3, 21, 21);
    fprintf('  Coons net size = %d x %d\n', size(Scoons.P,1), size(Scoons.P,2));

    fig3 = figure('Name','3 - Coons patch');
    Scoons.plot(24, 24, 'ShowCP', true, 'ShowIso', true);
    hold on;
    Cu0.plot(100, 'ShowCP', false, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
    Cu1.plot(100, 'ShowCP', false, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
    Cv0.plot(100, 'ShowCP', false, 'Color', [0.1 0.1 0.9], 'LineWidth', 2.0);
    Cv1.plot(100, 'ShowCP', false, 'Color', [0.1 0.1 0.9], 'LineWidth', 2.0);
    title('Coons patch');
    view(3);
catch ME
    fprintf('  Coons section failed: %s\n', ME.message);
end

%% ======================================================================
%  4. GORDON SURFACE
%% ======================================================================
fprintf('\n--- 4. Gordon surface ---\n');

try
    p1 = geom.NURBSCurve([0 0 0; 1 0.0 0.4; 2 0.0 0.3; 3 0 0], 3);
    p2 = geom.NURBSCurve([0 0.8 0.3; 1 0.8 0.8; 2 0.8 0.6; 3 0.8 0.2], 3);
    p3 = geom.NURBSCurve([0 1.6 0.1; 1 1.6 0.7; 2 1.6 0.9; 3 1.6 0.3], 3);
    p4 = geom.NURBSCurve([0 2.4 0.0; 1 2.4 0.2; 2 2.4 0.5; 3 2.4 0.1], 3);

    g1 = geom.NURBSCurve([0 0 0; 0 0.8 0.4; 0 1.6 0.1; 0 2.4 0], 3);
    g2 = geom.NURBSCurve([1 0.0 0.4; 1 0.8 0.8; 1 1.6 0.7; 1 2.4 0.2], 3);
    g3 = geom.NURBSCurve([2 0.0 0.3; 2 0.8 0.6; 2 1.6 0.9; 2 2.4 0.5], 3);
    g4 = geom.NURBSCurve([3 0.0 0; 3 0.8 0.2; 3 1.6 0.3; 3 2.4 0.1], 3);

    Sgordon = geom.NURBSSurface.gordon({p1,p2,p3,p4}, {g1,g2,g3,g4}, 3, 3, 25, 25);
    fprintf('  Gordon net size = %d x %d\n', size(Sgordon.P,1), size(Sgordon.P,2));

    fig4 = figure('Name','4 - Gordon surface');
    Sgordon.plot(28, 28, 'ShowCP', true, 'ShowIso', true);
    hold on;
    p1.plot(100, 'ShowCP', false, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
    p2.plot(100, 'ShowCP', false, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
    p3.plot(100, 'ShowCP', false, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
    p4.plot(100, 'ShowCP', false, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
    g1.plot(100, 'ShowCP', false, 'Color', [0.1 0.1 0.9], 'LineWidth', 2.0);
    g2.plot(100, 'ShowCP', false, 'Color', [0.1 0.1 0.9], 'LineWidth', 2.0);
    g3.plot(100, 'ShowCP', false, 'Color', [0.1 0.1 0.9], 'LineWidth', 2.0);
    g4.plot(100, 'ShowCP', false, 'Color', [0.1 0.1 0.9], 'LineWidth', 2.0);
    title('Gordon surface');
    view(3);
catch ME
    fprintf('  Gordon section failed: %s\n', ME.message);
end

%% ======================================================================
%  5. SURFACE OF REVOLUTION
%% ======================================================================
fprintf('\n--- 5. Surface of revolution ---\n');

try
    Cprof = geom.NURBSCurve([1.0 0.0 0.0;
                             1.2 0.0 0.4;
                             1.4 0.0 0.8;
                             1.6 0.0 1.2], 3);

    Srev = geom.NURBSSurface.revolve(Cprof, [0 0 0], [0 0 1], 2*pi, 3, 9);
    fprintf('  revolution net size = %d x %d\n', size(Srev.P,1), size(Srev.P,2));

    fig5 = figure('Name','5 - Surface of revolution');
    Srev.plot(30, 28, 'ShowCP', true, 'ShowIso', true);
    title('Surface of revolution');
    view(3);
catch ME
    fprintf('  Revolution section failed: %s\n', ME.message);
end

%% ======================================================================
%  6. SWEPT SURFACE
%% ======================================================================
fprintf('\n--- 6. Swept surface ---\n');

try
    % Profile curve in local coordinates
    Cprof2 = geom.NURBSCurve([0  -0.4  0.0;
                              0  -0.1  0.2;
                              0   0.2  0.2;
                              0   0.5  0.0], 3);

    % Spine curve
    Cspine = geom.NURBSCurve([0.0 0.0 0.0;
                              1.0 0.4 0.5;
                              2.2 0.2 1.1;
                              3.2 0.8 1.6], 3);

    Sswp = geom.NURBSSurface.sweep(Cprof2, Cspine, 3, 10, [0 0 1]);
    fprintf('  sweep net size = %d x %d\n', size(Sswp.P,1), size(Sswp.P,2));

    fig6 = figure('Name','6 - Swept surface');
    Sswp.plot(28, 22, 'ShowCP', true, 'ShowIso', true);
    hold on;
    Cspine.plot(200, 'ShowCP', false, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
    title('Swept surface');
    view(3);
catch ME
    fprintf('  Sweep section failed: %s\n', ME.message);
end

%% ======================================================================
%  SUMMARY
%% ======================================================================
fprintf('\n=== Surface Constructor Demo Complete ===\n');
fprintf('Figures:\n');
fprintf('  (1) Base surface\n');
fprintf('  (2) Loft / skinning\n');
fprintf('  (3) Coons patch\n');
fprintf('  (4) Gordon surface\n');
fprintf('  (5) Surface of revolution\n');
fprintf('  (6) Swept surface\n');