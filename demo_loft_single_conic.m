%% demo_loft_single_conic.m
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

% -------------------------------------------------------------------------
% User-editable Liming conic construction points
% -------------------------------------------------------------------------
P0 = [0.0, 0.0, 1.0];
P1 = [1.2, 0.0, 0.0];
T  = [0.0, 0.0, 0.0];
S  = [0.42, 0.0, 0.38];

fprintf('\n--- Build single Liming conic ---\n');

[C, meta] = geom.Loft.limingConic(P0, P1, T, S);

if isstruct(meta)
    if isfield(meta, 'parameter')
        fprintf(' solved/estimated shoulder parameter u = %.8f\n', meta.parameter);
    end
    if isfield(meta, 'weight')
        fprintf(' solved middle weight w             = %.8f\n', meta.weight);
    end
    if isfield(meta, 'pointError')
        fprintf(' metadata shoulder fit error        = %.3e\n', meta.pointError);
    end
end

if isstruct(meta) && isfield(meta, 'parameter')
    Psh = C.evaluate(meta.parameter);
    err = norm(Psh - S);
    fprintf(' verify ||C(u)-S||                  = %.3e\n', err);
end

% -------------------------------------------------------------------------
% Plot 3D construction
% -------------------------------------------------------------------------
figure('Name','Single Liming conic');
hold on; grid on; axis equal; view(3);

% Curve
C.plot([], 'Color', [0.85 0.15 0.15], 'LineWidth', 2.0, 'ShowCP', true);

% Construction polygon / points
plot3([P0(1) T(1)], [P0(2) T(2)], [P0(3) T(3)], 'k--', 'LineWidth', 1.0);
plot3([P1(1) T(1)], [P1(2) T(2)], [P1(3) T(3)], 'k--', 'LineWidth', 1.0);

plot3(P0(1), P0(2), P0(3), 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 7);
plot3(P1(1), P1(2), P1(3), 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 7);
plot3(T(1),  T(2),  T(3),  'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 7);
plot3(S(1),  S(2),  S(3),  'ko', 'MarkerFaceColor', 'm', 'MarkerSize', 7);

text(P0(1), P0(2), P0(3), '  P0');
text(P1(1), P1(2), P1(3), '  P1');
text(T(1),  T(2),  T(3),  '  T');
text(S(1),  S(2),  S(3),  '  S');

xlabel('X');
ylabel('Y');
zlabel('Z');
title('Single Liming conic from P0, P1, T, S');

fprintf('\nDone.\n');