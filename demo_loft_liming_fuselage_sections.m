%% demo_loft_liming_fuselage_sections.m
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


% -------------------------------------------------------------------------
% Plot 3D construction
% -------------------------------------------------------------------------
figure('Name','Single Liming conic');
hold on; grid on; axis equal; view(3);

% Curve
C.plot([], 'Color', [0.85 0.15 0.15], 'LineWidth', 2.0, 'ShowCP', true);


xlabel('X');
ylabel('Y');
zlabel('Z');
title('Single Liming conic from P0, P1, T, S');

fprintf('\nDone.\n');