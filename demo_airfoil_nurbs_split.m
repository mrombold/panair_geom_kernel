%% DEMO_AIRFOIL_NURBS_SPLIT.m
% Read an airfoil coordinate file, split upper/lower surfaces,
% fit separate NURBS curves, and plot the result.
%
% Expected airfoil file format:
%   first line: optional title
%   remaining lines: x y
%
% Typical ordering:
%   TE -> upper -> LE -> lower -> TE
%
% This script:
%   1. reads the file
%   2. finds the leading edge from minimum x
%   3. splits upper and lower point sets
%   4. fits one NURBS curve to each side
%   5. plots raw points and fitted curves
%
% Usage:
%   edit airfoilFile below, then run:
%       demo_airfoil_nurbs_split

clear; close all; clc;

%% ----------------------------------------------------------------------
% User settings
%% ----------------------------------------------------------------------
airfoilFile = ['airfoil_sharpTE.dat'];   % change to your file name
pFit        = 3;                 % NURBS degree for fitted curves
fitMethod   = 'centripetal';     % uniform | chord | centripetal | arc_length

%% ----------------------------------------------------------------------
% Path / package check
%% ----------------------------------------------------------------------
script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir)
    addpath(script_dir);
end

if ~exist('+geom', 'dir') && isempty(which('geom.NURBSCurve'))
    error(['Cannot find +geom package.\n' ...
           'Expected it near: %s\n' ...
           'Add the geometry kernel folder to your MATLAB path.'], ...
           fullfile(script_dir, '+geom'));
end

fprintf('geom package found at: %s\n', script_dir);

%% ----------------------------------------------------------------------
% Read airfoil file
%% ----------------------------------------------------------------------
[nameLine, xy] = readAirfoilFileLocal(airfoilFile);

fprintf('\n--- 1. File read ---\n');
fprintf('  file: %s\n', airfoilFile);
fprintf('  title: %s\n', strtrim(nameLine));
fprintf('  # coordinate rows: %d\n', size(xy,1));

if size(xy,1) < 10
    error('Too few airfoil coordinates were read.');
end

%% ----------------------------------------------------------------------
% Split upper / lower surfaces
%% ----------------------------------------------------------------------
[upperPts, lowerPts, leIdx] = splitAirfoilUpperLowerLocal(xy);

fprintf('\n--- 2. Surface split ---\n');
fprintf('  leading-edge index = %d\n', leIdx);
fprintf('  # upper points = %d\n', size(upperPts,1));
fprintf('  # lower points = %d\n', size(lowerPts,1));

% Reorder each side LE -> TE for fitting convenience
upperFitPts = flipud(upperPts);   % upper was TE->LE, so flip to LE->TE
lowerFitPts = lowerPts;           % lower already LE->TE after split

% Remove duplicate points if present
upperFitPts = uniqueTolRowsLocal(upperFitPts, 1e-12);
lowerFitPts = uniqueTolRowsLocal(lowerFitPts, 1e-12);

fprintf('  # upper fit points (deduped) = %d\n', size(upperFitPts,1));
fprintf('  # lower fit points (deduped) = %d\n', size(lowerFitPts,1));

%% ----------------------------------------------------------------------
% Build NURBS curves
%% ----------------------------------------------------------------------
fprintf('\n--- 3. NURBS fitting ---\n');

Cup = geom.NURBSCurve.globalInterp(upperFitPts, pFit, fitMethod);
Clo = geom.NURBSCurve.globalInterp(lowerFitPts, pFit, fitMethod);

fprintf('  upper curve degree = %d\n', Cup.p);
fprintf('  lower curve degree = %d\n', Clo.p);
fprintf('  upper #CP = %d\n', size(Cup.P,1));
fprintf('  lower #CP = %d\n', size(Clo.P,1));

%% ----------------------------------------------------------------------
% Evaluate fitted curves
%% ----------------------------------------------------------------------
nu = 400;
uUp = linspace(Cup.domain(1), Cup.domain(2), nu);
uLo = linspace(Clo.domain(1), Clo.domain(2), nu);

xyUpFit = Cup.evaluate(uUp);
xyLoFit = Clo.evaluate(uLo);

%% ----------------------------------------------------------------------
% Optional fit error at original sample points
%% ----------------------------------------------------------------------
uDataUp = geom.NURBSCurve.parameterizeData(upperFitPts, fitMethod, pFit);
uDataLo = geom.NURBSCurve.parameterizeData(lowerFitPts, fitMethod, pFit);

upRecon = Cup.evaluate(uDataUp);
loRecon = Clo.evaluate(uDataLo);

upErr = sqrt(mean(sum((upRecon(:,1:2) - upperFitPts(:,1:2)).^2, 2)));
loErr = sqrt(mean(sum((loRecon(:,1:2) - lowerFitPts(:,1:2)).^2, 2)));

fprintf('  upper RMS interpolation check = %.8e\n', upErr);
fprintf('  lower RMS interpolation check = %.8e\n', loErr);

%% ----------------------------------------------------------------------
% Plots
%% ----------------------------------------------------------------------
fprintf('\n--- 4. Plotting ---\n');

fig1 = figure('Name','1 - Raw airfoil points');
hold on; grid on; axis equal;
plot(xy(:,1), xy(:,2), 'k.-', 'LineWidth', 1.0, 'MarkerSize', 10);
plot(xy(leIdx,1), xy(leIdx,2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
xlabel('x'); ylabel('y');
title(sprintf('Raw airfoil points: %s', strtrim(nameLine)));
legend({'Raw points','Detected LE'}, 'Location','best');

fig2 = figure('Name','2 - Upper / lower split');
hold on; grid on; axis equal;
plot(upperPts(:,1), upperPts(:,2), 'r.-', 'LineWidth', 1.2, 'MarkerSize', 10);
plot(lowerPts(:,1), lowerPts(:,2), 'b.-', 'LineWidth', 1.2, 'MarkerSize', 10);
xlabel('x'); ylabel('y');
title('Upper / lower split');
legend({'Upper raw segment','Lower raw segment'}, 'Location','best');

fig3 = figure('Name','3 - NURBS fit');
hold on; grid on; axis equal;

plot(upperFitPts(:,1), upperFitPts(:,2), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
plot(lowerFitPts(:,1), lowerFitPts(:,2), 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b');

plot(xyUpFit(:,1), xyUpFit(:,2), '-', 'Color', [0.85 0.15 0.15], 'LineWidth', 2.0);
plot(xyLoFit(:,1), xyLoFit(:,2), '-', 'Color', [0.15 0.15 0.85], 'LineWidth', 2.0);

plot(Cup.P(:,1), Cup.P(:,2), '--', 'Color', [1.0 0.5 0.5], 'LineWidth', 1.0);
plot(Clo.P(:,1), Clo.P(:,2), '--', 'Color', [0.5 0.5 1.0], 'LineWidth', 1.0);

plot(Cup.P(:,1), Cup.P(:,2), 's', 'Color', [0.85 0.15 0.15], ...
    'MarkerFaceColor', 'w', 'MarkerSize', 5);
plot(Clo.P(:,1), Clo.P(:,2), 's', 'Color', [0.15 0.15 0.85], ...
    'MarkerFaceColor', 'w', 'MarkerSize', 5);

xlabel('x'); ylabel('y');
title(sprintf('Upper / lower NURBS curves (%s parameterization)', fitMethod));
legend({'Upper data','Lower data','Upper NURBS','Lower NURBS', ...
        'Upper CP poly','Lower CP poly','Upper CP','Lower CP'}, ...
       'Location','best');

fig4 = figure('Name','4 - Overlay with original ordering');
hold on; grid on; axis equal;
plot(xy(:,1), xy(:,2), 'k.', 'MarkerSize', 8);
plot(xyUpFit(:,1), xyUpFit(:,2), '-', 'Color', [0.85 0.15 0.15], 'LineWidth', 2.0);
plot(xyLoFit(:,1), xyLoFit(:,2), '-', 'Color', [0.15 0.15 0.85], 'LineWidth', 2.0);
xlabel('x'); ylabel('y');
title('Original points with separate upper/lower NURBS fits');
legend({'Original file points','Upper fit','Lower fit'}, 'Location','best');

%% ----------------------------------------------------------------------
% Summary
%% ----------------------------------------------------------------------
fprintf('\n=== Airfoil NURBS Demo Complete ===\n');
fprintf('Figures:\n');
fprintf('  (1) Raw airfoil points\n');
fprintf('  (2) Upper / lower split\n');
fprintf('  (3) Separate NURBS fits\n');
fprintf('  (4) Overlay with original points\n');

%% ======================================================================
% Local functions
%% ======================================================================

function [titleLine, xy] = readAirfoilFileLocal(filename)
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end

    c = onCleanup(@() fclose(fid));

    lines = {};
    while true
        t = fgetl(fid);
        if ~ischar(t)
            break;
        end
        lines{end+1,1} = t; %#ok<AGROW>
    end

    if isempty(lines)
        error('Airfoil file is empty.');
    end

    titleLine = lines{1};
    data = [];

    for i = 2:numel(lines)
        s = strtrim(lines{i});
        if isempty(s)
            continue;
        end

        vals = sscanf(s, '%f %f');
        if numel(vals) >= 2
            data(end+1,:) = vals(1:2).'; %#ok<AGROW>
        end
    end

    % If nothing parsed after line 1, try parsing all lines as numeric
    if isempty(data)
        for i = 1:numel(lines)
            s = strtrim(lines{i});
            if isempty(s)
                continue;
            end
            vals = sscanf(s, '%f %f');
            if numel(vals) >= 2
                data(end+1,:) = vals(1:2).'; %#ok<AGROW>
            end
        end
        if isempty(data)
            error('No numeric x-y coordinate pairs found in file.');
        end
    end

    xy = data;
end

function [upperPts, lowerPts, leIdx] = splitAirfoilUpperLowerLocal(xy)
    % Typical ordering:
    %   TE -> upper -> LE -> lower -> TE
    %
    % We detect LE as minimum x. If several points tie, take the first.
    [~, leIdx] = min(xy(:,1));

    if leIdx <= 1 || leIdx >= size(xy,1)
        error('Could not find a usable leading-edge split point.');
    end

    upperPts = xy(1:leIdx, :);      % TE -> LE
    lowerPts = xy(leIdx:end, :);    % LE -> TE
end

function A = uniqueTolRowsLocal(A, tol)
    if isempty(A)
        return;
    end

    keep = true(size(A,1),1);
    for i = 2:size(A,1)
        if norm(A(i,:) - A(i-1,:)) <= tol
            keep(i) = false;
        end
    end
    A = A(keep,:);
end