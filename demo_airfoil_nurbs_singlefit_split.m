%% DEMO_AIRFOIL_NURBS_SINGLEFIT_SPLIT.m
% Read an airfoil file in Selig-style ordering, fit one continuous NURBS
% curve, find the geometric leading edge from the fitted curve, split into
% upper/lower fitted curves, and plot the result.
%
% Expected airfoil ordering:
%   TE -> upper -> LE -> lower -> TE
%
% Main idea:
%   1) fit one continuous NURBS curve through the full airfoil point set
%   2) locate the true leading edge on the fitted curve by solving dx/du = 0
%   3) split the fitted curve at that parameter
%   4) reorient branches so both upper and lower run LE -> TE
%
% Notes:
%   - This does not force exact trailing-edge closure constraints.
%   - It is intended as a clean geometry demo for later lofting work.
%
% Usage:
%   edit airfoilFile below, then run:
%       demo_airfoil_nurbs_singlefit_split

clear; close all; clc;

%% ----------------------------------------------------------------------
% User settings
%% ----------------------------------------------------------------------
airfoilFile = 'airfoil.dat';   % change to your file name
pFit        = 3;                 % degree of continuous fitted airfoil curve
fitMethod   = 'centripetal';     % uniform | chord | centripetal | arc_length
nPlot       = 600;               % dense plot sampling
nSearch     = 2000;              % dense sampling for LE seed search
showArrows  = true;

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
% Raw leading-edge estimate from data
%% ----------------------------------------------------------------------
[~, leRawIdx] = min(xy(:,1));
leRawPt = xy(leRawIdx,:);

fprintf('\n--- 2. Raw leading-edge estimate ---\n');
fprintf('  raw LE index = %d\n', leRawIdx);
fprintf('  raw LE point = [%.8f %.8f]\n', leRawPt(1), leRawPt(2));

%% ----------------------------------------------------------------------
% Fit one continuous NURBS curve through the full airfoil
%% ----------------------------------------------------------------------
fprintf('\n--- 3. Continuous NURBS fit ---\n');

xy3 = [xy, zeros(size(xy,1),1)];
Cfull = geom.NURBSCurve.globalInterp(xy3, pFit, fitMethod);

fprintf('  full airfoil degree = %d\n', Cfull.p);
fprintf('  full airfoil #CP    = %d\n', size(Cfull.P,1));

uDense = linspace(Cfull.domain(1), Cfull.domain(2), nPlot);
xyFitDense3 = Cfull.evaluate(uDense);
xyFitDense = xyFitDense3(:,1:2);

uData = geom.NURBSCurve.parameterizeData(xy3, fitMethod, pFit);
xyRecon3 = Cfull.evaluate(uData);
rmsFull = sqrt(mean(sum((xyRecon3(:,1:2) - xy).^2, 2)));

fprintf('  full-curve RMS interpolation check = %.8e\n', rmsFull);

%% ----------------------------------------------------------------------
% Find geometric LE on fitted curve by slope condition dx/du = 0
%% ----------------------------------------------------------------------
fprintf('\n--- 4. Geometric leading edge from fitted curve ---\n');

[uLE, pLE3, dxduLE] = findLeadingEdgeBySlopeLocal(Cfull, nSearch);
pLE = pLE3(1:2);

fprintf('  fitted LE parameter = %.12f\n', uLE);
fprintf('  fitted LE point     = [%.8f %.8f]\n', pLE(1), pLE(2));
fprintf('  dx/du at LE         = %.8e\n', dxduLE);

%% ----------------------------------------------------------------------
% Split full fitted curve at the geometric LE
%% ----------------------------------------------------------------------
fprintf('\n--- 5. Split fitted curve at LE ---\n');

[Ca, Cb] = Cfull.split(uLE);

% Determine which branch is upper/lower using midpoint y-value
uaMid = 0.5 * (Ca.domain(1) + Ca.domain(2));
ubMid = 0.5 * (Cb.domain(1) + Cb.domain(2));

paMid = Ca.evaluate(uaMid);
pbMid = Cb.evaluate(ubMid);

if paMid(2) >= pbMid(2)
    Cup_raw = Ca;   % likely TE -> upper -> LE
    Clo_raw = Cb;   % likely LE -> lower -> TE
else
    Cup_raw = Cb;
    Clo_raw = Ca;
end

% Reorient so both go LE -> TE
% Upper branch from split will usually be TE -> LE, so reverse it.
Cup = Cup_raw.reverse();
Clo = Clo_raw;

fprintf('  upper branch degree = %d\n', Cup.p);
fprintf('  lower branch degree = %d\n', Clo.p);
fprintf('  upper #CP           = %d\n', size(Cup.P,1));
fprintf('  lower #CP           = %d\n', size(Clo.P,1));

%% ----------------------------------------------------------------------
% Sample split branches
%% ----------------------------------------------------------------------
uUp = linspace(Cup.domain(1), Cup.domain(2), nPlot);
uLo = linspace(Clo.domain(1), Clo.domain(2), nPlot);

xyUp3 = Cup.evaluate(uUp);
xyLo3 = Clo.evaluate(uLo);

xyUp = xyUp3(:,1:2);
xyLo = xyLo3(:,1:2);

% Endpoints for reporting
upLE = xyUp(1,:);
upTE = xyUp(end,:);
loLE = xyLo(1,:);
loTE = xyLo(end,:);

fprintf('  upper LE endpoint   = [%.8f %.8f]\n', upLE(1), upLE(2));
fprintf('  upper TE endpoint   = [%.8f %.8f]\n', upTE(1), upTE(2));
fprintf('  lower LE endpoint   = [%.8f %.8f]\n', loLE(1), loLE(2));
fprintf('  lower TE endpoint   = [%.8f %.8f]\n', loTE(1), loTE(2));

%% ----------------------------------------------------------------------
% Plot 1: raw points
%% ----------------------------------------------------------------------
fig1 = figure('Name','1 - Raw airfoil points');
hold on; grid on; axis equal;
plot(xy(:,1), xy(:,2), 'k.-', 'LineWidth', 1.0, 'MarkerSize', 10);
plot(leRawPt(1), leRawPt(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
xlabel('x'); ylabel('y');
title(sprintf('Raw airfoil points: %s', strtrim(nameLine)));
legend({'Raw points','Raw LE point'}, 'Location','best');

%% ----------------------------------------------------------------------
% Plot 2: full fitted curve + fitted LE
%% ----------------------------------------------------------------------
fig2 = figure('Name','2 - Full fitted airfoil');
hold on; grid on; axis equal;
plot(xy(:,1), xy(:,2), 'k.', 'MarkerSize', 8);
plot(xyFitDense(:,1), xyFitDense(:,2), '-', 'Color', [0.1 0.5 0.9], 'LineWidth', 2.0);
plot(Cfull.P(:,1), Cfull.P(:,2), '--', 'Color', [0.6 0.8 1.0], 'LineWidth', 1.0);
plot(Cfull.P(:,1), Cfull.P(:,2), 's', 'Color', [0.1 0.5 0.9], ...
    'MarkerFaceColor', 'w', 'MarkerSize', 5);
plot(pLE(1), pLE(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
xlabel('x'); ylabel('y');
title(sprintf('Continuous NURBS fit (%s parameterization)', fitMethod));
legend({'Raw data','Full fitted curve','Control polygon','Control points','Fitted LE'}, ...
       'Location','best');

%% ----------------------------------------------------------------------
% Plot 3: split upper/lower fitted curves
%% ----------------------------------------------------------------------
fig3 = figure('Name','3 - Split upper/lower fitted curves');
hold on; grid on; axis equal;

plot(xy(:,1), xy(:,2), 'k.', 'MarkerSize', 8);
plot(xyUp(:,1), xyUp(:,2), '-', 'Color', [0.9 0.15 0.15], 'LineWidth', 2.0);
plot(xyLo(:,1), xyLo(:,2), '-', 'Color', [0.15 0.15 0.9], 'LineWidth', 2.0);

plot(Cup.P(:,1), Cup.P(:,2), '--', 'Color', [1.0 0.6 0.6], 'LineWidth', 1.0);
plot(Clo.P(:,1), Clo.P(:,2), '--', 'Color', [0.6 0.6 1.0], 'LineWidth', 1.0);

plot(Cup.P(:,1), Cup.P(:,2), 's', 'Color', [0.9 0.15 0.15], ...
    'MarkerFaceColor', 'w', 'MarkerSize', 5);
plot(Clo.P(:,1), Clo.P(:,2), 's', 'Color', [0.15 0.15 0.9], ...
    'MarkerFaceColor', 'w', 'MarkerSize', 5);

plot(pLE(1), pLE(2), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8);

if showArrows
    addDirectionArrowLocal(xyUp, 0.20, [0.9 0.15 0.15]);
    addDirectionArrowLocal(xyLo, 0.20, [0.15 0.15 0.9]);
end

xlabel('x'); ylabel('y');
title('Split fitted airfoil branches (both LE -> TE)');
legend({'Raw points','Upper fit','Lower fit', ...
        'Upper CP poly','Lower CP poly','Upper CP','Lower CP','Geometric LE'}, ...
       'Location','best');

%% ----------------------------------------------------------------------
% Plot 4: branch endpoints and directions
%% ----------------------------------------------------------------------
fig4 = figure('Name','4 - Branch endpoints and directions');
hold on; grid on; axis equal;

plot(xyUp(:,1), xyUp(:,2), '-', 'Color', [0.9 0.15 0.15], 'LineWidth', 2.0);
plot(xyLo(:,1), xyLo(:,2), '-', 'Color', [0.15 0.15 0.9], 'LineWidth', 2.0);

plot(upLE(1), upLE(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
plot(loLE(1), loLE(2), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 7);
plot(upTE(1), upTE(2), 'rs', 'MarkerFaceColor', 'w', 'MarkerSize', 7);
plot(loTE(1), loTE(2), 'bs', 'MarkerFaceColor', 'w', 'MarkerSize', 7);

if showArrows
    addDirectionArrowLocal(xyUp, 0.30, [0.9 0.15 0.15]);
    addDirectionArrowLocal(xyLo, 0.30, [0.15 0.15 0.9]);
end

xlabel('x'); ylabel('y');
title('Upper and lower fitted curves, both oriented LE -> TE');
legend({'Upper fit','Lower fit','Upper LE','Lower LE','Upper TE','Lower TE'}, ...
       'Location','best');

%% ----------------------------------------------------------------------
% Summary
%% ----------------------------------------------------------------------
fprintf('\n=== Single-Fit Airfoil NURBS Demo Complete ===\n');
fprintf('Figures:\n');
fprintf('  (1) Raw airfoil points\n');
fprintf('  (2) Full fitted airfoil with geometric LE\n');
fprintf('  (3) Split upper/lower fitted curves\n');
fprintf('  (4) Branch endpoints and directions\n');

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

function [uLE, pLE, dxduLE] = findLeadingEdgeBySlopeLocal(C, nSearch)
    % Find the geometric leading edge by solving dx/du = 0 and selecting
    % the minimum-x solution.

    us = linspace(C.domain(1), C.domain(2), nSearch);
    pts = C.evaluate(us);
    d1  = C.derivative(us, 1);

    xvals = pts(:,1);
    dxdu  = d1(:,1);

    % Candidate seeds:
    %   - sign changes in dx/du
    %   - smallest sampled x
    cand = [];

    for i = 1:numel(us)-1
        if dxdu(i) == 0
            cand(end+1) = us(i); %#ok<AGROW>
        elseif dxdu(i) * dxdu(i+1) < 0
            cand(end+1) = 0.5 * (us(i) + us(i+1)); %#ok<AGROW>
        end
    end

    [~, idxMin] = min(xvals);
    cand(end+1) = us(idxMin); %#ok<AGROW>
    cand = unique(cand);

    if isempty(cand)
        uLE = us(idxMin);
        pLE = C.evaluate(uLE);
        dxduLE = C.derivative(uLE,1);
        dxduLE = dxduLE(1);
        return;
    end

    bestU = cand(1);
    bestX = inf;

    for k = 1:numel(cand)
        uk = refineDxDuRootLocal(C, cand(k));

        pk = C.evaluate(uk);
        xk = pk(1);

        if xk < bestX
            bestX = xk;
            bestU = uk;
        end
    end

    uLE = bestU;
    pLE = C.evaluate(uLE);
    dLE = C.derivative(uLE,1);
    dxduLE = dLE(1);
end

function u = refineDxDuRootLocal(C, u0)
    % Newton refine root of dx/du = 0 using x''(u).

    u = u0;
    for iter = 1:40
        d1 = C.derivative(u,1);
        d2 = C.derivative(u,2);

        f  = d1(1);
        df = d2(1);

        if abs(df) < 1e-14
            break;
        end

        un = u - f / df;
        un = max(C.domain(1), min(C.domain(2), un));

        if abs(un - u) < 1e-13
            u = un;
            break;
        end

        u = un;
    end
end

function addDirectionArrowLocal(xy, frac, colorIn)
    if nargin < 2 || isempty(frac), frac = 0.25; end
    if nargin < 3 || isempty(colorIn), colorIn = [0 0 0]; end

    n = size(xy,1);
    i0 = max(1, min(n-1, round(frac*(n-1))));
    i1 = min(n, i0 + max(3, round(0.03*n)));

    p0 = xy(i0,:);
    p1 = xy(i1,:);
    dp = p1 - p0;

    quiver(p0(1), p0(2), dp(1), dp(2), 0, ...
        'Color', colorIn, 'LineWidth', 1.5, 'MaxHeadSize', 1.5);
end