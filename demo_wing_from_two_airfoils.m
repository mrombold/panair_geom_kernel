%% DEMO_WING_FROM_TWO_AIRFOILS.m
% Build a simple wing from root and tip airfoils.
%
% Inputs:
%   - root airfoil file
%   - tip airfoil file
%   - root chord
%   - tip chord
%   - leading-edge sweep
%   - whole-wing incidence
%   - tip twist relative to root
%   - dihedral
%
% Airfoils are assumed to be defined in local 2D with x in [0,1].
%
% The script:
%   1) reads root and tip airfoil files
%   2) fits one continuous airfoil NURBS to each
%   3) finds geometric LE and splits into upper/lower fitted curves
%   4) orients upper/lower as LE -> TE
%   5) scales by chord
%   6) places root and tip in 3D using sweep / incidence / twist / dihedral
%   7) lofts upper and lower wing surfaces separately

clear; close all; clc;

%% ----------------------------------------------------------------------
% User settings
%% ----------------------------------------------------------------------
rootAirfoilFile = 'airfoil_sharpTE.dat';
tipAirfoilFile  = 'airfoil_sharpTE.dat';

pFit       = 3;              % airfoil fit degree
fitMethod  = 'centripetal';  % uniform | chord | centripetal | arc_length

span       = 10.0;           % full span if symmetric? here we build one semispan side length directly
semispan   = 5.0;            % tip station in +Y

rootChord  = 2.0;
tipChord   = 1.0;

leSweepDeg = 20.0;           % leading-edge sweep angle
incDeg     = 2.0;            % whole-wing incidence
tipTwistDeg = -3.0;          % tip relative to root
dihedralDeg = 5.0;           % wing dihedral

loftDegree = 1;              % only two sections, so degree must be 1
nPlotCurve = 500;

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
% Read and fit root / tip airfoils
%% ----------------------------------------------------------------------
fprintf('\n--- 1. Read and fit root / tip airfoils ---\n');

rootSection = buildAirfoilSectionLocal(rootAirfoilFile, pFit, fitMethod);
tipSection  = buildAirfoilSectionLocal(tipAirfoilFile,  pFit, fitMethod);

fprintf('  root file: %s\n', rootAirfoilFile);
fprintf('  tip  file: %s\n', tipAirfoilFile);
fprintf('  root fitted LE = [%.8f %.8f]\n', rootSection.LE(1), rootSection.LE(2));
fprintf('  tip  fitted LE = [%.8f %.8f]\n', tipSection.LE(1),  tipSection.LE(2));

%% ----------------------------------------------------------------------
% Build root and tip transforms
%% ----------------------------------------------------------------------
fprintf('\n--- 2. Build section placement transforms ---\n');

sweepRad    = deg2rad(leSweepDeg);
incRad      = deg2rad(incDeg);
tipTwistRad = deg2rad(tipTwistDeg);
dihedRad    = deg2rad(dihedralDeg);

% Root leading edge at origin
rootLE = [0, 0, 0];

% Tip leading edge from sweep + dihedral
tipLE = [semispan * tan(sweepRad), ...
         semispan, ...
         semispan * tan(dihedRad)];

fprintf('  root LE = [%.6f %.6f %.6f]\n', rootLE);
fprintf('  tip  LE = [%.6f %.6f %.6f]\n', tipLE);

Troot = makeWingSectionTransformLocal(rootChord, incRad, rootLE);
Ttip  = makeWingSectionTransformLocal(tipChord, incRad + tipTwistRad, tipLE);

%% ----------------------------------------------------------------------
% Transform section curves into 3D
%% ----------------------------------------------------------------------
fprintf('\n--- 3. Transform section curves into 3D ---\n');

Cup_root = rootSection.upper.transform(Troot);
Clo_root = rootSection.lower.transform(Troot);

Cup_tip  = tipSection.upper.transform(Ttip);
Clo_tip  = tipSection.lower.transform(Ttip);

fprintf('  transformed root upper #CP = %d\n', size(Cup_root.P,1));
fprintf('  transformed tip  upper #CP = %d\n', size(Cup_tip.P,1));

%% ----------------------------------------------------------------------
% Loft upper and lower wing surfaces
%% ----------------------------------------------------------------------
fprintf('\n--- 4. Loft upper and lower wing surfaces ---\n');

Supper = geom.NURBSSurface.loft({Cup_root, Cup_tip}, loftDegree, 'chord');
Slower = geom.NURBSSurface.loft({Clo_root, Clo_tip}, loftDegree, 'chord');

fprintf('  upper surface net size = %d x %d\n', size(Supper.P,1), size(Supper.P,2));
fprintf('  lower surface net size = %d x %d\n', size(Slower.P,1), size(Slower.P,2));

%% ----------------------------------------------------------------------
% Sample curves for plotting
%% ----------------------------------------------------------------------
uRup = linspace(Cup_root.domain(1), Cup_root.domain(2), nPlotCurve);
uRlo = linspace(Clo_root.domain(1), Clo_root.domain(2), nPlotCurve);
uTup = linspace(Cup_tip.domain(1),  Cup_tip.domain(2),  nPlotCurve);
uTlo = linspace(Clo_tip.domain(1),  Clo_tip.domain(2),  nPlotCurve);

Rup = Cup_root.evaluate(uRup);
Rlo = Clo_root.evaluate(uRlo);
Tup = Cup_tip.evaluate(uTup);
Tlo = Clo_tip.evaluate(uTlo);

%% ----------------------------------------------------------------------
% Plot 1: 2D root and tip fitted sections
%% ----------------------------------------------------------------------
fig1 = figure('Name','1 - Root and tip fitted airfoils');
hold on; grid on; axis equal;

plot(rootSection.fullFit(:,1), rootSection.fullFit(:,2), '-', ...
    'Color', [0.85 0.15 0.15], 'LineWidth', 1.8);
plot(tipSection.fullFit(:,1), tipSection.fullFit(:,2), '-', ...
    'Color', [0.15 0.15 0.85], 'LineWidth', 1.8);

plot(rootSection.LE(1), rootSection.LE(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
plot(tipSection.LE(1), tipSection.LE(2), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 7);

xlabel('x'); ylabel('y');
title('Root and tip fitted airfoils in local 2D');
legend({'Root full fit','Tip full fit','Root LE','Tip LE'}, 'Location','best');

%% ----------------------------------------------------------------------
% Plot 2: transformed root / tip section curves in 3D
%% ----------------------------------------------------------------------
fig2 = figure('Name','2 - Root and tip sections in 3D');
hold on; grid on; axis equal;

plot3(Rup(:,1), Rup(:,2), Rup(:,3), '-', 'Color', [0.9 0.15 0.15], 'LineWidth', 2.0);
plot3(Rlo(:,1), Rlo(:,2), Rlo(:,3), '-', 'Color', [0.15 0.15 0.9], 'LineWidth', 2.0);
plot3(Tup(:,1), Tup(:,2), Tup(:,3), '-', 'Color', [0.9 0.5 0.15], 'LineWidth', 2.0);
plot3(Tlo(:,1), Tlo(:,2), Tlo(:,3), '-', 'Color', [0.15 0.7 0.7], 'LineWidth', 2.0);

plot3(rootLE(1), rootLE(2), rootLE(3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
plot3(tipLE(1), tipLE(2), tipLE(3), 'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 7);

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Root and tip airfoil sections positioned in 3D');
legend({'Root upper','Root lower','Tip upper','Tip lower','Root LE','Tip LE'}, ...
       'Location','best');
view(3);

%% ----------------------------------------------------------------------
% Plot 3: upper and lower lofted wing surfaces
%% ----------------------------------------------------------------------
fig3 = figure('Name','3 - Lofted wing upper/lower surfaces');
hold on; grid on; axis equal;

Supper.plot(32, 18, 'ShowCP', true, 'ShowIso', true, ...
    'FaceColor', [0.85 0.25 0.25], 'Alpha', 0.78);
Slower.plot(32, 18, 'ShowCP', true, 'ShowIso', true, ...
    'FaceColor', [0.25 0.35 0.90], 'Alpha', 0.78);

plot3(Rup(:,1), Rup(:,2), Rup(:,3), 'k-', 'LineWidth', 1.2);
plot3(Rlo(:,1), Rlo(:,2), Rlo(:,3), 'k-', 'LineWidth', 1.2);
plot3(Tup(:,1), Tup(:,2), Tup(:,3), 'k-', 'LineWidth', 1.2);
plot3(Tlo(:,1), Tlo(:,2), Tup(:,3)*0 + Tlo(:,3), 'k-', 'LineWidth', 1.2);

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Lofted wing surfaces (upper and lower)');
view(3);

%% ----------------------------------------------------------------------
% Plot 4: wing with leading-edge line and quarter-chord reference
%% ----------------------------------------------------------------------
fig4 = figure('Name','4 - Wing reference lines');
hold on; grid on; axis equal;

Supper.plot(28, 16, 'ShowCP', false, 'ShowIso', false, ...
    'FaceColor', [0.90 0.30 0.30], 'Alpha', 0.72);
Slower.plot(28, 16, 'ShowCP', false, 'ShowIso', false, ...
    'FaceColor', [0.30 0.40 0.95], 'Alpha', 0.72);

% LE line
plot3([rootLE(1), tipLE(1)], [rootLE(2), tipLE(2)], [rootLE(3), tipLE(3)], ...
    'k--', 'LineWidth', 1.5);

% Quarter-chord points
rootQC = rootLE + [0.25*rootChord*cos(incRad), 0, 0.25*rootChord*sin(incRad)];
tipQC  = tipLE  + [0.25*tipChord*cos(incRad + tipTwistRad), 0, 0.25*tipChord*sin(incRad + tipTwistRad)];

plot3([rootQC(1), tipQC(1)], [rootQC(2), tipQC(2)], [rootQC(3), tipQC(3)], ...
    'm-.', 'LineWidth', 1.5);

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Wing with LE and quarter-chord reference lines');
legend({'Upper surface','Lower surface','LE line','Quarter-chord line'}, ...
       'Location','best');
view(3);

%% ----------------------------------------------------------------------
% Summary
%% ----------------------------------------------------------------------
fprintf('\n=== Wing Loft Demo Complete ===\n');
fprintf('Inputs:\n');
fprintf('  semispan        = %.6f\n', semispan);
fprintf('  root chord      = %.6f\n', rootChord);
fprintf('  tip chord       = %.6f\n', tipChord);
fprintf('  LE sweep (deg)  = %.6f\n', leSweepDeg);
fprintf('  incidence (deg) = %.6f\n', incDeg);
fprintf('  tip twist (deg) = %.6f\n', tipTwistDeg);
fprintf('  dihedral (deg)  = %.6f\n', dihedralDeg);

fprintf('\nFigures:\n');
fprintf('  (1) Root and tip fitted airfoils in local 2D\n');
fprintf('  (2) Root and tip sections in 3D\n');
fprintf('  (3) Lofted upper/lower wing surfaces\n');
fprintf('  (4) Wing reference lines\n');

%% ======================================================================
% Local functions
%% ======================================================================

function section = buildAirfoilSectionLocal(filename, pFit, fitMethod)
    [titleLine, xy] = readAirfoilFileLocal(filename); %#ok<NASGU>

    xy3 = [xy, zeros(size(xy,1),1)];
    Cfull = geom.NURBSCurve.globalInterp(xy3, pFit, fitMethod);

    [uLE, pLE3] = findLeadingEdgeBySlopeLocal(Cfull, 2000); %#ok<ASGLU>
    LE = pLE3(1:2);

    [Ca, Cb] = Cfull.split(uLE);

    uaMid = 0.5 * (Ca.domain(1) + Ca.domain(2));
    ubMid = 0.5 * (Cb.domain(1) + Cb.domain(2));

    paMid = Ca.evaluate(uaMid);
    pbMid = Cb.evaluate(ubMid);

    if paMid(2) >= pbMid(2)
        Cup_raw = Ca;
        Clo_raw = Cb;
    else
        Cup_raw = Cb;
        Clo_raw = Ca;
    end

    % Reorient both LE -> TE
    Cup = Cup_raw.reverse();
    Clo = Clo_raw;

    uFull = linspace(Cfull.domain(1), Cfull.domain(2), 600);
    fullFit3 = Cfull.evaluate(uFull);

    section.full  = Cfull;
    section.upper = Cup;
    section.lower = Clo;
    section.LE    = LE;
    section.fullFit = fullFit3(:,1:2);
end

function T = makeWingSectionTransformLocal(chord, alpha, LE_xyz)
    % Local airfoil coordinates:
    %   x = chordwise
    %   y = thickness direction in 2D file
    %
    % Map local 2D airfoil into global 3D as:
    %   local x -> global X/Z plane via incidence rotation about Y
    %   local y -> global Z before incidence
    %
    % So the section stays in a plane of constant global Y.

    c = cos(alpha);
    s = sin(alpha);

    % Start from scaling
    S = eye(4);
    S(1,1) = chord;
    S(2,2) = chord;
    S(3,3) = chord;

    % Embed 2D airfoil into X-Z by mapping local y -> global z.
    % Since NURBSCurve.transform applies 4x4 to [x y z 1],
    % and airfoil has z=0, we use a matrix that moves local y into Z.
    E = eye(4);
    E(2,2) = 0;
    E(3,2) = 1;
    E(3,3) = 0;

    % Incidence rotation about global Y
    Ry = eye(4);
    Ry(1,1) =  c;
    Ry(1,3) =  s;
    Ry(3,1) = -s;
    Ry(3,3) =  c;

    % Translation to section LE location
    Tr = eye(4);
    Tr(1:3,4) = LE_xyz(:);

    T = Tr * Ry * E * S;
end

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
    us = linspace(C.domain(1), C.domain(2), nSearch);
    pts = C.evaluate(us);
    d1  = C.derivative(us, 1);

    xvals = pts(:,1);
    dxdu  = d1(:,1);

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