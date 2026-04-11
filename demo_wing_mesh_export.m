%% DEMO_WING_MESH_EXPORT.m
% Build a simple wing from root and tip airfoils, loft upper/lower surfaces,
% mesh them, and export structured meshes through geom.MeshWriter.

clear; close all; clc;

%% ----------------------------------------------------------------------
% User settings
%% ----------------------------------------------------------------------
rootAirfoilFile = 'airfoil_sharpTE.dat';
tipAirfoilFile  = 'airfoil_sharpTE.dat';

pFit       = 3;
fitMethod  = 'centripetal';

semispan    = 5.0;
rootChord   = 2.0;
tipChord    = 1.0;
leSweepDeg  = 20.0;
incDeg      = 2.0;
tipTwistDeg = -3.0;
dihedralDeg = 5.0;

loftDegree = 1;   % two sections only

% Structured mesh density
nu = 41;   % chordwise
nv = 21;   % spanwise

% Export switches
doCSV = true;
doVTK = true;
doSTL = true;
doWGS = true;

% Output names
csvUpper = 'wing_upper_mesh.csv';
csvLower = 'wing_lower_mesh.csv';

vtkUpper = 'wing_upper.vtk';
vtkLower = 'wing_lower.vtk';

stlUpper = 'wing_upper_mesh.stl';
stlLower = 'wing_lower_mesh.stl';

wgsFile  = 'wing_upper_lower.wgs';

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
% Section placement transforms
%% ----------------------------------------------------------------------
fprintf('\n--- 2. Build section placement transforms ---\n');

sweepRad    = deg2rad(leSweepDeg);
incRad      = deg2rad(incDeg);
tipTwistRad = deg2rad(tipTwistDeg);
dihedRad    = deg2rad(dihedralDeg);

rootLE = [0, 0, 0];
tipLE  = [semispan * tan(sweepRad), ...
          semispan, ...
          semispan * tan(dihedRad)];

fprintf('  root LE = [%.6f %.6f %.6f]\n', rootLE);
fprintf('  tip  LE = [%.6f %.6f %.6f]\n', tipLE);

Troot = makeWingSectionTransformLocal(rootChord, incRad, rootLE);
Ttip  = makeWingSectionTransformLocal(tipChord, incRad + tipTwistRad, tipLE);

%% ----------------------------------------------------------------------
% Transform upper / lower section curves into 3D
%% ----------------------------------------------------------------------
fprintf('\n--- 3. Transform section curves into 3D ---\n');

Cup_root = rootSection.upper.transform(Troot);
Clo_root = rootSection.lower.transform(Troot);

Cup_tip  = tipSection.upper.transform(Ttip);
Clo_tip  = tipSection.lower.transform(Ttip);

fprintf('  root upper #CP = %d\n', size(Cup_root.P,1));
fprintf('  tip  upper #CP = %d\n', size(Cup_tip.P,1));
fprintf('  root lower #CP = %d\n', size(Clo_root.P,1));
fprintf('  tip  lower #CP = %d\n', size(Clo_tip.P,1));

%% ----------------------------------------------------------------------
% Loft upper / lower wing surfaces
%% ----------------------------------------------------------------------
fprintf('\n--- 4. Loft upper and lower wing surfaces ---\n');

Supper = geom.NURBSSurface.loft({Cup_root, Cup_tip}, loftDegree, 'chord');
Slower = geom.NURBSSurface.loft({Clo_root, Clo_tip}, loftDegree, 'chord');

fprintf('  upper surface net size = %d x %d\n', size(Supper.P,1), size(Supper.P,2));
fprintf('  lower surface net size = %d x %d\n', size(Slower.P,1), size(Slower.P,2));

%% ----------------------------------------------------------------------
% Mesh upper / lower surfaces
%% ----------------------------------------------------------------------
fprintf('\n--- 5. Mesh upper and lower wing surfaces ---\n');

meshUpper = Supper.isoMesh(nu, nv);
meshLower = Slower.isoMesh(nu, nv);

fprintf('  upper mesh: nu = %d, nv = %d, quads = %d\n', ...
    meshUpper.nu, meshUpper.nv, size(meshUpper.connectivity,1));
fprintf('  lower mesh: nu = %d, nv = %d, quads = %d\n', ...
    meshLower.nu, meshLower.nv, size(meshLower.connectivity,1));

TEu = [meshUpper.X(end,:).', meshUpper.Y(end,:).', meshUpper.Z(end,:).'];
TEl = [meshLower.X(end,:).', meshLower.Y(end,:).', meshLower.Z(end,:).'];
dTE = vecnorm(TEu - TEl, 2, 2);

fprintf('  TE mesh mismatch: max = %.3e, RMS = %.3e\n', ...
    max(dTE), sqrt(mean(dTE.^2)));

%% ----------------------------------------------------------------------
% Plot 1: section placement
%% ----------------------------------------------------------------------
uPlot = linspace(Cup_root.domain(1), Cup_root.domain(2), 400);
Rup = Cup_root.evaluate(uPlot);
Rlo = Clo_root.evaluate(uPlot);
Tup = Cup_tip.evaluate(uPlot);
Tlo = Clo_tip.evaluate(uPlot);

fig1 = figure('Name','1 - Root and tip sections in 3D');
hold on; grid on; axis equal;
plot3(Rup(:,1), Rup(:,2), Rup(:,3), '-', 'Color', [0.90 0.15 0.15], 'LineWidth', 2.0);
plot3(Rlo(:,1), Rlo(:,2), Rlo(:,3), '-', 'Color', [0.15 0.15 0.90], 'LineWidth', 2.0);
plot3(Tup(:,1), Tup(:,2), Tup(:,3), '-', 'Color', [0.90 0.55 0.15], 'LineWidth', 2.0);
plot3(Tlo(:,1), Tlo(:,2), Tlo(:,3), '-', 'Color', [0.15 0.75 0.75], 'LineWidth', 2.0);
plot3(rootLE(1), rootLE(2), rootLE(3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
plot3(tipLE(1), tipLE(2), tipLE(3), 'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 7);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Root and tip sections in 3D');
legend({'Root upper','Root lower','Tip upper','Tip lower','Root LE','Tip LE'}, ...
       'Location','best');
view(3);

%% ----------------------------------------------------------------------
% Plot 2: lofted wing surfaces
%% ----------------------------------------------------------------------
fig2 = figure('Name','2 - Lofted wing surfaces');
hold on; grid on; axis equal;
Supper.plot(32, 18, 'ShowCP', true, 'ShowIso', true, ...
    'FaceColor', [0.85 0.25 0.25], 'Alpha', 0.78);
Slower.plot(32, 18, 'ShowCP', true, 'ShowIso', true, ...
    'FaceColor', [0.25 0.35 0.90], 'Alpha', 0.78);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Lofted upper and lower wing surfaces');
view(3);

%% ----------------------------------------------------------------------
% Plot 3: structured surface meshes
%% ----------------------------------------------------------------------
fig3 = figure('Name','3 - Structured wing meshes');
hold on; grid on; axis equal;
surf(meshUpper.X, meshUpper.Y, meshUpper.Z, ...
    'FaceColor', [0.90 0.35 0.35], 'FaceAlpha', 0.65, ...
    'EdgeColor', 'k', 'EdgeAlpha', 0.25);
surf(meshLower.X, meshLower.Y, meshLower.Z, ...
    'FaceColor', [0.35 0.45 0.95], 'FaceAlpha', 0.65, ...
    'EdgeColor', 'k', 'EdgeAlpha', 0.25);
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('Structured meshes (nu=%d, nv=%d)', nu, nv));
view(3);

%% ----------------------------------------------------------------------
% Plot 4: TE lines from mesh
%% ----------------------------------------------------------------------
fig4 = figure('Name','4 - Trailing-edge check');
hold on; grid on; axis equal;
plot3(TEu(:,1), TEu(:,2), TEu(:,3), 'r-', 'LineWidth', 2.0);
plot3(TEl(:,1), TEl(:,2), TEl(:,3), 'b--', 'LineWidth', 1.8);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Upper / lower trailing-edge curves from structured mesh');
legend({'Upper TE','Lower TE'}, 'Location','best');
view(3);

%% ----------------------------------------------------------------------
% Exports
%% ----------------------------------------------------------------------
fprintf('\n--- 6. Export meshes ---\n');

if doCSV
    geom.MeshWriter.toCSV(meshUpper, csvUpper);
    geom.MeshWriter.toCSV(meshLower, csvLower);
end

if doVTK
    geom.MeshWriter.toVTK(meshUpper, vtkUpper, 'IncludeNormals', true);
    geom.MeshWriter.toVTK(meshLower, vtkLower, 'IncludeNormals', true);
end

if doSTL
    geom.MeshWriter.toSTL(meshUpper, stlUpper);
    geom.MeshWriter.toSTL(meshLower, stlLower);
end

if doWGS
    netUpper = geom.MeshWriter.meshToNetwork(meshUpper, 'WingUpper', 0, 0);
    netLower = geom.MeshWriter.meshToNetwork(meshLower, 'WingLower', 0, 0);
    geom.MeshWriter.toWGS({netUpper, netLower}, wgsFile);
end

fprintf('\n=== Wing Mesh Export Demo Complete ===\n');
fprintf('Figures:\n');
fprintf('  (1) Root and tip sections in 3D\n');
fprintf('  (2) Lofted upper/lower wing surfaces\n');
fprintf('  (3) Structured surface meshes\n');
fprintf('  (4) Trailing-edge check\n');

%% ======================================================================
% Local functions
%% ======================================================================

function section = buildAirfoilSectionLocal(filename, pFit, fitMethod)
    [~, xy] = readAirfoilFileLocal(filename);

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

    section.full  = Cfull;
    section.upper = Cup_raw.reverse(); % LE -> TE
    section.lower = Clo_raw;           % LE -> TE
    section.LE    = LE;
end

function T = makeWingSectionTransformLocal(chord, alpha, LE_xyz)
    c = cos(alpha);
    s = sin(alpha);

    S = eye(4);
    S(1,1) = chord;
    S(2,2) = chord;
    S(3,3) = chord;

    E = eye(4);
    E(2,2) = 0;
    E(3,2) = 1;
    E(3,3) = 0;

    Ry = eye(4);
    Ry(1,1) =  c;
    Ry(1,3) =  s;
    Ry(3,1) = -s;
    Ry(3,3) =  c;

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