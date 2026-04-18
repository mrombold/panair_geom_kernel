%% demo_tip_closeout_same_as_repo.m
% Build and mesh ONLY the wing-tip closeout surface using the SAME airfoil
% read/fit/split logic already used in the repo demos.
%
% This does NOT modify +geom/NURBSCurve.m.
%
% Workflow matches the current repo:
%   read Selig-style airfoil file
%   -> fit one continuous NURBS curve
%   -> find geometric LE from fitted curve by dx/du = 0
%   -> split fitted curve at LE
%   -> reverse upper branch so both branches go LE -> TE
%   -> transform to 3D tip section
%   -> build ruled closeout surface
%   -> mesh and plot only that surface
%
% Run from the repo root:
%   demo_tip_closeout_same_as_repo

clear; close all; clc;

%% ------------------------------------------------------------------------
% User settings
%% ------------------------------------------------------------------------
airfoilFile  = 'airfoil_sharpTE.dat';   % same style as repo demo
pFit         = 3;
fitMethod    = 'centripetal';
tipLE        = [0.0, 5.0, 0.0];
tipChord     = 1.0;
tipTwistDeg  = 0.0;
dihedralDeg  = 0.0;

nuMesh       = 41;
nvMesh       = 5;
nPlotCurve   = 500;

%% ------------------------------------------------------------------------
% Path / package check
%% ------------------------------------------------------------------------
script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir)
    addpath(script_dir);
end
if ~exist('+geom', 'dir') && isempty(which('geom.NURBSCurve'))
    error(['Cannot find +geom package.\n' ...
           'Run this from the repo root or add the geometry kernel to path.']);
end
fprintf('geom package found at: %s\n', script_dir);

%% ------------------------------------------------------------------------
% Read and fit airfoil using same helper pattern as repo demos
%% ------------------------------------------------------------------------
fprintf('\n--- 1. Read and fit airfoil ---\n');
section = buildAirfoilSectionLocal(airfoilFile, pFit, fitMethod);
fprintf(' airfoil file: %s\n', airfoilFile);
fprintf(' fitted LE = [%.8f %.8f]\n', section.LE(1), section.LE(2));

%% ------------------------------------------------------------------------
% Transform split branches into 3D tip section
%% ------------------------------------------------------------------------
fprintf('\n--- 2. Transform section to tip station ---\n');
Cup_tip = transformSectionCurveLocal(section.Cup, tipLE, tipChord, tipTwistDeg, dihedralDeg);
Clo_tip = transformSectionCurveLocal(section.Clo, tipLE, tipChord, tipTwistDeg, dihedralDeg);

%% ------------------------------------------------------------------------
% Create only the tip closeout surface
%% ------------------------------------------------------------------------
fprintf('\n--- 3. Build tip closeout surface ---\n');
curvesC = geom.NURBSSurface.makeCompatibleCurves({Cup_tip, Clo_tip});
Cup_tip = curvesC{1};
Clo_tip = curvesC{2};
Stip = geom.NURBSSurface.ruled(Cup_tip, Clo_tip);

uvmid = [mean(Stip.domainU), mean(Stip.domainV)];
nmid  = Stip.normal(uvmid(1), uvmid(2));
fprintf(' midpoint normal before orientation check = [%.6f %.6f %.6f]\n', ...
    nmid(1), nmid(2), nmid(3));

if nmid(2) < 0
    Stip = geom.NURBSSurface.ruled(Clo_tip, Cup_tip);
    fprintf(' tip closeout flipped so midpoint normal points roughly +Y.\n');
end

%% ------------------------------------------------------------------------
% Mesh only the tip closeout
%% ------------------------------------------------------------------------
fprintf('\n--- 4. Mesh tip closeout surface ---\n');
meshTip = Stip.isoMesh(nuMesh, nvMesh, ...
    'SpacingU', 'cosine', ...
    'SpacingV', 'linear');

fprintf(' nu = %d, nv = %d\n', nuMesh, nvMesh);
fprintf(' nodes = %d\n', numel(meshTip.X));
fprintf(' quads = %d\n', size(meshTip.connectivity,1));

%% ------------------------------------------------------------------------
% Plot curves used to make the closeout
%% ------------------------------------------------------------------------
uUp = linspace(section.Cup.domain(1), section.Cup.domain(2), nPlotCurve);
uLo = linspace(section.Clo.domain(1), section.Clo.domain(2), nPlotCurve);
xyUp3 = section.Cup.evaluate(uUp);
xyLo3 = section.Clo.evaluate(uLo);

figure('Name','1 - Split airfoil branches used for tip');
hold on; grid on; axis equal;
plot(section.xy(:,1), section.xy(:,2), 'k.', 'MarkerSize', 8);
plot(xyUp3(:,1), xyUp3(:,2), 'r-', 'LineWidth', 2.0);
plot(xyLo3(:,1), xyLo3(:,2), 'b-', 'LineWidth', 2.0);
plot(section.LE(1), section.LE(2), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
xlabel('x'); ylabel('y');
title('Split fitted airfoil branches (same repo logic)');
legend({'Raw points','Upper fit','Lower fit','Geometric LE'}, 'Location','best');

%% ------------------------------------------------------------------------
% Plot only the tip closeout mesh
%% ------------------------------------------------------------------------
figure('Name','2 - Tip closeout mesh only');
hold on; grid on; axis equal; view(3);

surf(meshTip.X, meshTip.Y, meshTip.Z, zeros(size(meshTip.X)), ...
    'FaceColor', [0.85 0.90 1.00], ...
    'FaceAlpha', 0.70, ...
    'EdgeColor', [0.15 0.25 0.55], ...
    'LineWidth', 0.8);

plot3(meshTip.X(:,1),   meshTip.Y(:,1),   meshTip.Z(:,1),   'r-', 'LineWidth', 2.0);
plot3(meshTip.X(:,end), meshTip.Y(:,end), meshTip.Z(:,end), 'b-', 'LineWidth', 2.0);
plot3(meshTip.X(1,:),   meshTip.Y(1,:),   meshTip.Z(1,:),   'k-', 'LineWidth', 2.0);
plot3(meshTip.X(end,:), meshTip.Y(end,:), meshTip.Z(end,:), 'k-', 'LineWidth', 2.0);

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Wing-tip closeout surface mesh only');
camlight headlight; lighting gouraud;

skipU = max(1, floor(nuMesh/10));
skipV = max(1, floor(nvMesh/8));
quiver3(meshTip.X(1:skipU:end,1:skipV:end), ...
        meshTip.Y(1:skipU:end,1:skipV:end), ...
        meshTip.Z(1:skipU:end,1:skipV:end), ...
        0.06 * meshTip.normals(1:skipU:end,1:skipV:end,1), ...
        0.06 * meshTip.normals(1:skipU:end,1:skipV:end,2), ...
        0.06 * meshTip.normals(1:skipU:end,1:skipV:end,3), ...
        0, 'Color', [0.0 0.45 0.0]);

%% ------------------------------------------------------------------------
% Optional export
%% ------------------------------------------------------------------------
if exist('geom.MeshWriter', 'class') || ~isempty(which('geom.MeshWriter'))
    if ~exist('exports', 'dir'), mkdir('exports'); end
    try
        geom.MeshWriter.toVTK({meshTip}, 'exports/tip_closeout_same_as_repo.vtk');
        geom.MeshWriter.toSTL({meshTip}, 'exports/tip_closeout_same_as_repo.stl');
        fprintf(' exported VTK/STL for tip closeout only.\n');
    catch ME
        fprintf(' export skipped/failed: %s\n', ME.message);
    end
end

fprintf('\nDone.\n');

%% ========================================================================
% Local functions copied in the same style as the repo demos
%% ========================================================================
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
    uUp   = linspace(Cup.domain(1),   Cup.domain(2),   600);
    uLo   = linspace(Clo.domain(1),   Clo.domain(2),   600);

    xyFull3 = Cfull.evaluate(uFull);
    xyUp3   = Cup.evaluate(uUp);
    xyLo3   = Clo.evaluate(uLo);

    section.title = titleLine;
    section.xy    = xy;
    section.Cfull = Cfull;
    section.Cup   = Cup;
    section.Clo   = Clo;
    section.LE    = LE;
    section.xyFit = xyFull3(:,1:2);
    section.xyUp  = xyUp3(:,1:2);
    section.xyLo  = xyLo3(:,1:2);
end

function C3 = transformSectionCurveLocal(C2, LE3, chord, twistDeg, dihedralDeg)
    P2 = C2.P;

    % Same 2D->3D convention used in repo wing demos:
    % local airfoil (x,y) -> global (X,Z) at constant span station Y
    X = chord * P2(:,1);
    Y = zeros(size(X));
    Z = chord * P2(:,2);

    th = deg2rad(twistDeg);
    Xt =  cos(th)*X + sin(th)*Z;
    Yt =  Y;
    Zt = -sin(th)*X + cos(th)*Z;

    ph = deg2rad(dihedralDeg);
    Xd = Xt;
    Yd =  cos(ph)*Yt - sin(ph)*Zt;
    Zd =  sin(ph)*Yt + cos(ph)*Zt;

    P3 = [Xd + LE3(1), Yd + LE3(2), Zd + LE3(3)];
    C3 = geom.NURBSCurve(P3, C2.p, C2.U, C2.W);
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
        if ~ischar(t), break; end
        lines{end+1,1} = t; %#ok<AGROW>
    end

    if isempty(lines)
        error('Airfoil file is empty.');
    end

    titleLine = lines{1};
    data = [];

    for i = 2:numel(lines)
        s = strtrim(lines{i});
        if isempty(s), continue; end
        vals = sscanf(s, '%f %f');
        if numel(vals) >= 2
            data(end+1,:) = vals(1:2).'; %#ok<AGROW>
        end
    end

    if isempty(data)
        for i = 1:numel(lines)
            s = strtrim(lines{i});
            if isempty(s), continue; end
            vals = sscanf(s, '%f %f');
            if numel(vals) >= 2
                data(end+1,:) = vals(1:2).'; %#ok<AGROW>
            end
        end
    end

    if isempty(data)
        error('No numeric x-y coordinate pairs found in file.');
    end

    xy = data;
end

function [uLE, pLE, dxduLE] = findLeadingEdgeBySlopeLocal(C, nSearch)
    us   = linspace(C.domain(1), C.domain(2), nSearch);
    pts  = C.evaluate(us);
    d1   = C.derivative(us, 1);
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
        dLE = C.derivative(uLE,1);
        dxduLE = dLE(1);
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
