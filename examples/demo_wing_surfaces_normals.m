%% demo_wing_surfaces_normals.m
% Clean demo: build wing surfaces + ruled tip surface
% Visualize surface normals only (no meshing, no seams)

clear; close all; clc;

%% ------------------------------------------------------------------------
% Settings
%% ------------------------------------------------------------------------
rootAirfoilFile = 'airfoil_sharpTE.dat';
tipAirfoilFile  = 'airfoil_sharpTE.dat';

pFit = 3;
fitMethod = 'centripetal';

% Geometry
semispan    = 5.0;
rootChord   = 2.0;
tipChord    = 1.0;

leSweepDeg  = 20.0;
incDeg      = 2.0;
tipTwistDeg = -3.0;
dihedralDeg = 5.0;

%% ------------------------------------------------------------------------
% 1) Airfoil sections
%% ------------------------------------------------------------------------
fprintf('\n--- Airfoil fit ---\n');

rootSection = buildAirfoilSectionLocal(rootAirfoilFile, pFit, fitMethod);
tipSection  = buildAirfoilSectionLocal(tipAirfoilFile,  pFit, fitMethod);

%% ------------------------------------------------------------------------
% 2) Placement transforms
%% ------------------------------------------------------------------------
fprintf('\n--- Section placement ---\n');

sweepRad    = deg2rad(leSweepDeg);
incRad      = deg2rad(incDeg);
tipTwistRad = deg2rad(tipTwistDeg);
dihedRad    = deg2rad(dihedralDeg);

rootLE = [0, 0, 0];
tipLE  = [semispan * tan(sweepRad), semispan, semispan * tan(dihedRad)];

Troot = makeWingSectionTransformLocal(rootChord, incRad, rootLE);
Ttip  = makeWingSectionTransformLocal(tipChord, incRad + tipTwistRad, tipLE);

%% ------------------------------------------------------------------------
% 3) Transform curves
%% ------------------------------------------------------------------------
Cup_root = rootSection.upper.transform(Troot);
Clo_root = rootSection.lower.transform(Troot);

Cup_tip  = tipSection.upper.transform(Ttip);
Clo_tip  = tipSection.lower.transform(Ttip);

%% ------------------------------------------------------------------------
% 4) Build surfaces
%% ------------------------------------------------------------------------
fprintf('\n--- Build surfaces ---\n');

Supper = geom.NURBSSurface.loft({Cup_root, Cup_tip}, 1, 'chord');
Slower = geom.NURBSSurface.loft({Clo_root, Clo_tip}, 1, 'chord');

% Ensure consistent outward normals
[Supper, ~, ~] = orientSurfaceNormalLocal(Supper, 'up', 'WingUpper');
[Slower, ~, ~] = orientSurfaceNormalLocal(Slower, 'down', 'WingLower');

%% ------------------------------------------------------------------------
% 5) Build ruled tip surface (KEY PART)
%% ------------------------------------------------------------------------
fprintf('\n--- Build ruled tip surface ---\n');
curvesC = geom.NURBSSurface.makeCompatibleCurves({Cup_tip, Clo_tip});
Cup_tip = curvesC{1};
Clo_tip = curvesC{2};

Stip = geom.NURBSSurface.ruled(Cup_tip, Clo_tip);

uvmid = [mean(Stip.domainU), mean(Stip.domainV)];
nmid  = Stip.normal(uvmid(1), uvmid(2));

if nmid(2) < 0
    Stip = geom.NURBSSurface.ruled(Clo_tip, Cup_tip);
end

%% ------------------------------------------------------------------------
% 6) Plot surfaces + normals
%% ------------------------------------------------------------------------
figure('Name','Wing surfaces + normals');
hold on; axis equal; grid on; view(3);

plotSurfaceWithNormals(Supper, [0.9 0.3 0.3]);
plotSurfaceWithNormals(Slower, [0.3 0.3 0.9]);
plotSurfaceWithNormals(Stip,   [0.9 0.9 0.3]);

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Wing surfaces with normals');

legend({'Upper','Lower','Tip'});


%% ------------------------------------------------------------------------
% 7) Build meshes using topology system
%% ------------------------------------------------------------------------
fprintf('\n--- Build meshes ---\n');

% Wrap surfaces
U = geom.TopologySurface('upper', Supper);
L = geom.TopologySurface('lower', Slower);
T = geom.TopologySurface('tip',   Stip);

% Build topology
topo = geom.TopologyMesh();
topo.addSurface(U);
topo.addSurface(L);
topo.addSurface(T);

% Declare connections (THIS is the key)
topo.connect(U, 'jN', T, 'i1');   % upper tip edge
topo.connect(L, 'jN', T, 'iN');   % lower tip edge


% Build meshes
res.upper = [25, 10];
res.lower = [25, 10];
res.tip   = [25,4];   % 👈 thin through thickness

topo.buildMeshes(res);

meshUpper = topo.meshes.upper;
meshLower = topo.meshes.lower;
meshTip   = topo.meshes.tip;






%% ------------------------------------------------------------------------
% 8) Plot meshes
%% ------------------------------------------------------------------------
figure('Name','Wing meshes');
hold on; axis equal; grid on; view(3);

plotMesh(meshUpper, [1 0 0]);
plotMesh(meshLower, [0 0 1]);
plotMesh(meshTip,   [1 1 0]);

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Structured mesh');


figure('Name','TIP DEBUG');
clf; hold on; axis equal; grid on; view(3);

plotMesh(meshTip, [0 0 0]);   % wireframe

title('Tip Mesh Debug');
xlabel('X'); ylabel('Y'); zlabel('Z');




















%% ========================================================================
% Helpers (minimal)
%% ========================================================================


function plotSurfaceWithNormals(S, color)
    nu = 15;
    nv = 15;

    u = linspace(S.domainU(1), S.domainU(2), nu);
    v = linspace(S.domainV(1), S.domainV(2), nv);

    X = zeros(nu,nv);
    Y = zeros(nu,nv);
    Z = zeros(nu,nv);

    for i=1:nu
        for j=1:nv
            P = S.evaluate(u(i), v(j));
            X(i,j)=P(1); Y(i,j)=P(2); Z(i,j)=P(3);
        end
    end

    surf(X,Y,Z,'FaceColor',color,'FaceAlpha',0.6,'EdgeColor','none');

    % normals
    skip = 3;
    for i=1:skip:nu
        for j=1:skip:nv
            SKL = S.derivatives(u(i), v(j), 1);
            Su = SKL{2,1};
            Sv = SKL{1,2};
            n = cross(Su,Sv);
            n = n / norm(n);

            quiver3(X(i,j),Y(i,j),Z(i,j), ...
                n(1),n(2),n(3),0.5,'k');
        end
    end
end

function section = buildAirfoilSectionLocal(filename, pFit, fitMethod)
    [~, xy] = readAirfoilFileLocal(filename);
    xy3 = [xy, zeros(size(xy,1),1)];
    Cfull = geom.NURBSCurve.globalInterp(xy3, pFit, fitMethod);

    [uLE, pLE3] = findLeadingEdgeBySlopeLocal(Cfull, 2000); 
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
    end

    if isempty(data)
        error('No numeric x-y coordinate pairs found in file.');
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
    cand(end+1) = us(idxMin);
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


function T = makeWingSectionTransformLocal(chord, alpha, LE_xyz)
    c = cos(alpha);
    s = sin(alpha);

    S = eye(4);
    S(1,1) = chord;
    S(2,2) = chord;
    S(3,3) = chord;

    % Embed airfoil x-y into global X-Z at constant span station Y
    E = eye(4);
    E(2,2) = 0;
    E(3,2) = 1;
    E(3,3) = 0;

    Ry = eye(4);
    Ry(1,1) = c;
    Ry(1,3) = s;
    Ry(3,1) = -s;
    Ry(3,3) = c;

    Tr = eye(4);
    Tr(1:3,4) = LE_xyz(:);

    T = Tr * Ry * E * S;
end


function [Sout, labelOut, nOut] = orientSurfaceNormalLocal(Sin, expectedNormal, label)
    variants = {Sin, Sin.flipNormals()};
    labels   = {'none', 'flipNormals'};

    bestIdx = 1;
    bestPenalty = inf;
    for k = 1:numel(variants)
        n = sampleSurfaceNormalLocal(variants{k});
        p = normalPenaltyLocal(n, expectedNormal);
        if p < bestPenalty
            bestPenalty = p;
            bestIdx = k;
        end
    end

    Sout = variants{bestIdx};
    labelOut = labels{bestIdx};
    nOut = sampleSurfaceNormalLocal(Sout);

    if normalPenaltyLocal(nOut, expectedNormal) > 0
        error('orientSurfaceNormalLocal: could not satisfy requested normal for %s.', label);
    end
end

function n = sampleSurfaceNormalLocal(S)
    % Sample near the middle of the parametric domain and use first partials.
    udom = S.domainU;
    vdom = S.domainV;

    u = 0.5 * (udom(1) + udom(2));
    v = 0.5 * (vdom(1) + vdom(2));

    % Get derivative tensor up to first order.
    SKL = S.derivatives(u, v, 1);

    % SKL{k+1,l+1} corresponds to d^(k+l)S / du^k dv^l
    Su = SKL{2,1};   % dS/du
    Sv = SKL{1,2};   % dS/dv

    n = cross(Su, Sv);
    nn = norm(n);

    if nn < 1e-14
        % fallback: try a slightly perturbed point
        du = 1e-3 * max(udom(2) - udom(1), 1);
        dv = 1e-3 * max(vdom(2) - vdom(1), 1);
        u2 = min(max(u + du, udom(1)), udom(2));
        v2 = min(max(v + dv, vdom(1)), vdom(2));

        SKL = S.derivatives(u2, v2, 1);
        Su = SKL{2,1};
        Sv = SKL{1,2};

        n = cross(Su, Sv);
        nn = norm(n);
    end

    if nn < 1e-14
        error('sampleSurfaceNormalLocal: could not compute a reliable surface normal.');
    end

    n = n / nn;
end


function p = normalPenaltyLocal(n, expectedNormal)
    switch lower(strtrim(expectedNormal))
        case 'up'
            p = double(n(3) < 0);
        case 'down'
            p = double(n(3) > 0);
        case 'outboard'
            p = double(n(2) < 0);
        otherwise
            error('Unknown expected normal direction: %s', expectedNormal);
    end
end


function plotMesh(mesh, color)

    surf(mesh.X, mesh.Y, mesh.Z, ...
        'FaceColor', 'none', ...
        'EdgeColor', color, ...
        'LineWidth', 1.0);

end


