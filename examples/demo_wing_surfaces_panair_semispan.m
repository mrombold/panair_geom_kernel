%% demo_wing_surfaces_panair_semispan.m
% Surface-based semispan wing demo for PANAIR export.
% Uses:
%   - lofted upper/lower wing surfaces
%   - ruled tip surface from compatible tip curves
%   - direct meshing from surfaces
%   - TE averaging + planar wake
%   - WGS + AUX export
%
% Run from repo root:
%   demo_wing_surfaces_panair_semispan

clear; close all; clc;

%% ------------------------------------------------------------------------
% User settings
%% ------------------------------------------------------------------------
rootAirfoilFile = 'airfoil_sharpTE.dat';
tipAirfoilFile  = 'airfoil_sharpTE.dat';
pFit = 3;
fitMethod = 'centripetal';

% Semispan wing geometry
semispan    = 5.0;
rootChord   = 2.0;
tipChord    = 1.0;
leSweepDeg  = 20.0;
incDeg      = 2.0;
tipTwistDeg = -3.0;
dihedralDeg = 5.0;
loftDegree  = 1;

% Mesh density
nuWing = 41;    % chordwise
nvWing = 20;    % spanwise
nvTip  = 5;     % thickness-wise on tip closeout
nuWake = 2;
nvWake = nvWing;

spacingUWing = 'cosine';
spacingVWing = 'linear';
spacingUTip  = spacingUWing;   % must match wing chordwise spacing
spacingVTip  = 'linear';

% Wake
wakeLength = 50 * rootChord;
wakeXDir   = [1, 0, 0];

% Analysis conditions for AUX file
mach  = 0.20;
alpha = 2.0;
cbar  = 0.5 * (rootChord + tipChord);
span  = 2.0 * semispan;
sref  = 2.0 * semispan * 0.5 * (rootChord + tipChord);
xref  = 0.25 * cbar;
zref  = 0.0;

% Output
outDir  = 'exports';
base    = 'wing_surfaces_semispan_refined';
wgsFile = fullfile(outDir, [base '.wgs']);
auxFile = fullfile(outDir, [base '.aux']);

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
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% ------------------------------------------------------------------------
% 1) Read and fit root / tip airfoils
%% ------------------------------------------------------------------------
fprintf('\n--- 1. Read and fit root / tip airfoils ---\n');
rootSection = buildAirfoilSectionLocal(rootAirfoilFile, pFit, fitMethod);
tipSection  = buildAirfoilSectionLocal(tipAirfoilFile,  pFit, fitMethod);
fprintf(' root file: %s\n', rootAirfoilFile);
fprintf(' tip file: %s\n',  tipAirfoilFile);
fprintf(' root fitted LE = [%.8f %.8f]\n', rootSection.LE(1), rootSection.LE(2));
fprintf(' tip fitted LE = [%.8f %.8f]\n',  tipSection.LE(1),  tipSection.LE(2));

%% ------------------------------------------------------------------------
% 2) Section placement
%% ------------------------------------------------------------------------
fprintf('\n--- 2. Build semispan section placement transforms ---\n');
sweepRad    = deg2rad(leSweepDeg);
incRad      = deg2rad(incDeg);
tipTwistRad = deg2rad(tipTwistDeg);
dihedRad    = deg2rad(dihedralDeg);

rootLE = [0, 0, 0];
tipLE  = [semispan * tan(sweepRad), semispan, semispan * tan(dihedRad)];
fprintf(' root LE = [%.6f %.6f %.6f]\n', rootLE);
fprintf(' tip LE = [%.6f %.6f %.6f]\n', tipLE);

Troot = makeWingSectionTransformLocal(rootChord, incRad, rootLE);
Ttip  = makeWingSectionTransformLocal(tipChord, incRad + tipTwistRad, tipLE);

%% ------------------------------------------------------------------------
% 3) Transform section curves into 3D
%% ------------------------------------------------------------------------
fprintf('\n--- 3. Transform section curves into 3D ---\n');
Cup_root = rootSection.upper.transform(Troot);
Clo_root = rootSection.lower.transform(Troot);
Cup_tip  = tipSection.upper.transform(Ttip);
Clo_tip  = tipSection.lower.transform(Ttip);

%% ------------------------------------------------------------------------
% 4) Build upper / lower wing surfaces
%% ------------------------------------------------------------------------
fprintf('\n--- 4. Build upper / lower wing surfaces ---\n');
Supper = geom.NURBSSurface.loft({Cup_root, Cup_tip}, loftDegree, 'chord');
Slower = geom.NURBSSurface.loft({Clo_root, Clo_tip}, loftDegree, 'chord');

[Supper, upperSurfLabel, nSU] = orientSurfaceNormalLocal(Supper, 'up',   'WingUpper');
[Slower, lowerSurfLabel, nSL] = orientSurfaceNormalLocal(Slower, 'down', 'WingLower');

fprintf(' surface normals: WingUpper=%s [%.3f %.3f %.3f], WingLower=%s [%.3f %.3f %.3f]\n', ...
    upperSurfLabel, nSU(1), nSU(2), nSU(3), ...
    lowerSurfLabel, nSL(1), nSL(2), nSL(3));

%% ------------------------------------------------------------------------
% 5) Build ruled wing-tip closeout surface
%% ------------------------------------------------------------------------
fprintf('\n--- 5. Build ruled tip surface ---\n');
curvesC = geom.NURBSSurface.makeCompatibleCurves({Cup_tip, Clo_tip});
Cup_tip_c = curvesC{1};
Clo_tip_c = curvesC{2};

Stip = geom.NURBSSurface.ruled(Cup_tip_c, Clo_tip_c);

uvmid = [mean(Stip.domainU), mean(Stip.domainV)];
nmid  = Stip.normal(uvmid(1), uvmid(2));
fprintf(' midpoint tip normal before orientation check = [%.6f %.6f %.6f]\n', ...
    nmid(1), nmid(2), nmid(3));

if nmid(2) < 0
    Stip = geom.NURBSSurface.ruled(Clo_tip_c, Cup_tip_c);
    fprintf(' tip closeout flipped so midpoint normal points roughly +Y.\n');
end

%% ------------------------------------------------------------------------
% 6) Mesh surfaces directly from geometry
%% ------------------------------------------------------------------------
fprintf('\n--- 6. Mesh surfaces ---\n');
meshUpper = Supper.isoMesh(nuWing, nvWing, ...
    'SpacingU', spacingUWing, 'SpacingV', spacingVWing);
meshLower = Slower.isoMesh(nuWing, nvWing, ...
    'SpacingU', spacingUWing, 'SpacingV', spacingVWing);

% For ruled tip surface:
%   U = chordwise along the airfoil branches
%   V = thickness-wise between upper/lower branches
meshTip = Stip.isoMesh(nuWing, nvTip, ...
    'SpacingU', spacingUTip, ...
    'SpacingV', spacingVTip);

fprintf(' upper mesh quads = %d\n', size(meshUpper.connectivity,1));
fprintf(' lower mesh quads = %d\n', size(meshLower.connectivity,1));
fprintf(' tip   mesh quads = %d\n', size(meshTip.connectivity,1));

%% ------------------------------------------------------------------------
% 7) Validate shared tip edges (non-destructive)
%% ------------------------------------------------------------------------
fprintf('\n--- 7. Validate shared tip joins ---\n');
Pup = getMeshEdgeLocal(meshUpper, 'jN');
Plo = getMeshEdgeLocal(meshLower, 'jN');

% In the working ruled-tip demo, the upper/lower airfoil branches appear on
% the tip mesh j1 / jN edges respectively.
Tj1 = getMeshEdgeLocal(meshTip, 'j1');
TjN = getMeshEdgeLocal(meshTip, 'jN');

cmpUp = compareEdgeOrientationLocal(Pup, Tj1);
cmpLo = compareEdgeOrientationLocal(Plo, TjN);

fprintf(' WingUpper jN <-> WingTip j1 : err = %.3e (%s)\n', ...
    cmpUp.bestError, cmpUp.orientation);
fprintf(' WingLower jN <-> WingTip jN : err = %.3e (%s)\n', ...
    cmpLo.bestError, cmpLo.orientation);

%% ------------------------------------------------------------------------
% 8) Identify LE / TE roles on wing surfaces
%% ------------------------------------------------------------------------
fprintf('\n--- 8. Detect LE / TE roles ---\n');
[teEdgeUpper, leEdgeUpper] = detectChordwiseEdgeRolesLocal(meshUpper);
[teEdgeLower, leEdgeLower] = detectChordwiseEdgeRolesLocal(meshLower);

fprintf(' detected chordwise roles: WingUpper TE=%s LE=%s; WingLower TE=%s LE=%s\n', ...
    teEdgeUpper, leEdgeUpper, teEdgeLower, leEdgeLower);

%% ------------------------------------------------------------------------
% 9) Build semispan wake from unified trailing edge
%% ------------------------------------------------------------------------
fprintf('\n--- 9. Build semispan wake network ---\n');
TEu = getMeshEdgeLocal(meshUpper, teEdgeUpper);
TEl_raw = getMeshEdgeLocal(meshLower, teEdgeLower);

teCmp = compareEdgeOrientationLocal(TEu, TEl_raw);
TEl   = orientEdgeLikeLocal(TEl_raw, teCmp.orientation);

dTE = vecnorm(TEu - TEl, 2, 2);
fprintf(' upper/lower TE mismatch before averaging: max = %.3e, RMS = %.3e\n', ...
    max(dTE), sqrt(mean(dTE.^2)));

TE  = 0.5 * (TEu + TEl);

% Intentionally unify upper/lower trailing-edge coordinates for wake attach
meshUpper = setMeshEdgeLocal(meshUpper, teEdgeUpper, TE);
meshLower = setMeshEdgeLocal(meshLower, teEdgeLower, TE);

wakeDir = wakeXDir(:).' / norm(wakeXDir);
TE2     = TE + wakeLength * wakeDir;

meshWake = struct();
meshWake.X = [TE(:,1), TE2(:,1)].';
meshWake.Y = [TE(:,2), TE2(:,2)].';
meshWake.Z = [TE(:,3), TE2(:,3)].';
meshWake.nu = size(meshWake.X,1);
meshWake.nv = size(meshWake.X,2);
meshWake.connectivity = buildConnectivityLocal(meshWake.nu, meshWake.nv);
meshWake.normals = [];

% Match the old semispan demo orientation/layout
meshWake = transposeStructuredMeshLocal(meshWake);
fprintf(' wake mesh quads = %d\n', size(meshWake.connectivity,1));

%% ------------------------------------------------------------------------
% 10) Final normals
%% ------------------------------------------------------------------------
fprintf('\n--- 10. Final network normals ---\n');
nU = geom.MeshWriter.averageNetworkNormal(meshUpper);
nL = geom.MeshWriter.averageNetworkNormal(meshLower);
nT = geom.MeshWriter.averageNetworkNormal(meshTip);

fprintf(' WingUpper normal = [%.3f %.3f %.3f]\n', nU);
fprintf(' WingLower normal = [%.3f %.3f %.3f]\n', nL);
fprintf(' WingTip   normal = [%.3f %.3f %.3f]\n', nT);

%% ------------------------------------------------------------------------
% 11) Plot geometry and meshes
%% ------------------------------------------------------------------------
figure('Name','Semispan wing surfaces for PANAIR');
hold on; grid on; axis equal; view(3);

surf(meshUpper.X, meshUpper.Y, meshUpper.Z, ...
    'FaceColor', [0.90 0.35 0.35], 'FaceAlpha', 0.65, ...
    'EdgeColor', 'k', 'EdgeAlpha', 0.20);

surf(meshLower.X, meshLower.Y, meshLower.Z, ...
    'FaceColor', [0.35 0.45 0.95], 'FaceAlpha', 0.65, ...
    'EdgeColor', 'k', 'EdgeAlpha', 0.20);

surf(meshTip.X, meshTip.Y, meshTip.Z, ...
    'FaceColor', [0.95 0.90 0.35], 'FaceAlpha', 0.78, ...
    'EdgeColor', 'k', 'EdgeAlpha', 0.30);

surf(meshWake.X, meshWake.Y, meshWake.Z, ...
    'FaceColor', [0.45 0.95 0.55], 'FaceAlpha', 0.35, ...
    'EdgeColor', 'k', 'EdgeAlpha', 0.25);

plot3([0 0], [0 semispan], [0 0], 'k--', 'LineWidth', 1.0);

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Semispan wing networks for PANAIR-style export');
legend({'Upper','Lower','Tip','Wake','Symmetry line'}, 'Location','best');

%% ------------------------------------------------------------------------
% 12) Export LaWGS and AUX
%% ------------------------------------------------------------------------
fprintf('\n--- 12. Export LaWGS and AUX ---\n');
netUpper = geom.MeshWriter.meshToNetwork(meshUpper, 'WingUpper', 0, 0);
netLower = geom.MeshWriter.meshToNetwork(meshLower, 'WingLower', 0, 0);
netTip   = geom.MeshWriter.meshToNetwork(meshTip,   'WingTip',   0, 0);
netWake  = geom.MeshWriter.meshToNetwork(meshWake,  'WingWake',  0, 0);

nets = {netUpper, netLower, netTip, netWake};

fprintf(' Running network edge audit before WGS export...\n');
report = geom.MeshWriter.auditWGSNetworks(nets, 'Tolerance', 1e-10, 'Verbose', true); %#ok<NASGU>

% Use the ruled-tip edge mapping:
%   WingUpper jN <-> WingTip j1
%   WingLower jN <-> WingTip jN
%
% Leave LE/TE joins as bookkeeping concepts; only assert what matters for export.
topologyRules = {
    'WingUpper', 'jN',       'regular',   'TipCloseout', 'WingTip',  'j1';
    'WingLower', 'jN',       'regular',   'TipCloseout', 'WingTip',  'jN';
    'WingUpper', teEdgeUpper,'sharp_te',  'TrailingEdge','',         '';
    'WingLower', teEdgeLower,'sharp_te',  'TrailingEdge','',         '';
    'WingWake',  'j1',       'sharp_te',  'TrailingEdge','',         ''};

try
    geom.MeshWriter.assertTopologyRules(nets, topologyRules, ...
        'Tolerance', 1e-10, 'Verbose', true);
catch ME
    fprintf(' Topology rule check warning: %s\n', ME.message);
end

geom.MeshWriter.toWGS(nets, wgsFile);
writePaninAuxLocal(auxFile, wgsFile, nets, ...
    mach, alpha, cbar, span, sref, xref, zref);

fprintf(' wrote WGS: %s\n', wgsFile);
fprintf(' wrote AUX: %s\n', auxFile);
fprintf('\nDone.\n');
fprintf('Next step:\n');
fprintf(' Place panin and these files together, run panin, and enter:\n');
fprintf(' %s\n', [base '.aux']);

%% ========================================================================
% Local functions
%% ========================================================================
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

function conn = buildConnectivityLocal(nu, nv)
    nquads = (nu - 1) * (nv - 1);
    conn = zeros(nquads, 4);
    q = 1;
    for j = 1:(nv - 1)
        for i = 1:(nu - 1)
            n1 = sub2ind([nu, nv], i,   j);
            n2 = sub2ind([nu, nv], i+1, j);
            n3 = sub2ind([nu, nv], i+1, j+1);
            n4 = sub2ind([nu, nv], i,   j+1);
            conn(q,:) = [n1 n2 n3 n4];
            q = q + 1;
        end
    end
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
    udom = S.domainU;
    vdom = S.domainV;

    u = 0.5 * (udom(1) + udom(2));
    v = 0.5 * (vdom(1) + vdom(2));

    SKL = S.derivatives(u, v, 1);
    Su = SKL{2,1};
    Sv = SKL{1,2};

    n = cross(Su, Sv);
    nn = norm(n);

    if nn < 1e-14
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

function [teEdge, leEdge] = detectChordwiseEdgeRolesLocal(mesh)
    e1 = getMeshEdgeLocal(mesh, 'i1');
    eN = getMeshEdgeLocal(mesh, 'iN');

    x1 = mean(e1(:,1));
    xN = mean(eN(:,1));

    if xN >= x1
        teEdge = 'iN';
        leEdge = 'i1';
    else
        teEdge = 'i1';
        leEdge = 'iN';
    end
end

function meshT = transposeStructuredMeshLocal(mesh)
    meshT = mesh;
    meshT.X = mesh.X.';
    meshT.Y = mesh.Y.';
    meshT.Z = mesh.Z.';
    meshT.nu = size(meshT.X, 1);
    meshT.nv = size(meshT.X, 2);
    meshT.connectivity = buildConnectivityLocal(meshT.nu, meshT.nv);
    meshT.normals = [];
end

function P = getMeshEdgeLocal(mesh, edgeName)
    switch lower(strtrim(edgeName))
        case {'i1','imin','row1'}
            P = [mesh.X(1,:).',   mesh.Y(1,:).',   mesh.Z(1,:).'];
        case {'in','imax','rowend'}
            P = [mesh.X(end,:).', mesh.Y(end,:).', mesh.Z(end,:).'];
        case {'j1','jmin','col1'}
            P = [mesh.X(:,1),     mesh.Y(:,1),     mesh.Z(:,1)];
        case {'jn','jmax','colend'}
            P = [mesh.X(:,end),   mesh.Y(:,end),   mesh.Z(:,end)];
        otherwise
            error('Unknown mesh edge name: %s', edgeName);
    end
end

function mesh = setMeshEdgeLocal(mesh, edgeName, P)
    if size(P,2) ~= 3
        error('Edge array must be n-by-3.');
    end
    switch lower(strtrim(edgeName))
        case {'i1','imin','row1'}
            if size(P,1) ~= mesh.nv
                error('setMeshEdgeLocal: edge length mismatch for %s.', edgeName);
            end
            mesh.X(1,:) = P(:,1).';
            mesh.Y(1,:) = P(:,2).';
            mesh.Z(1,:) = P(:,3).';
        case {'in','imax','rowend'}
            if size(P,1) ~= mesh.nv
                error('setMeshEdgeLocal: edge length mismatch for %s.', edgeName);
            end
            mesh.X(end,:) = P(:,1).';
            mesh.Y(end,:) = P(:,2).';
            mesh.Z(end,:) = P(:,3).';
        case {'j1','jmin','col1'}
            if size(P,1) ~= mesh.nu
                error('setMeshEdgeLocal: edge length mismatch for %s.', edgeName);
            end
            mesh.X(:,1) = P(:,1);
            mesh.Y(:,1) = P(:,2);
            mesh.Z(:,1) = P(:,3);
        case {'jn','jmax','colend'}
            if size(P,1) ~= mesh.nu
                error('setMeshEdgeLocal: edge length mismatch for %s.', edgeName);
            end
            mesh.X(:,end) = P(:,1);
            mesh.Y(:,end) = P(:,2);
            mesh.Z(:,end) = P(:,3);
        otherwise
            error('Unknown mesh edge name: %s', edgeName);
    end
end

function writePaninAuxLocal(filename, wgsFile, networks, mach, alpha, cbar, span, sref, xref, zref)
    [~, name, ext] = fileparts(wgsFile);
    wgsBase = [name ext];

    if isnumeric(alpha) && isscalar(alpha)
        alphaStr = sprintf('%.8g', alpha);
    elseif isnumeric(alpha)
        alphaStr = sprintf(' %.8g', alpha);
        alphaStr = strtrim(alphaStr);
    else
        error('alpha must be numeric.');
    end

    if ~iscell(networks)
        tmp = networks;
        networks = cell(1, numel(tmp));
        for k = 1:numel(tmp)
            networks{k} = tmp(k);
        end
    end

    bounVals = zeros(1, numel(networks));
    for k = 1:numel(networks)
        net = networks{k};
        if isfield(net, 'boun') && ~isempty(net.boun)
            bounVals(k) = net.boun;
            continue;
        end

        netName = '';
        if isfield(net, 'name') && ~isempty(net.name)
            netName = char(net.name);
        end
        if contains(upper(strtrim(netName)), 'WAKE')
            bounVals(k) = 18;
        else
            bounVals(k) = 1;
        end
    end
    bounStr = strtrim(sprintf('%d ', bounVals));

    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not open AUX file for writing: %s', filename);
    end
    c = onCleanup(@() fclose(fid));

    fprintf(fid, '// Auxiliary file for %s\n', upper(wgsBase));
    fprintf(fid, 'WGS %s\n', upper(wgsBase));
    fprintf(fid, 'MACH %.8g\n', mach);
    fprintf(fid, 'ALPHA %s\n', alphaStr);
    fprintf(fid, 'CBAR %.8g\n', cbar);
    fprintf(fid, 'SPAN %.8g\n', span);
    fprintf(fid, 'SREF %.8g\n', sref);
    fprintf(fid, 'XREF %.8g\n', xref);
    fprintf(fid, 'ZREF %.8g\n', zref);
    fprintf(fid, 'BOUN %s\n', bounStr);
end

function cmp = compareEdgeOrientationLocal(P, Q)
    if size(P,2) ~= 3 || size(Q,2) ~= 3
        error('compareEdgeOrientationLocal: P and Q must be n-by-3.');
    end
    if size(P,1) ~= size(Q,1)
        error('compareEdgeOrientationLocal: edge lengths must match.');
    end

    dSame = vecnorm(P - Q, 2, 2);
    dRev  = vecnorm(P - flipud(Q), 2, 2);

    cmp.sameError = max(dSame);
    cmp.revError  = max(dRev);

    if cmp.revError < cmp.sameError
        cmp.orientation = 'reversed';
        cmp.bestError = cmp.revError;
    else
        cmp.orientation = 'same';
        cmp.bestError = cmp.sameError;
    end
end

function Qout = orientEdgeLikeLocal(Qin, orientation)
    switch lower(strtrim(orientation))
        case 'same'
            Qout = Qin;
        case 'reversed'
            Qout = flipud(Qin);
        otherwise
            error('orientEdgeLikeLocal: orientation must be ''same'' or ''reversed''.');
    end
end