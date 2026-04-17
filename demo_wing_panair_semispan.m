%% demo_wing_panair_semispan.m
% Build a semispan wing on the XZ symmetry plane, create upper/lower/tip/wake
% structured networks, then write:
% - a LaWGS geometry file (.wgs)
% - a simple panin auxiliary file (.aux)
%
% This is the next step toward a pyPanair / panin / panair workflow.
%
% Notes:
% - Root section is placed on y = 0 (the XZ symmetry plane)
% - Only the starboard semispan is meshed/exported
% - A simple planar wake is shed from the trailing edge
% - Boundary type 1 is used for solid-wall networks
% - Boundary type 18 is used for the wake network
%
% Run from the repo root:
%   demo_wing_panair_semispan

clear; close all; clc;

%% ------------------------------------------------------------------------
% User settings
%% ------------------------------------------------------------------------
rootAirfoilFile = 'airfoil_sharpTE.dat';
tipAirfoilFile  = 'airfoil_sharpTE.dat';
pFit = 3;
fitMethod = 'centripetal';

% Semispan wing geometry
semispan    = 5.0;   % exported half-span only; root lies on symmetry plane y=0
rootChord   = 2.0;
tipChord    = 1.0;
leSweepDeg  = 20.0;
incDeg      = 2.0;
tipTwistDeg = -3.0;
dihedralDeg = 5.0;
loftDegree  = 1;     % root + tip only

% Mesh density
nuWing = 41; % chordwise points
nvWing = 21; % spanwise points
nChordTip = 25; % chordwise points used on the tip cap (can be coarser than wing)
nuTip=nChordTip;
nvTip     = 5;  % upper->lower resolution on the tip cap
nuWake = 2;  % two streamwise wake lines
nvWake = nvWing; % same spanwise discretization as wing TE

spacingUWing = 'cosine';
spacingVWing = 'linear';
spacingUTip  = 'uniform'; % use milder end clustering on the tip cap
spacingVTip  = 'linear';

% Wake
wakeLength = 50 * rootChord; % PANAIR tutorial suggests about 25-50 cbar wake
wakeXDir   = [1, 0, 0];      % simple streamwise wake

% Analysis conditions for AUX file
mach  = 0.20;
alpha = 2.0;
cbar  = 0.5 * (rootChord + tipChord);
span  = 2.0 * semispan;                               % full reference span
sref  = 2.0 * semispan * 0.5 * (rootChord + tipChord); % full wing area
xref  = 0.25 * cbar;
zref  = 0.0;

% Output
outDir  = 'exports';
base    = 'wing_semispan';
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

rootLE = [0, 0, 0]; % root on XZ plane for symmetry
tipLE  = [semispan * tan(sweepRad), semispan, semispan * tan(dihedRad)];
fprintf(' root LE = [%.6f %.6f %.6f]\n', rootLE);
fprintf(' tip LE = [%.6f %.6f %.6f]\n', tipLE);

Troot = makeWingSectionTransformLocal(rootChord, incRad, rootLE);
Ttip  = makeWingSectionTransformLocal(tipChord, incRad + tipTwistRad, tipLE);

%% ------------------------------------------------------------------------
% 3) Transform split section curves into 3D
%% ------------------------------------------------------------------------
fprintf('\n--- 3. Transform section curves into 3D ---\n');
Cup_root = rootSection.upper.transform(Troot);
Clo_root = rootSection.lower.transform(Troot);
Cup_tip  = tipSection.upper.transform(Ttip);
Clo_tip  = tipSection.lower.transform(Ttip);

%% ------------------------------------------------------------------------
% 4) Loft wing surfaces
%% ------------------------------------------------------------------------
fprintf('\n--- 4. Build upper / lower wing surfaces ---\n');
Supper = geom.NURBSSurface.loft({Cup_root, Cup_tip}, loftDegree, 'chord');
Slower = geom.NURBSSurface.loft({Clo_root, Clo_tip}, loftDegree, 'chord');

%% ------------------------------------------------------------------------
% 5) Mesh wing surfaces and build a ruled-surface tip with reduced chordwise density
%% ------------------------------------------------------------------------
fprintf('5. Mesh wing surfaces and build ruled-surface tip cap ---');
meshUpper = Supper.isoMesh(nuWing, nvWing, ...
    'SpacingU', spacingUWing, 'SpacingV', spacingVWing);
meshLower = Slower.isoMesh(nuWing, nvWing, ...
    'SpacingU', spacingUWing, 'SpacingV', spacingVWing);

% Return to the original ruled-surface tip construction, but use a coarser
% chordwise tip mesh so the collapsed LE/TE ends do not create ultra-tiny
% first/last panel strips.
tipCurves = geom.NURBSSurface.makeCompatibleCurves({Cup_tip, Clo_tip});
Cup_tip_c = tipCurves{1};
Clo_tip_c = tipCurves{2};
Stip = geom.NURBSSurface.ruled(Cup_tip_c, Clo_tip_c);

meshTip = Stip.isoMesh(nuTip, nvTip, 'SpacingU', spacingUTip, 'SpacingV', spacingVTip);

auditRepeatedPointLinesLocal(meshTip, 'WingTip');
auditDegeneratePanelsLocal(meshTip, 'WingTip');

fprintf(' upper mesh quads = %d', size(meshUpper.connectivity,1));
fprintf(' lower mesh quads = %d', size(meshLower.connectivity,1));
fprintf(' tip mesh quads   = %d', size(meshTip.connectivity,1));

%% ------------------------------------------------------------------------
% 6) Build a simple wake sheet from the common trailing edge
%% ------------------------------------------------------------------------
fprintf('\n--- 6. Build semispan wake network ---\n');
TEu = [meshUpper.X(end,:).', meshUpper.Y(end,:).', meshUpper.Z(end,:).'];
TEl = [meshLower.X(end,:).', meshLower.Y(end,:).', meshLower.Z(end,:).'];
dTE = vecnorm(TEu - TEl, 2, 2);
fprintf(' TE coincidence: max = %.3e, RMS = %.3e\n', max(dTE), sqrt(mean(dTE.^2)));

TE      = 0.5 * (TEu + TEl); % robust common trailing edge
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

meshWake = transposeStructuredMeshLocal(meshWake);

fprintf(' wake mesh quads = %d\n', size(meshWake.connectivity,1));

%% ------------------------------------------------------------------------
% 7) Plot geometry and networks
%% ------------------------------------------------------------------------
figure('Name','1 - Semispan wing surfaces and wake');
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
% 8) Export LaWGS and AUX
%% ------------------------------------------------------------------------
fprintf('\n--- 8. Export LaWGS and AUX ---\n');
netUpper = geom.MeshWriter.meshToNetwork(meshUpper, 'WingUpper', 0, 0);
netLower = geom.MeshWriter.meshToNetwork(meshLower, 'WingLower', 0, 0);
netTip   = geom.MeshWriter.meshToNetwork(meshTip,   'WingTip',   0, 0);
netWake  = geom.MeshWriter.meshToNetwork(meshWake,  'WingWake',  0, 0);

nets = {netUpper, netLower, netTip, netWake};
geom.MeshWriter.toWGS(nets, wgsFile);
writePaninAuxLocal(auxFile, wgsFile, nets, mach, alpha, cbar, span, sref, xref, zref);

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
    % Row-major quad connectivity matching a structured nu x nv node grid.
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

function meshTip = buildTipMeshFromEdgesLocal(Pupper, Plower, nClose)
    % Build a structured tip closeout directly from corresponding upper/lower
    % edge samples. Pupper and Plower are [nChord x 3] arrays in the same
    % chordwise order. The resulting network is [nClose x nChord], so any
    % LE/TE collapse appears on network edges rather than repeated-point lines.

    if size(Pupper,1) ~= size(Plower,1)
        error('Upper/lower tip edges must have the same number of chordwise points.');
    end
    if size(Pupper,2) ~= 3 || size(Plower,2) ~= 3
        error('Tip edge arrays must be n-by-3.');
    end
    if nClose < 2
        error('Tip closeout needs at least two rows.');
    end

    nChord = size(Pupper,1);
    t = linspace(0, 1, nClose).';

    meshTip = struct();
    meshTip.X = zeros(nClose, nChord);
    meshTip.Y = zeros(nClose, nChord);
    meshTip.Z = zeros(nClose, nChord);
    for i = 1:nClose
        Pi = (1 - t(i)) * Pupper + t(i) * Plower;
        meshTip.X(i,:) = Pi(:,1).';
        meshTip.Y(i,:) = Pi(:,2).';
        meshTip.Z(i,:) = Pi(:,3).';
    end

    meshTip.nu = nClose;
    meshTip.nv = nChord;
    meshTip.connectivity = buildConnectivityLocal(meshTip.nu, meshTip.nv);
    meshTip.normals = [];

    % Keep the midpoint normal roughly +Y when possible.
    if size(meshTip.X,1) > 1 && size(meshTip.X,2) > 1
        i = max(1, floor(size(meshTip.X,1)/2));
        j = max(1, floor(size(meshTip.X,2)/2));
        p00 = [meshTip.X(i,j),   meshTip.Y(i,j),   meshTip.Z(i,j)];
        p10 = [meshTip.X(i+1,j), meshTip.Y(i+1,j), meshTip.Z(i+1,j)];
        p01 = [meshTip.X(i,j+1), meshTip.Y(i,j+1), meshTip.Z(i,j+1)];
        n = cross(p10 - p00, p01 - p00);
        if n(2) < 0
            meshTip.X = flipud(meshTip.X);
            meshTip.Y = flipud(meshTip.Y);
            meshTip.Z = flipud(meshTip.Z);
        end
    end
end

function auditDegeneratePanelsLocal(mesh, label)
    if nargin < 2
        label = 'mesh';
    end

    tol = 1e-12;
    bad = 0;
    amin = inf;
    for i = 1:(size(mesh.X,1)-1)
        for j = 1:(size(mesh.X,2)-1)
            p11 = [mesh.X(i,j),     mesh.Y(i,j),     mesh.Z(i,j)];
            p12 = [mesh.X(i,j+1),   mesh.Y(i,j+1),   mesh.Z(i,j+1)];
            p21 = [mesh.X(i+1,j),   mesh.Y(i+1,j),   mesh.Z(i+1,j)];
            p22 = [mesh.X(i+1,j+1), mesh.Y(i+1,j+1), mesh.Z(i+1,j+1)];
            a1 = 0.5 * norm(cross(p12-p11, p21-p11));
            a2 = 0.5 * norm(cross(p22-p12, p21-p12));
            a = a1 + a2;
            amin = min(amin, a);
            if a < tol
                bad = bad + 1;
            end
        end
    end

    fprintf(' %s min panel area = %.3e; degenerate panels (< %.1e) = %d\n', ...
        label, amin, tol, bad);
end

function Pout = resamplePolylineLocal(Pin, nOut, spacingMode)
    if nargin < 3 || isempty(spacingMode)
        spacingMode = 'uniform';
    end
    if size(Pin,1) < 2
        error('Need at least two input points to resample a polyline.');
    end
    if nOut < 2
        error('Need at least two output points to resample a polyline.');
    end

    d = vecnorm(diff(Pin,1,1), 2, 2);
    s = [0; cumsum(d)];
    if s(end) <= 0
        error('Polyline has zero total arclength.');
    end
    s = s / s(end);

    switch lower(strtrim(spacingMode))
        case {'uniform','linear'}
            t = linspace(0,1,nOut).';
        case 'cosine'
            th = linspace(0, pi, nOut).';
            t = 0.5 * (1 - cos(th));
        otherwise
            error('Unsupported spacing mode for tip resampling: %s', spacingMode);
    end

    Pout = zeros(nOut, 3);
    for k = 1:3
        Pout(:,k) = interp1(s, Pin(:,k), t, 'pchip');
    end
    Pout(1,:)   = Pin(1,:);
    Pout(end,:) = Pin(end,:);
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

function auditRepeatedPointLinesLocal(mesh, name)
    if nargin < 2
        name = 'mesh';
    end

    repeatedCount = 0;
    for i = 1:size(mesh.X,1)
        P = [mesh.X(i,:).', mesh.Y(i,:).', mesh.Z(i,:).'];
        if max(vecnorm(P - P(1,:), 2, 2)) < 1e-12
            repeatedCount = repeatedCount + 1;
            fprintf(' WARNING: %s line %d is a repeated-point line.\n', name, i);
        end
    end

    if repeatedCount == 0
        fprintf(' %s has no repeated-point lines.\n', name);
    else
        fprintf(' %s repeated-point lines: %d\n', name, repeatedCount);
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
