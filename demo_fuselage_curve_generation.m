%% demo_fuselage_curve_generation.m
% Curve-only half-fuselage demo for the panair_geom_kernel.
%
% Goal:
%   - generate top / bottom longitudinal profile curves on the y=0 plane
%   - generate half cross-section curves at several x-stations
%   - plot raw construction points, fitted NURBS curves, and optional control polygons
%   - make it easy to tweak fuselage shaping before lofting surfaces
%
% Run from repo root:
%   demo_fuselage_curve_generation

clear; close all; clc;

%% ------------------------------------------------------------------------
% User settings
%% ------------------------------------------------------------------------
L = 12.0;                 % fuselage length
xStations = [0.0 0.4 1.0 2.0 4.0 6.5 8.5 10.0 11.2 12.0];

% Longitudinal shaping (all lengths in model units)
%   zTop(x) = zc(x) + h(x)
%   zBot(x) = zc(x) - h(x)
maxHalfWidth  = 1.10;
maxHalfHeight = 1.10;

% Vertical centerline shift along body (camber-ish body center motion)
centerShiftAmp = 0.2;

% Section squareness / fullness
% n = 2   -> ellipse
% n > 2   -> squarer / fuller sided superellipse
% n < 2   -> pointier / diamond-ish
superNose  = 2.0;
superMid   = 2.4;
superTail  = 2.1;

% Curve fitting / plotting
pProfile = 3;
pSection = 3;
fitMethod = 'centripetal';
sectionSamples = 81;      % points used to define each half-section
profilePlotPts = 300;
sectionPlotPts = 160;
showControlPolygons = true;
showConstructionPts = true;
sectionScaleXYZ = [1 1 1]; % optional anisotropic display scaling of section points

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
fprintf('geom package found at: %s\n', pwd);

%% ------------------------------------------------------------------------
% 1) Build longitudinal profile construction points
%% ------------------------------------------------------------------------
fprintf('\n--- 1. Build longitudinal profile points ---\n');

nStations = numel(xStations);
Ptop = zeros(nStations, 3);
Pbot = zeros(nStations, 3);
Pctr = zeros(nStations, 3);

for k = 1:nStations
    x = xStations(k);
    xi = x / L;

    w = halfWidthEnvelopeLocal(xi, maxHalfWidth);
    h = halfHeightEnvelopeLocal(xi, maxHalfHeight);
    zc = centerShiftLocal(xi, centerShiftAmp);

    Ptop(k,:) = [x, 0, zc + h];
    Pbot(k,:) = [x, 0, zc - h];
    Pctr(k,:) = [x, 0, zc];

    fprintf(' station %2d: x = %7.3f, half-width = %.4f, half-height = %.4f, zc = %.4f\n', ...
        k, x, w, h, zc);
end

%% ------------------------------------------------------------------------
% 2) Fit top / bottom / center longitudinal curves
%% ------------------------------------------------------------------------
fprintf('\n--- 2. Fit longitudinal NURBS curves ---\n');
Ctop = geom.NURBSCurve.globalInterp(Ptop, pProfile, fitMethod);
Cbot = geom.NURBSCurve.globalInterp(Pbot, pProfile, fitMethod);
Cctr = geom.NURBSCurve.globalInterp(Pctr, pProfile, fitMethod);

fprintf(' built Ctop, Cbot, Cctr\n');

%% ------------------------------------------------------------------------
% 3) Build half cross-section curves
%% ------------------------------------------------------------------------
fprintf('\n--- 3. Build half cross-section curves ---\n');

Csec = cell(nStations,1);
secPts = cell(nStations,1);
secMeta = struct('x', cell(nStations,1), 'w', [], 'h', [], 'zc', [], 'n', []);

for k = 1:nStations
    x = xStations(k);
    xi = x / L;

    w = halfWidthEnvelopeLocal(xi, maxHalfWidth);
    h = halfHeightEnvelopeLocal(xi, maxHalfHeight);
    zc = centerShiftLocal(xi, centerShiftAmp);
    nsec = sectionExponentLocal(xi, superNose, superMid, superTail);

    Pk = makeHalfSectionPointsLocal(x, w, h, zc, nsec, sectionSamples);
    Ck = geom.NURBSCurve.globalInterp(Pk, pSection, fitMethod);

    secPts{k} = Pk;
    Csec{k} = Ck;
    secMeta(k).x = x;
    secMeta(k).w = w;
    secMeta(k).h = h;
    secMeta(k).zc = zc;
    secMeta(k).n = nsec;

    fprintf(' section %2d built: x = %7.3f, n = %.3f\n', k, x, nsec);
end

%% ------------------------------------------------------------------------
% 4) Plot curves
%% ------------------------------------------------------------------------
fprintf('\n--- 4. Plot curves ---\n');
figure('Name','Half-fuselage curve generation');
hold on; grid on; axis equal; view(3);

% Plot longitudinal curves
uTop = linspace(Ctop.domain(1), Ctop.domain(2), profilePlotPts);
uBot = linspace(Cbot.domain(1), Cbot.domain(2), profilePlotPts);
uCtr = linspace(Cctr.domain(1), Cctr.domain(2), profilePlotPts);
Qtop = Ctop.evaluate(uTop);
Qbot = Cbot.evaluate(uBot);
Qctr = Cctr.evaluate(uCtr);

plot3(Qtop(:,1), Qtop(:,2), Qtop(:,3), 'r-', 'LineWidth', 2.0);
plot3(Qbot(:,1), Qbot(:,2), Qbot(:,3), 'b-', 'LineWidth', 2.0);
plot3(Qctr(:,1), Qctr(:,2), Qctr(:,3), 'k--', 'LineWidth', 1.0);

if showConstructionPts
    plot3(Ptop(:,1), Ptop(:,2), Ptop(:,3), 'ro', 'MarkerSize', 5, 'LineWidth', 1.0);
    plot3(Pbot(:,1), Pbot(:,2), Pbot(:,3), 'bo', 'MarkerSize', 5, 'LineWidth', 1.0);
end

if showControlPolygons
    plot3(Ctop.P(:,1), Ctop.P(:,2), Ctop.P(:,3), 'r:', 'LineWidth', 1.0);
    plot3(Cbot.P(:,1), Cbot.P(:,2), Cbot.P(:,3), 'b:', 'LineWidth', 1.0);
end

% Plot cross-sections
for k = 1:nStations
    Ck = Csec{k};
    uk = linspace(Ck.domain(1), Ck.domain(2), sectionPlotPts);
    Qk = Ck.evaluate(uk);
    Qk = Qk .* sectionScaleXYZ;

    plot3(Qk(:,1), Qk(:,2), Qk(:,3), '-', 'LineWidth', 1.4);

    if showConstructionPts
        Pk = secPts{k} .* sectionScaleXYZ;
        plot3(Pk(:,1), Pk(:,2), Pk(:,3), '.', 'MarkerSize', 8);
    end

    if showControlPolygons
        Pctrl = Ck.P .* sectionScaleXYZ;
        plot3(Pctrl(:,1), Pctrl(:,2), Pctrl(:,3), ':', 'LineWidth', 0.8);
    end
end

% Symmetry plane reference line
plot3([0 L], [0 0], [0 0], 'k-.', 'LineWidth', 1.0);

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Half-fuselage defining curves');
legend({'Top profile','Bottom profile','Centerline'}, 'Location','best');

%% ------------------------------------------------------------------------
% 5) Plot section-only view for shape tuning
%% ------------------------------------------------------------------------
figure('Name','Half-section shapes in local yz view');
hold on; grid on; axis equal;

for k = 1:nStations
    Ck = Csec{k};
    uk = linspace(Ck.domain(1), Ck.domain(2), sectionPlotPts);
    Qk = Ck.evaluate(uk);
    plot(Qk(:,2), Qk(:,3), 'LineWidth', 1.5);

    if showConstructionPts
        Pk = secPts{k};
        plot(Pk(:,2), Pk(:,3), '.');
    end
end

xlabel('Y'); ylabel('Z');
title('Half cross-sections in local section view');
plot([0 0], ylim, 'k--');

%% ------------------------------------------------------------------------
% 6) Optional tangent check at symmetry-plane endpoints
%% ------------------------------------------------------------------------
fprintf('\n--- 5. Endpoint tangent check at y = 0 plane ---\n');
for k = 1:nStations
    Ck = Csec{k};
    dA = Ck.derivative(Ck.domain(1), 1);
    dB = Ck.derivative(Ck.domain(2), 1);
    fprintf([' section %2d: |dy/ds| start = %.3e, end = %.3e  ' ...
             '(small means near-tangent to symmetry plane)\n'], ...
             k, abs(dA(2)), abs(dB(2)));
end

fprintf('\nDone.\n');
fprintf('Suggested next experiments:\n');
fprintf('  - change halfWidthEnvelopeLocal / halfHeightEnvelopeLocal\n');
fprintf('  - change centerShiftLocal\n');
fprintf('  - change superNose / superMid / superTail\n');
fprintf('  - reduce sectionSamples if you want lighter construction curves\n');
