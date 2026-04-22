%% demo_fuselage_surfaces_semispan.m
% Surface-based half-fuselage demo (symmetry plane at Y = 0).
% Similar in style to demo_wing_surfaces_panair_semispan, but for a body.
%
% Builds:
%   - upper half-fuselage surface from longitudinal profile + sections
%   - lower half-fuselage surface from longitudinal profile + sections
%   - optional nose closeout surface on the symmetry plane
%   - optional tail closeout surface on the symmetry plane
%
% Geometry only for now; no panel meshing/export yet.
%
% Run from repo root:
%   demo_fuselage_surfaces_semispan

clear; close all; clc;

%% -----------------------------------------------------------------------
% User settings
%% -----------------------------------------------------------------------
% Overall dimensions
L_fuse = 12.0;     % fuselage length
R_max  = 1.10;     % maximum body radius

% Longitudinal station layout
xNoseEnd = 2.0;
xCylEnd  = 8.5;

% NURBS / loft settings
pProfile = 3;
pSection = 3;
loftDegree = 3;

% Section set (nose -> tail)
xStations = [ ...
    0.0, ...
    0.4, ...
    1.0, ...
    2.0, ...
    4.0, ...
    6.5, ...
    8.5, ...
    10.0, ...
    11.2, ...
    12.0];

% Vertical shape controls at each section
%   zcFrac = vertical center shift / R_max
%   hFrac  = total section height / (2*R_max)
%   wFrac  = total section width  / (2*R_max)
zcFrac = [0.00, 0.00, 0.01, 0.02, 0.03, 0.03, 0.02, 0.01, 0.00, 0.00];
hFrac  = [0.00, 0.30, 0.68, 0.94, 1.00, 1.00, 0.98, 0.82, 0.42, 0.00];
wFrac  = [0.00, 0.22, 0.58, 0.90, 1.00, 1.00, 0.98, 0.76, 0.34, 0.00];

% Cross-section shape exponent:
%   n = 2   -> ellipse
%   n > 2   -> boxier/superellipse-like
%   n < 2   -> pointier
nShape  = [2.0, 2.0, 2.2, 2.4, 2.6, 2.6, 2.5, 2.3, 2.1, 2.0];

% Sampling for plotting only
nProfilePlot = 300;
nSectionPlot = 160;

%% -----------------------------------------------------------------------
% Path / package check
%% -----------------------------------------------------------------------
script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir)
    addpath(script_dir);
end
if ~exist('+geom', 'dir') && isempty(which('geom.NURBSCurve'))
    error(['Cannot find +geom package.\n' ...
           'Run this from the repo root or add the geometry kernel to path.']);
end
fprintf('geom package found at: %s\n', script_dir);

%% -----------------------------------------------------------------------
% 1) Build target upper / lower longitudinal profile points
%% -----------------------------------------------------------------------
fprintf('\n--- 1. Build longitudinal profile curves ---\n');

xProf = [0.0, 0.2, 0.8, xNoseEnd, 4.0, xCylEnd, 10.0, 11.2, L_fuse];
zTop  = [0.0, 0.18, 0.55, 0.95, 1.00, 0.99, 0.72, 0.28, 0.0] * R_max;
zBot  = [0.0,-0.16,-0.48,-0.88,-0.97,-0.96,-0.70,-0.24, 0.0] * R_max;

Ptop = [xProf(:), zeros(numel(xProf),1), zTop(:)];
Pbot = [xProf(:), zeros(numel(xProf),1), zBot(:)];

Ctop = geom.NURBSCurve.globalInterp(Ptop, pProfile, 'chord');
Cbot = geom.NURBSCurve.globalInterp(Pbot, pProfile, 'chord');

fprintf(' built upper and lower symmetry-plane profile curves\n');

%% -----------------------------------------------------------------------
% 2) Build body cross-sections from section controls
%% -----------------------------------------------------------------------
fprintf('\n--- 2. Build fuselage cross-sections ---\n');

nSec = numel(xStations);
Cup = cell(1, nSec);
Clo = cell(1, nSec);
sections = struct('x', cell(1,nSec), 'zc', cell(1,nSec), 'a', cell(1,nSec), 'b', cell(1,nSec));

for k = 1:nSec
    xk = xStations(k);
    zck = zcFrac(k) * R_max;
    bk  = hFrac(k)  * R_max;   % semi-height proxy
    ak  = wFrac(k)  * R_max;   % semi-width  proxy
    nk  = nShape(k);

    [Cup{k}, Clo{k}] = buildHalfSectionCurvesLocal(xk, ak, bk, zck, nk, pSection);

    sections(k).x  = xk;
    sections(k).zc = zck;
    sections(k).a  = ak;
    sections(k).b  = bk;

    fprintf(' section %2d: x = %7.3f, half-width = %.4f, half-height = %.4f\n', ...
        k, xk, ak, bk);
end

%% -----------------------------------------------------------------------
% 3) Force exact profile/section consistency on symmetry line
%% -----------------------------------------------------------------------
fprintf('\n--- 3. Match section apex points to profile curves ---\n');

for k = 1:nSec
    xk = xStations(k);

    % Find corresponding profile points by nearest x in parameter sweep.
    pTop = evaluateCurveAtTargetXLocal(Ctop, xk, 2000);
    pBot = evaluateCurveAtTargetXLocal(Cbot, xk, 2000);

    Cup{k} = snapCurveEndsLocal(Cup{k}, pTop, 'first_last');
    Clo{k} = snapCurveEndsLocal(Clo{k}, pBot, 'first_last');
end

fprintf(' snapped upper/lower section end points to longitudinal profiles\n');

%% -----------------------------------------------------------------------
% 4) Loft upper and lower half-body surfaces
%% -----------------------------------------------------------------------
fprintf('\n--- 4. Loft upper / lower fuselage surfaces ---\n');

Supper = geom.NURBSSurface.loft(Cup, loftDegree, 'chord');
Slower = geom.NURBSSurface.loft(Clo, loftDegree, 'chord');

[Supper, upperSurfLabel, nSU] = orientSurfaceNormalLocal(Supper, 'up',   'FuselageUpper');
[Slower, lowerSurfLabel, nSL] = orientSurfaceNormalLocal(Slower, 'down', 'FuselageLower');

fprintf(' surface normals: FuselageUpper=%s [%.3f %.3f %.3f], FuselageLower=%s [%.3f %.3f %.3f]\n', ...
    upperSurfLabel, nSU(1), nSU(2), nSU(3), ...
    lowerSurfLabel, nSL(1), nSL(2), nSL(3));

%% -----------------------------------------------------------------------
% 5) Optional nose and tail symmetry closeout surfaces
%% -----------------------------------------------------------------------
fprintf('\n--- 5. Build nose / tail closeout surfaces on symmetry plane ---\n');

% These are ruled surfaces connecting upper and lower section curves at the
% degenerate nose and tail stations. They help visualize closure but are
% still just geometry objects for now.
curvesNose = geom.NURBSSurface.makeCompatibleCurves({Cup{1}, Clo{1}});
SNose = geom.NURBSSurface.ruled(curvesNose{1}, curvesNose{2});

curvesTail = geom.NURBSSurface.makeCompatibleCurves({Cup{end}, Clo{end}});
STail = geom.NURBSSurface.ruled(curvesTail{1}, curvesTail{2});

fprintf(' built nose and tail closeout surfaces\n');

%% -----------------------------------------------------------------------
% 6) Plot profiles, sections, and surfaces
%% -----------------------------------------------------------------------
fprintf('\n--- 6. Plot fuselage geometry ---\n');

figure('Name','Half-fuselage surfaces');
hold on; grid on; axis equal; view(3);

plotCurveLocal(Ctop, nProfilePlot, 'k-', 1.8);
plotCurveLocal(Cbot, nProfilePlot, 'k-', 1.8);

for k = 1:nSec
    plotCurveLocal(Cup{k}, nSectionPlot, '-', 0.9);
    plotCurveLocal(Clo{k}, nSectionPlot, '-', 0.9);
end

plotSurfaceLocal(Supper, 60, 30, [0.88 0.35 0.35], 0.72);
plotSurfaceLocal(Slower, 60, 30, [0.35 0.45 0.92], 0.72);
plotSurfaceLocal(SNose,  16, 10, [0.92 0.85 0.35], 0.78);
plotSurfaceLocal(STail,  16, 10, [0.92 0.85 0.35], 0.78);

plot3([0 L_fuse], [0 0], [0 0], 'k--', 'LineWidth', 1.0);

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Half-fuselage NURBS surfaces (symmetry plane at Y = 0)');
legend({'Upper profile','Lower profile','Sections','Upper surface','Lower surface', ...
        'Nose closeout','Tail closeout','Symmetry line'}, 'Location','best');

%% -----------------------------------------------------------------------
% 7) Print suggested next steps
%% -----------------------------------------------------------------------
fprintf('\nDone.\n');
fprintf('Suggested next step:\n');
fprintf('  - mesh Supper and Slower with isoMesh(...)\n');
fprintf('  - optionally mesh the nose/tail closeouts\n');
fprintf('  - export WGS networks for a half-body PANAIR model\n');

%% =======================================================================
% Local functions
%% =======================================================================
function [Cup, Clo] = buildHalfSectionCurvesLocal(x0, a, b, zc, nShape, p)
    % Builds upper and lower half-section curves in the Y-Z plane at x=x0.
    % Curves run from symmetry-plane crown/keel to outer side point and back
    % to symmetry plane so they are suitable for separate upper/lower lofts.

    if a < 1e-12 || b < 1e-12
        P = [x0 0 zc; x0 0 zc; x0 0 zc];
        Cup = geom.NURBSCurve.globalInterp(P, min(p,2), 'chord');
        Clo = geom.NURBSCurve.globalInterp(P, min(p,2), 'chord');
        return;
    end

    nTheta = 61;

    % Upper half: crown (Y=0,+b) -> side (+a,0) -> symmetry crown again is
    % not what we want. We instead define a single open curve from crown to
    % side so loft sections remain consistent by using all sections with the
    % same endpoint ordering. Same for lower.
    thetaUp = linspace(pi/2, 0, nTheta);
    thetaLo = linspace(-pi/2, 0, nTheta);

    [yu, zu] = superellipseQuarterLocal(a, b, zc, nShape, thetaUp);
    [yl, zl] = superellipseQuarterLocal(a, b, zc, nShape, thetaLo);

    Pup = [x0 * ones(numel(yu),1), yu(:), zu(:)];
    Plo = [x0 * ones(numel(yl),1), yl(:), zl(:)];

    Cup = geom.NURBSCurve.globalInterp(Pup, p, 'centripetal');
    Clo = geom.NURBSCurve.globalInterp(Plo, p, 'centripetal');
end

function [y, z] = superellipseQuarterLocal(a, b, zc, n, theta)
    c = cos(theta);
    s = sin(theta);

    y = a * sign(c) .* abs(c).^(2./n);
    z = zc + b * sign(s) .* abs(s).^(2./n);
end

function P = evaluateCurveAtTargetXLocal(C, xTarget, nSearch)
    us = linspace(C.domain(1), C.domain(2), nSearch);
    Pts = C.evaluate(us);
    [~, idx] = min(abs(Pts(:,1) - xTarget));
    P = Pts(idx,:);
end

function C2 = snapCurveEndsLocal(C, Ptarget, mode)
    P = C.P;
    W = C.W;

    switch lower(strtrim(mode))
        case 'first_last'
            % For quarter-section curves, the first control point sits on the
            % symmetry plane (top or bottom), so snap only that end.
            P(1,:) = Ptarget(:).';
        otherwise
            error('Unknown snap mode: %s', mode);
    end

    C2 = geom.NURBSCurve(P, C.p, C.U, W);
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
        otherwise
            error('Unknown expected normal direction: %s', expectedNormal);
    end
end

function plotCurveLocal(C, n, style, lw)
    u = linspace(C.domain(1), C.domain(2), n);
    P = C.evaluate(u);
    plot3(P(:,1), P(:,2), P(:,3), style, 'LineWidth', lw);
end

function plotSurfaceLocal(S, nu, nv, faceColor, faceAlpha)
    u = linspace(S.domainU(1), S.domainU(2), nu);
    v = linspace(S.domainV(1), S.domainV(2), nv);
    [U, V] = meshgrid(u, v);

    X = zeros(size(U));
    Y = zeros(size(U));
    Z = zeros(size(U));

    for jj = 1:size(U,1)
        for ii = 1:size(U,2)
            P = S.evaluate(U(jj,ii), V(jj,ii));
            X(jj,ii) = P(1);
            Y(jj,ii) = P(2);
            Z(jj,ii) = P(3);
        end
    end

    surf(X, Y, Z, ...
        'FaceColor', faceColor, ...
        'FaceAlpha', faceAlpha, ...
        'EdgeColor', 'k', ...
        'EdgeAlpha', 0.18);
end
