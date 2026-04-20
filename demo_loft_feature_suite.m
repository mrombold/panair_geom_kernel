%% demo_loft_feature_suite.m
% Simple checkout script for public features in geom.Loft.
%
% This is intentionally lightweight:
%  - small deterministic geometries
%  - numeric pass/fail checks
%  - one optional summary figure at the end
%
% Run from the repo root, or add the repo root to the MATLAB path first.

clear; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir)
    addpath(script_dir);
end

if ~exist('+geom', 'dir') && isempty(which('geom.Loft'))
    error(['Cannot find +geom package.\n' ...
        'Run this from the repo root or add the geometry kernel to path.']);
end

fprintf('geom package found at: %s\n', pwd);
fprintf('\n=== Loft Feature Suite ===\n');

TOL_TIGHT = 1.0e-10;
TOL_MED   = 1.0e-8;
checksPassed = 0;
checksTotal  = 0;

%% 1) sampleCurveAtStation
fprintf('\n--- 1. sampleCurveAtStation ---\n');
Cline = geom.NURBSCurve([0 0 0; 2 1 0], 1);
[u1, p1, info1] = geom.Loft.sampleCurveAtStation(Cline, 1.0, 'Axis', 'x');
expected1 = [1.0, 0.5, 0.0];
checksTotal = checksTotal + 3;
checksPassed = checksPassed + expectSmall(norm(p1 - expected1), TOL_TIGHT, ...
    'sampleCurveAtStation point');
checksPassed = checksPassed + expectSmall(abs(u1 - 0.5), TOL_TIGHT, ...
    'sampleCurveAtStation parameter');
checksPassed = checksPassed + expectTrue(isstruct(info1) && isfield(info1,'residual') && abs(info1.residual) < TOL_TIGHT, ...
    'sampleCurveAtStation residual info');
fprintf(' u = %.12f\n', u1);
fprintf(' point = [% .6f % .6f % .6f]\n', p1);

%% 2) intersectCurveWithPlane
fprintf('\n--- 2. intersectCurveWithPlane ---\n');
[u2, p2, info2] = geom.Loft.intersectCurveWithPlane(Cline, [1 0 0], 1.0);
checksTotal = checksTotal + 2;
checksPassed = checksPassed + expectSmall(norm(p2 - expected1), TOL_TIGHT, ...
    'intersectCurveWithPlane point');
checksPassed = checksPassed + expectTrue(abs(info2.residual) < TOL_TIGHT && abs(u2-u1) < TOL_TIGHT, ...
    'intersectCurveWithPlane consistency');

%% 3) makePlaneFrame / toPlane2D / fromPlane2D
fprintf('\n--- 3. Plane frame helpers ---\n');
origin = [1 2 3];
normal = [0 0 2];
xHint  = [1 0 0];
frame = geom.Loft.makePlaneFrame(origin, normal, xHint);
P3 = [1 2 3; 2 2 3; 1 3 3; 0.5 1.5 3];
uv = geom.Loft.toPlane2D(P3, frame);
P3_rt = geom.Loft.fromPlane2D(uv, frame);
checksTotal = checksTotal + 3;
checksPassed = checksPassed + expectSmall(norm(frame.R' * frame.R - eye(3), 'fro'), 1e-12, ...
    'makePlaneFrame orthonormality');
checksPassed = checksPassed + expectSmall(norm(P3_rt - P3, 'fro'), TOL_TIGHT, ...
    'toPlane2D/fromPlane2D round trip');
checksPassed = checksPassed + expectSmall(abs(det(frame.R) - 1.0), 1e-12, ...
    'makePlaneFrame right-handed rotation');

%% 4) lineIntersection2D
fprintf('\n--- 4. lineIntersection2D ---\n');
[pInt, t1, t2] = geom.Loft.lineIntersection2D([0 0], [1 1], [1 0], [-1 1]);
expectedInt = [0.5 0.5];
checksTotal = checksTotal + 3;
checksPassed = checksPassed + expectSmall(norm(pInt - expectedInt), TOL_TIGHT, ...
    'lineIntersection2D point');
checksPassed = checksPassed + expectSmall(t1 - 0.5, TOL_TIGHT, ...
    'lineIntersection2D t1');
checksPassed = checksPassed + expectSmall(t2 - 0.5, TOL_TIGHT, ...
    'lineIntersection2D t2');

%% 5) limingConic
fprintf('\n--- 5. limingConic ---\n');
P0 = [0.0, 0.0, 1.0];
P1 = [1.2, 0.0, 0.0];
T  = [0.0, 0.0, 0.0];
S  = [0.42, 0.0, 0.38];
[Clim, metaLim] = geom.Loft.limingConic(P0, P1, T, S);
P0chk = Clim.evaluate(0.0);
P1chk = Clim.evaluate(1.0);
Pschk = Clim.evaluate(metaLim.parameter);
checksTotal = checksTotal + 4;
checksPassed = checksPassed + expectSmall(norm(P0chk - P0), TOL_TIGHT, ...
    'limingConic start point');
checksPassed = checksPassed + expectSmall(norm(P1chk - P1), TOL_TIGHT, ...
    'limingConic end point');
checksPassed = checksPassed + expectSmall(norm(Pschk - S), max(TOL_MED, 10*metaLim.fitToleranceUsed), ...
    'limingConic shoulder point');
checksPassed = checksPassed + expectTrue(isfield(metaLim,'weight') && metaLim.weight > 0, ...
    'limingConic positive weight');
fprintf(' solved u = %.8f, w = %.8f, shoulder err = %.3e\n', ...
    metaLim.parameter, metaLim.weight, norm(Pschk - S));

%% 6) buildFuselageSectionAtStation
fprintf('\n--- 6. buildFuselageSectionAtStation ---\n');
station = 1.5;
secSpec = makeSimpleSectionSpec();
section = geom.Loft.buildFuselageSectionAtStation(station, ...
    'UpperProfile',   secSpec.UpperProfile, ...
    'LowerProfile',   secSpec.LowerProfile, ...
    'MaxBreadth',     secSpec.MaxBreadth, ...
    'UpperShoulder',  secSpec.UpperShoulder, ...
    'LowerShoulder',  secSpec.LowerShoulder, ...
    'Axis',           'x');

checksTotal = checksTotal + 5;
checksPassed = checksPassed + expectSmall(norm(section.points.top        - [station 0 secSpec.zTop]), TOL_TIGHT, ...
    'buildFuselageSectionAtStation top sample');
checksPassed = checksPassed + expectSmall(norm(section.points.bottom     - [station 0 secSpec.zBot]), TOL_TIGHT, ...
    'buildFuselageSectionAtStation bottom sample');
checksPassed = checksPassed + expectSmall(norm(section.points.maxBreadth - [station secSpec.yMax 0]), TOL_TIGHT, ...
    'buildFuselageSectionAtStation max breadth sample');
checksPassed = checksPassed + expectSmall(norm(section.upperCurve.evaluate(section.upperMeta.parameter) - section.points.upperShoulder), ...
    max(TOL_MED, 10*section.upperMeta.fitToleranceUsed), ...
    'buildFuselageSectionAtStation upper shoulder fit');
checksPassed = checksPassed + expectSmall(norm(section.lowerCurve.evaluate(section.lowerMeta.parameter) - section.points.lowerShoulder), ...
    max(TOL_MED, 10*section.lowerMeta.fitToleranceUsed), ...
    'buildFuselageSectionAtStation lower shoulder fit');

%% 7) buildSectionFamily
fprintf('\n--- 7. buildSectionFamily ---\n');
stations = [0.0 1.5 3.0];
sections = geom.Loft.buildSectionFamily(stations, ...
    'UpperProfile',   secSpec.UpperProfile, ...
    'LowerProfile',   secSpec.LowerProfile, ...
    'MaxBreadth',     secSpec.MaxBreadth, ...
    'UpperShoulder',  secSpec.UpperShoulder, ...
    'LowerShoulder',  secSpec.LowerShoulder, ...
    'Axis',           'x');

checksTotal = checksTotal + 3;
checksPassed = checksPassed + expectTrue(iscell(sections) && numel(sections) == numel(stations), ...
    'buildSectionFamily output shape');
checksPassed = checksPassed + expectSmall(sections{1}.stationValue - stations(1), TOL_TIGHT, ...
    'buildSectionFamily first station');
checksPassed = checksPassed + expectSmall(sections{end}.stationValue - stations(end), TOL_TIGHT, ...
    'buildSectionFamily last station');

%% 8) loftSections
fprintf('\n--- 8. loftSections ---\n');
upperCurves = cellfun(@(s) s.upperCurve, sections, 'UniformOutput', false);
lowerCurves = cellfun(@(s) s.lowerCurve, sections, 'UniformOutput', false);
Supper = geom.Loft.loftSections(upperCurves, 2, 'chord');
Slower = geom.Loft.loftSections(lowerCurves, 2, 'chord');

pcU = Supper.evaluate(0.35, 0.5);
pcL = Slower.evaluate(0.35, 0.5);
checksTotal = checksTotal + 4;
checksPassed = checksPassed + expectTrue(isa(Supper, 'geom.NURBSSurface'), ...
    'loftSections upper surface type');
checksPassed = checksPassed + expectTrue(isa(Slower, 'geom.NURBSSurface'), ...
    'loftSections lower surface type');
checksPassed = checksPassed + expectTrue(all(isfinite(pcU)), ...
    'loftSections upper evaluate');
checksPassed = checksPassed + expectTrue(all(isfinite(pcL)), ...
    'loftSections lower evaluate');
fprintf(' upper surface sample = [% .6f % .6f % .6f]\n', pcU);
fprintf(' lower surface sample = [% .6f % .6f % .6f]\n', pcL);

%% Summary
fprintf('\n=== Summary ===\n');
fprintf('Passed %d / %d checks.\n', checksPassed, checksTotal);
if checksPassed ~= checksTotal
    error('demo_loft_feature_suite:FailedChecks', ...
        'One or more Loft feature checks failed.');
end
fprintf('All Loft feature checks passed.\n');

%% Optional visualization
figure('Name', 'Loft feature suite overview');
subplot(1,2,1); hold on; grid on; axis equal; view(3);
plotCurve3D(Cline, 'k-', 1.5);
plotCurve3D(Clim, 'r-', 2.0);
plot3(P0(1), P0(2), P0(3), 'go', 'MarkerFaceColor', 'g');
plot3(P1(1), P1(2), P1(3), 'bo', 'MarkerFaceColor', 'b');
plot3(T(1), T(2), T(3), 'ko', 'MarkerFaceColor', 'y');
plot3(S(1), S(2), S(3), 'mo', 'MarkerFaceColor', 'm');
title('Curve-side Loft checks');
xlabel('x'); ylabel('y'); zlabel('z');
legend({'station test curve','Liming conic','P0','P1','T','S'}, 'Location','best');

subplot(1,2,2); hold on; grid on; axis equal; view(3);
plotCurve3D(section.upperCurve, 'b-', 1.5);
plotCurve3D(section.lowerCurve, 'r-', 1.5);
plot3(section.points.top(1), section.points.top(2), section.points.top(3), 'ko', 'MarkerFaceColor', 'g');
plot3(section.points.bottom(1), section.points.bottom(2), section.points.bottom(3), 'ko', 'MarkerFaceColor', 'g');
plot3(section.points.maxBreadth(1), section.points.maxBreadth(2), section.points.maxBreadth(3), 'ko', 'MarkerFaceColor', 'c');
plot3(section.points.upperShoulder(1), section.points.upperShoulder(2), section.points.upperShoulder(3), 'ko', 'MarkerFaceColor', 'm');
plot3(section.points.lowerShoulder(1), section.points.lowerShoulder(2), section.points.lowerShoulder(3), 'ko', 'MarkerFaceColor', 'm');
title('Section construction checks');
xlabel('x'); ylabel('y'); zlabel('z');
legend({'upper section','lower section','top/bottom','max breadth','shoulders'}, 'Location','best');


function spec = makeSimpleSectionSpec()
    % A deliberately simple family of guide curves.
    %
    % Each guide curve is monotonic in x, so station sampling is easy and
    % deterministic. Shoulder guide curves are built from exact rational
    % quadratic section curves so buildFuselageSectionAtStation gets a clean,
    % self-consistent test case.

    spec.zTop = 1.0;
    spec.zBot = -1.0;
    spec.yMax = 1.2;
    xEnds = [0.0; 3.0];

    Ptop = [1.5, 0.0, spec.zTop];
    Pbot = [1.5, 0.0, spec.zBot];
    Pmax = [1.5, spec.yMax, 0.0];
    Tint = [1.5, 0.0, 0.0];

    Uq = [0 0 0 1 1 1];
    Wq = [1; 0.70; 1];
    Cup = geom.NURBSCurve([Ptop; Tint; Pmax], 2, Uq, Wq);
    Clo = geom.NURBSCurve([Pmax; Tint; Pbot], 2, Uq, Wq);

    us = 0.45;
    Sup = Cup.evaluate(us);
    Slo = Clo.evaluate(us);

    spec.UpperProfile  = geom.NURBSCurve([xEnds, [0.0; 0.0], [spec.zTop; spec.zTop]], 1);
    spec.LowerProfile  = geom.NURBSCurve([xEnds, [0.0; 0.0], [spec.zBot; spec.zBot]], 1);
    spec.MaxBreadth    = geom.NURBSCurve([xEnds, [spec.yMax; spec.yMax], [0.0; 0.0]], 1);
    spec.UpperShoulder = geom.NURBSCurve([xEnds, [Sup(2); Sup(2)], [Sup(3); Sup(3)]], 1);
    spec.LowerShoulder = geom.NURBSCurve([xEnds, [Slo(2); Slo(2)], [Slo(3); Slo(3)]], 1);
end

function ok = expectSmall(val, tol, label)
    ok = isfinite(val) && abs(val) <= tol;
    if ok
        fprintf(' PASS  %-42s  err = %.3e\n', label, abs(val));
    else
        fprintf(' FAIL  %-42s  err = %.3e  tol = %.3e\n', label, abs(val), tol);
    end
end

function ok = expectTrue(cond, label)
    ok = logical(cond);
    if ok
        fprintf(' PASS  %s\n', label);
    else
        fprintf(' FAIL  %s\n', label);
    end
end

function plotCurve3D(C, style, lw)
    if nargin < 2 || isempty(style), style = 'b-'; end
    if nargin < 3 || isempty(lw), lw = 1.5; end
    u = linspace(C.domain(1), C.domain(2), 200);
    P = C.evaluate(u);
    plot3(P(:,1), P(:,2), P(:,3), style, 'LineWidth', lw);
end
