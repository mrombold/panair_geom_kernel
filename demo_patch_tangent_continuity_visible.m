%% demo_patch_tangent_continuity_visible
% A more visible continuity demo:
%   - plots transverse section cuts before/after
%   - plots edge tangent-ribbon vectors before/after
%   - keeps the geometry simple enough to see the G1 change

clear; close all; clc;

if ~exist('+geom', 'dir') && isempty(which('geom.NURBSSurface'))
    error('Run this from the repo root or add the geometry kernel to the MATLAB path.');
end
if isempty(which('geom.Patch'))
    error('geom.Patch not found. Add Patch.m to +geom first.');
end

fprintf('=== Visible Patch Tangent Continuity Demo ===\n');

xStations = [0 4 8 12];
yRowsA = [0.0 0.45 0.90 1.30];
yRowsB = [1.30 1.80 2.25 2.70];
sharedY = yRowsA(end);

QA = zeros(numel(yRowsA), numel(xStations), 3);
QB = zeros(numel(yRowsB), numel(xStations), 3);

for j = 1:numel(yRowsA)
    y = yRowsA(j);
    for i = 1:numel(xStations)
        x = xStations(i);
        z = 1.00 - 0.18*(y-0.25)^2 + 0.08*cos(pi*x/xStations(end));
        QA(j,i,:) = [x, y, z];
    end
end

for j = 1:numel(yRowsB)
    y = yRowsB(j);
    for i = 1:numel(xStations)
        x = xStations(i);

        if abs(y - sharedY) < 1e-12
            z = 1.00 - 0.18*(y-0.25)^2 + 0.08*cos(pi*x/xStations(end));
        else
            dy = y - sharedY;
            % Deliberately make the upper patch leave the shared edge with a
            % visibly different slope and curvature.
            z = 0.92 - 0.30*dy - 0.12*dy.^2 + 0.05*cos(pi*x/xStations(end) + 0.35);
        end

        QB(j,i,:) = [x, y, z];
    end
end

SA = geom.NURBSSurface.globalInterpNet(QA, 3, 3, 'centripetal', 'centripetal');
SB = geom.NURBSSurface.globalInterpNet(QB, 3, 3, 'centripetal', 'centripetal');

patchA = geom.Patch('PatchA', SA);
patchB = geom.Patch('PatchB', SB);

repBefore = patchB.g1ReportWithPatch('u0', patchA, 'u1');
repAfter  = patchB.enforceG1WithPatch('u0', patchA, 'u1', ...
    'Lambda', 1.0, 'AverageEdge', true);

fprintf('Before G1: max position err          = %.3e\n', repBefore.maxPositionError);
fprintf('Before G1: max cross-plane mismatch  = %.3e\n', repBefore.maxCrossPlaneMismatch);
fprintf('After  G1: max position err          = %.3e\n', repAfter.maxPositionError);
fprintf('After  G1: max cross-plane mismatch  = %.3e\n', repAfter.maxCrossPlaneMismatch);

figure('Name','Visible G1 continuity demo','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

%% 1) 3D before
nexttile; hold on; axis equal; grid on; view(3);
title('Before G1');
plotSurfaceWire(SA, 31, 21, [0.1 0.55 0.1]);
plotSurfaceWire(SB, 31, 21, [0.75 0.25 0.25]);
plotEdge(SA, 'u1', 'k-', 2.0);
plotEdge(SB, 'u0', 'r--', 1.5);
plotTangentRibbon(SA, 'u1', 9, 0.8, [0 0 0]);
plotTangentRibbon(SB, 'u0', 9, 0.8, [0.85 0 0]);
xlabel('x'); ylabel('y'); zlabel('z');

%% 2) 3D after
nexttile; hold on; axis equal; grid on; view(3);
title('After G1');
plotSurfaceWire(SA, 31, 21, [0.1 0.55 0.1]);
plotSurfaceWire(patchB.S, 31, 21, [0.0 0.45 0.85]);
plotEdge(SA, 'u1', 'k-', 2.0);
plotEdge(patchB.S, 'u0', 'b-', 1.5);
plotTangentRibbon(SA, 'u1', 9, 0.8, [0 0 0]);
plotTangentRibbon(patchB.S, 'u0', 9, 0.8, [0 0 0.85]);
xlabel('x'); ylabel('y'); zlabel('z');

%% 3) Section cuts before/after at mid-span x
nexttile; hold on; axis equal; grid on;
title('Section cut near mid-span');
xCut = 6.0;
PbA  = sampleCurveAtX(SA, xCut);
PbB0 = sampleCurveAtX(SB, xCut);
PbB1 = sampleCurveAtX(patchB.S, xCut);

plot(PbA(:,2),  PbA(:,3),  'k-', 'LineWidth', 2.0);
plot(PbB0(:,2), PbB0(:,3), 'r--', 'LineWidth', 1.5);
plot(PbB1(:,2), PbB1(:,3), 'b-', 'LineWidth', 1.5);
xlabel('y'); ylabel('z');
legend({'Patch A cut','Patch B before','Patch B after'}, 'Location','best');

%% 4) Cross-edge derivative mismatch along the join
nexttile; hold on; grid on;
title('Cross-edge tangent mismatch along join');

ts = linspace(0,1,31);
mBefore = zeros(size(ts));
mAfter  = zeros(size(ts));

for k = 1:numel(ts)
    [ua, va] = edgeUV(SA, 'u1', ts(k));
    [ub, vb] = edgeUV(SB, 'u0', ts(k));

    [~, dA_before] = edgeCrossDerivative(SA, 'u1', ua, va);
    [~, dB_before] = edgeCrossDerivative(SB, 'u0', ub, vb);
    [~, dB_after ] = edgeCrossDerivative(patchB.S, 'u0', ub, vb);

    mBefore(k) = mismatchMetric(dA_before, dB_before);
    mAfter(k)  = mismatchMetric(dA_before, dB_after);
end

plot(ts, mBefore, 'r--', 'LineWidth', 1.5);
plot(ts, mAfter,  'b-',  'LineWidth', 1.8);
xlabel('edge parameter');
ylabel('mismatch metric');

fprintf('\nDone.\n');

%% ------------------------------------------------------------------------
function plotSurfaceWire(S, nu, nv, colorRGB)
    u = linspace(S.domainU(1), S.domainU(2), nu);
    v = linspace(S.domainV(1), S.domainV(2), nv);

    for j = 1:numel(v)
        P = zeros(numel(u),3);
        for i = 1:numel(u)
            P(i,:) = S.evaluate(u(i), v(j));
        end
        plot3(P(:,1), P(:,2), P(:,3), '-', 'Color', colorRGB, 'LineWidth', 0.75);
    end
    for i = 1:numel(u)
        P = zeros(numel(v),3);
        for j = 1:numel(v)
            P(j,:) = S.evaluate(u(i), v(j));
        end
        plot3(P(:,1), P(:,2), P(:,3), '-', 'Color', colorRGB, 'LineWidth', 0.75);
    end
end

function plotEdge(S, edgeId, style, lw)
    uv = geom.Patch.edgeSampleUV(S, edgeId, 101);
    P = zeros(size(uv,1),3);
    for ii = 1:size(uv,1)
        P(ii,:) = S.evaluate(uv(ii,1), uv(ii,2));
    end
    plot3(P(:,1), P(:,2), P(:,3), style, 'LineWidth', lw);
end

function plotTangentRibbon(S, edgeId, nPts, scale, colorRGB)
    ts = linspace(0,1,nPts);
    for k = 1:numel(ts)
        [u,v] = edgeUV(S, edgeId, ts(k));
        P = S.evaluate(u,v);
        [~, dCross] = edgeCrossDerivative(S, edgeId, u, v);
        d = dCross / max(norm(dCross), 1e-14);
        Q = P + scale*d;
        plot3([P(1) Q(1)], [P(2) Q(2)], [P(3) Q(3)], '-', ...
            'Color', colorRGB, 'LineWidth', 1.5);
    end
end

function P = sampleCurveAtX(S, xTarget)
    % Robust approximate constant-x cut by scanning u for each v.
    nV = 121;
    nU = 241;
    ts = linspace(S.domainV(1), S.domainV(2), nV);
    us = linspace(S.domainU(1), S.domainU(2), nU);

    P = zeros(nV,3);
    for k = 1:nV
        v = ts(k);
        X = zeros(size(us));
        XYZ = zeros(numel(us),3);
        for i = 1:numel(us)
            XYZ(i,:) = S.evaluate(us(i), v);
            X(i) = XYZ(i,1);
        end
        [~, idx] = min(abs(X - xTarget));
        P(k,:) = XYZ(idx,:);
    end
end

function [u,v] = edgeUV(S, edgeId, t)
    du = S.domainU; dv = S.domainV;
    switch lower(edgeId)
        case 'u0', u = du(1); v = dv(1) + t*(dv(2)-dv(1));
        case 'u1', u = du(2); v = dv(1) + t*(dv(2)-dv(1));
        case 'v0', u = du(1) + t*(du(2)-du(1)); v = dv(1);
        case 'v1', u = du(1) + t*(du(2)-du(1)); v = dv(2);
        otherwise, error('Unknown edge id.');
    end
end

function [tAlong, dCross] = edgeCrossDerivative(S, edgeId, u, v)
    D = S.derivatives(u, v, 1);
    switch lower(edgeId)
        case {'u0','u1'}
            tAlong = toRowVec(squeeze(D(1,2,:)));
            dCross = toRowVec(squeeze(D(2,1,:)));
        case {'v0','v1'}
            tAlong = toRowVec(squeeze(D(2,1,:)));
            dCross = toRowVec(squeeze(D(1,2,:)));
        otherwise
            error('Unknown edge id.');
    end
end

function v = toRowVec(x)
    v = localFlattenToDouble(x);
    v = v(:).';
end

function out = localFlattenToDouble(x)
    if isa(x, 'sym')
        out = double(x);
        return;
    end

    if isnumeric(x) || islogical(x)
        out = double(x);
        return;
    end

    if iscell(x)
        parts = cell(size(x));
        nTotal = 0;
        for ii = 1:numel(x)
            parts{ii} = localFlattenToDouble(x{ii});
            nTotal = nTotal + numel(parts{ii});
        end
        out = zeros(nTotal,1);
        k = 1;
        for ii = 1:numel(parts)
            p = parts{ii}(:);
            out(k:k+numel(p)-1) = p;
            k = k + numel(p);
        end
        return;
    end

    try
        out = double(x);
    catch ME
        error('Could not convert derivative output to numeric form: %s', ME.message);
    end
end

function m = mismatchMetric(a, b)
    na = norm(a); nb = norm(b);
    if na < 1e-14 || nb < 1e-14
        m = NaN;
        return;
    end
    ah = a/na;
    bh = b/nb;
    m = norm(ah + bh);  % opposite vectors -> ~0 for opposing patch normals
end
