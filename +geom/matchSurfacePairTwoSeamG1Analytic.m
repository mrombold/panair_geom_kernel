function [S1out, S2out, info] = matchSurfacePairTwoSeamG1Analytic(S1, S2, varargin)
% geom.matchSurfacePairTwoSeamG1Analytic
% Simultaneous two-seam analytic G1 matcher for two compatible tensor-product
% NURBS surfaces. Intended for cases like symmetric cowling half-patches where
% the same two surfaces share two seam edges (e.g. upper and lower centerline).
%
% This method builds one least-squares solve for BOTH seam edges at once,
% using the existing analytic surface derivatives in the kernel.
%
% Name/value options
%   'EdgePairs'            : 2x2 cellstr of edge names. If omitted, auto-detect
%                            the best two distinct matching seam pairs.
%   'ReverseAlongEdge'     : 2x1 logical or 'auto' (default 'auto')
%   'BoundaryTol'          : shared-edge compatibility tol (default 1e-8)
%   'AverageBoundary'      : average corresponding seam control rows before
%                            solving (default false)
%   'UseAverageMagnitude'  : use average projected derivative magnitude on both
%                            sides at each sample (default true)
%   'MasterScale'          : scale on S1 target magnitude (default 1)
%   'SlaveScale'           : scale on S2 target magnitude (default 1)
%   'LambdaRegularization' : Tikhonov regularization on control movement
%                            (default 1e-10)

pa = inputParser;
addParameter(pa, 'EdgePairs', []);
addParameter(pa, 'ReverseAlongEdge', 'auto');
addParameter(pa, 'BoundaryTol', 1e-8);
addParameter(pa, 'AverageBoundary', false);
addParameter(pa, 'UseAverageMagnitude', true);
addParameter(pa, 'MasterScale', 1.0);
addParameter(pa, 'SlaveScale', 1.0);
addParameter(pa, 'LambdaRegularization', 1e-10);
parse(pa, varargin{:});
opt = pa.Results;

if isempty(opt.EdgePairs)
    pairs = localFindBestTwoEdgePairs(S1, S2);
    edgePairs = {pairs(1).edge1, pairs(1).edge2; pairs(2).edge1, pairs(2).edge2};
    reverseFlags = [pairs(1).reverse; pairs(2).reverse];
else
    edgePairs = opt.EdgePairs;
    if ~iscell(edgePairs) || ~isequal(size(edgePairs), [2 2])
        error('geom:matchSurfacePairTwoSeamG1Analytic', ...
            'EdgePairs must be a 2x2 cell array.');
    end
    if ischar(opt.ReverseAlongEdge) || isstring(opt.ReverseAlongEdge)
        if lower(string(opt.ReverseAlongEdge)) ~= "auto"
            error('geom:matchSurfacePairTwoSeamG1Analytic', ...
                'ReverseAlongEdge must be ''auto'' or a 2x1 logical vector.');
        end
        reverseFlags = localDetectReverseFlags(S1, S2, edgePairs);
    else
        reverseFlags = logical(opt.ReverseAlongEdge(:));
        if numel(reverseFlags) ~= 2
            error('geom:matchSurfacePairTwoSeamG1Analytic', ...
                'ReverseAlongEdge must have two logical values.');
        end
    end
end

pair(1).m = localEdgeAccessors(S1, edgePairs{1,1});
pair(1).s = localEdgeAccessors(S2, edgePairs{1,2});
pair(1).reverse = reverseFlags(1);

pair(2).m = localEdgeAccessors(S1, edgePairs{2,1});
pair(2).s = localEdgeAccessors(S2, edgePairs{2,2});
pair(2).reverse = reverseFlags(2);

for k = 1:2
    B1 = pair(k).m.getB(S1.P);
    B2 = pair(k).s.getB(S2.P);
    if size(B1,1) ~= size(B2,1)
        error('geom:matchSurfacePairTwoSeamG1Analytic', ...
            'Edge pair %d incompatible control counts: %d vs %d.', ...
            k, size(B1,1), size(B2,1));
    end
    if pair(k).reverse
        B2cmp = flipud(B2);
    else
        B2cmp = B2;
    end
    mismatch = max(vecnorm(B1 - B2cmp, 2, 2));
    if ~opt.AverageBoundary && mismatch > opt.BoundaryTol
        error('geom:matchSurfacePairTwoSeamG1Analytic', ...
            'Seam pair %d boundary mismatch exceeds tolerance: %.3e', k, mismatch);
    end
    pair(k).mismatch = mismatch;
    if opt.AverageBoundary
        pair(k).Btarget = 0.5*(B1 + B2cmp);
    else
        pair(k).Btarget = B1;
    end
end

P1 = S1.P;
P2 = S2.P;
for k = 1:2
    P1 = pair(k).m.setB(P1, pair(k).Btarget);
    if pair(k).reverse
        P2 = pair(k).s.setB(P2, flipud(pair(k).Btarget));
    else
        P2 = pair(k).s.setB(P2, pair(k).Btarget);
    end
end

maps = localBuildUnknownMaps(S1, S2, pair);
[A, b, sampleInfo] = localBuildSystem(S1, S2, pair, maps, opt);

lam = opt.LambdaRegularization;
AtA = A.'*A;
rhs = A.'*b;
if lam > 0
    x = (AtA + lam*speye(size(AtA))) \ rhs;
else
    x = AtA \ rhs;
end

P1new = P1;
P2new = P2;
for q = 1:numel(maps.entries)
    ent = maps.entries(q);
    delta = x(3*(q-1)+(1:3)).';
    if ent.surface == 1
        P1new(ent.i, ent.j, :) = reshape(squeeze(P1new(ent.i, ent.j, :)).' + delta, [1 1 3]);
    else
        P2new(ent.i, ent.j, :) = reshape(squeeze(P2new(ent.i, ent.j, :)).' + delta, [1 1 3]);
    end
end

S1out = geom.NURBSSurface(P1new, S1.p, S1.q, S1.U, S1.V, S1.W);
S2out = geom.NURBSSurface(P2new, S2.p, S2.q, S2.U, S2.V, S2.W);

info = struct();
info.edgePairs = edgePairs;
info.reverseAlongEdge = reverseFlags;
info.boundaryMismatch = [pair(1).mismatch; pair(2).mismatch];
info.numUnknownPoints = numel(maps.entries);
info.numScalarUnknowns = size(A,2);
info.numScalarConstraints = size(A,1);
info.residualNorm = norm(A*x - b);
info.sampleInfo = sampleInfo;
info.maps = maps;
end

function [A, b, sampleInfo] = localBuildSystem(S1, S2, pair, maps, opt)
rows = {};
rhs = [];
sampleInfo = struct('seam', {}, 't', {}, 'C', {}, 'T', {}, 'DmProj', {}, 'DsProj', {}, 'dTarget', {}, 'Lm', {}, 'Ls', {});

for k = 1:2
    tVals = localGreville(pair(k).m.alongKnots, pair(k).m.alongDegree, size(pair(k).Btarget,1));
    if pair(k).reverse
        tVals2 = fliplr(tVals);
    else
        tVals2 = tVals;
    end

    for j = 1:numel(tVals)
        t1 = tVals(j);
        t2 = tVals2(j);

        [C1, T1, D1] = localSeamQuantities(S1, pair(k).m.edgeName, t1);
        [C2, T2, D2] = localSeamQuantities(S2, pair(k).s.edgeName, t2);

        T = T1 + T2;
        if norm(T) < 1e-12, T = T1; end
        if norm(T) < 1e-12, T = T2; end
        if norm(T) < 1e-12, T = [1 0 0]; end
        T = T / norm(T);

        D1p = D1 - dot(D1,T)*T;
        D2p = D2 - dot(D2,T)*T;
        L1 = norm(D1p);
        L2 = norm(D2p);

        if L1 < 1e-12 && L2 < 1e-12
            continue;
        end
        if L1 < 1e-12
            u1 = -D2p / max(L2, eps);
        else
            u1 = D1p / L1;
        end
        if L2 < 1e-12
            u2 = -D1p / max(L1, eps);
        else
            u2 = -D2p / L2;
        end

        d = u1 + u2;
        if norm(d) < 1e-12
            if L1 >= L2 && L1 > 1e-12
                d = D1p / L1;
            else
                d = -D2p / max(L2, eps);
            end
        else
            d = d / norm(d);
        end

        if opt.UseAverageMagnitude
            Lm = 0.5*(L1 + L2) * opt.MasterScale;
            Ls = 0.5*(L1 + L2) * opt.SlaveScale;
        else
            Lm = L1 * opt.MasterScale;
            Ls = L2 * opt.SlaveScale;
        end

        D1t =  Lm * d;
        D2t = -Ls * d;

        r1 = zeros(3, 3*numel(maps.entries));
        r2 = zeros(3, 3*numel(maps.entries));

        [i1, j1] = pair(k).m.interiorSubscripts(j, size(pair(k).Btarget,1), size(S1.P));
        idx1 = maps.lookup([1, i1, j1]);

        jSlaveAlong = localSlaveIndex(j, size(pair(k).Btarget,1), pair(k).reverse);
        [i2, j2] = pair(k).s.interiorSubscripts(jSlaveAlong, size(pair(k).Btarget,1), size(S2.P));
        idx2 = maps.lookup([2, i2, j2]);

        cols1 = 3*(idx1-1)+(1:3);
        cols2 = 3*(idx2-1)+(1:3);

        r1(:, cols1) = pair(k).m.inwardFactor * eye(3);
        r2(:, cols2) = pair(k).s.inwardFactor * eye(3);

        rows{end+1} = r1; %#ok<AGROW>
        rhs = [rhs; (D1t - D1).']; %#ok<AGROW>
        rows{end+1} = r2; %#ok<AGROW>
        rhs = [rhs; (D2t - D2).']; %#ok<AGROW>

        si.seam = k;
        si.t = t1;
        si.C = 0.5*(C1 + C2);
        si.T = T;
        si.DmProj = D1p;
        si.DsProj = D2p;
        si.dTarget = d;
        si.Lm = Lm;
        si.Ls = Ls;
        sampleInfo(end+1) = si; %#ok<AGROW>
    end
end

if isempty(rows)
    A = zeros(0, 3*numel(maps.entries));
    b = zeros(0,1);
else
    A = vertcat(rows{:});
    b = rhs;
end
end

function jj = localSlaveIndex(j, n, reverseFlag)
if reverseFlag
    jj = n - j + 1;
else
    jj = j;
end
end

function maps = localBuildUnknownMaps(S1, S2, pair)
entries = struct('surface', {}, 'i', {}, 'j', {});
lookupMap = containers.Map();

for k = 1:2
    nAlong = size(pair(k).Btarget,1);
    for j = 1:nAlong
        [i1, j1] = pair(k).m.interiorSubscripts(j, nAlong, size(S1.P));
        key1 = sprintf('1_%d_%d', i1, j1);
        if ~isKey(lookupMap, key1)
            entries(end+1).surface = 1; %#ok<AGROW>
            entries(end).i = i1;
            entries(end).j = j1;
            lookupMap(key1) = numel(entries);
        end

        jSlaveAlong = localSlaveIndex(j, nAlong, pair(k).reverse);
        [i2, j2] = pair(k).s.interiorSubscripts(jSlaveAlong, nAlong, size(S2.P));
        key2 = sprintf('2_%d_%d', i2, j2);
        if ~isKey(lookupMap, key2)
            entries(end+1).surface = 2; %#ok<AGROW>
            entries(end).i = i2;
            entries(end).j = j2;
            lookupMap(key2) = numel(entries);
        end
    end
end

maps = struct();
maps.entries = entries;
maps.lookupMap = lookupMap;
maps.lookup = @(sub) lookupMap(sprintf('%d_%d_%d', sub(1), sub(2), sub(3)));
end

function [C, T, Din] = localSeamQuantities(S, edgeName, t)
edgeName = lower(string(edgeName));
switch edgeName
    case "u0"
        u = S.domainU(1); v = t;
        D = S.derivatives(u, v, 1);
        C = D{1,1}; Su = D{2,1}; Sv = D{1,2};
        T = Sv; Din = Su;
    case "u1"
        u = S.domainU(2); v = t;
        D = S.derivatives(u, v, 1);
        C = D{1,1}; Su = D{2,1}; Sv = D{1,2};
        T = Sv; Din = -Su;
    case "v0"
        u = t; v = S.domainV(1);
        D = S.derivatives(u, v, 1);
        C = D{1,1}; Su = D{2,1}; Sv = D{1,2};
        T = Su; Din = Sv;
    case "v1"
        u = t; v = S.domainV(2);
        D = S.derivatives(u, v, 1);
        C = D{1,1}; Su = D{2,1}; Sv = D{1,2};
        T = Su; Din = -Sv;
    otherwise
        error('geom:matchSurfacePairTwoSeamG1Analytic', 'Unknown edge "%s".', char(edgeName));
end
end

function t = localGreville(K, p, nCtrl)
t = zeros(1, nCtrl);
for i = 1:nCtrl
    if p == 0
        t(i) = K(i+1);
    else
        t(i) = mean(K(i+1:i+p));
    end
end
end

function flags = localDetectReverseFlags(S1, S2, edgePairs)
flags = false(2,1);
for k = 1:2
    B1 = localGetEdgeCP(S1, edgePairs{k,1});
    B2 = localGetEdgeCP(S2, edgePairs{k,2});
    sameErr = max(vecnorm(B1 - B2, 2, 2));
    revErr  = max(vecnorm(B1 - flipud(B2), 2, 2));
    flags(k) = revErr < sameErr;
end
end

function pairs = localFindBestTwoEdgePairs(S1, S2)
edges = {'u0','u1','v0','v1'};
cand = struct('edge1', {}, 'edge2', {}, 'reverse', {}, 'mismatch', {});
for i = 1:numel(edges)
    for j = 1:numel(edges)
        B1 = localGetEdgeCP(S1, edges{i});
        B2 = localGetEdgeCP(S2, edges{j});
        if size(B1,1) ~= size(B2,1), continue; end
        sameErr = max(vecnorm(B1 - B2, 2, 2));
        revErr  = max(vecnorm(B1 - flipud(B2), 2, 2));
        c.edge1 = edges{i};
        c.edge2 = edges{j};
        c.reverse = revErr < sameErr;
        c.mismatch = min(sameErr, revErr);
        cand(end+1) = c; %#ok<AGROW>
    end
end
[~, order] = sort([cand.mismatch], 'ascend');
cand = cand(order);
pairs = cand(1);
used1 = string(cand(1).edge1);
used2 = string(cand(1).edge2);
for k = 2:numel(cand)
    if string(cand(k).edge1) ~= used1 && string(cand(k).edge2) ~= used2
        pairs(2) = cand(k); %#ok<AGROW>
        return;
    end
end
error('Could not identify two distinct seam edge pairs automatically.');
end

function B = localGetEdgeCP(S, edgeName)
switch lower(edgeName)
    case 'u0'
        B = squeeze(S.P(1,:,:));
    case 'u1'
        B = squeeze(S.P(end,:,:));
    case 'v0'
        B = squeeze(S.P(:,1,:));
    case 'v1'
        B = squeeze(S.P(:,end,:));
    otherwise
        error('Unknown edge');
end
end

function info = localEdgeAccessors(S, edgeName)
edgeName = lower(string(edgeName));
info.edgeName = char(edgeName);
switch edgeName
    case "u0"
        info.getB = @(P) squeeze(P(1,:,:));
        info.getI = @(P) squeeze(P(2,:,:));
        info.setB = @(Pin, B) localSetRow(Pin, 1, B);
        info.setI = @(Pin, I) localSetRow(Pin, 2, I);
        info.inwardFactor = S.p / (S.U(S.p+2) - S.U(2));
        info.alongDegree = S.q;
        info.alongKnots = S.V;
        info.interiorSubscripts = @(j, nAlong, szP) deal(2, j);
    case "u1"
        info.getB = @(P) squeeze(P(end,:,:));
        info.getI = @(P) squeeze(P(end-1,:,:));
        info.setB = @(Pin, B) localSetRow(Pin, size(Pin,1), B);
        info.setI = @(Pin, I) localSetRow(Pin, size(Pin,1)-1, I);
        info.inwardFactor = S.p / (S.U(end) - S.U(end-S.p-1));
        info.alongDegree = S.q;
        info.alongKnots = S.V;
        info.interiorSubscripts = @(j, nAlong, szP) deal(szP(1)-1, j);
    case "v0"
        info.getB = @(P) squeeze(P(:,1,:));
        info.getI = @(P) squeeze(P(:,2,:));
        info.setB = @(Pin, B) localSetCol(Pin, 1, B);
        info.setI = @(Pin, I) localSetCol(Pin, 2, I);
        info.inwardFactor = S.q / (S.V(S.q+2) - S.V(2));
        info.alongDegree = S.p;
        info.alongKnots = S.U;
        info.interiorSubscripts = @(j, nAlong, szP) deal(j, 2);
    case "v1"
        info.getB = @(P) squeeze(P(:,end,:));
        info.getI = @(P) squeeze(P(:,end-1,:));
        info.setB = @(Pin, B) localSetCol(Pin, size(Pin,2), B);
        info.setI = @(Pin, I) localSetCol(Pin, size(Pin,2)-1, I);
        info.inwardFactor = S.q / (S.V(end) - S.V(end-S.q-1));
        info.alongDegree = S.p;
        info.alongKnots = S.U;
        info.interiorSubscripts = @(j, nAlong, szP) deal(j, szP(2)-1);
    otherwise
        error('geom:matchSurfacePairTwoSeamG1Analytic', 'Unknown edge "%s".', char(edgeName));
end
end

function P = localSetRow(P, i, B)
P(i,:,:) = reshape(B, [1 size(B,1) 3]);
end

function P = localSetCol(P, j, B)
P(:,j,:) = reshape(B, [size(B,1) 1 3]);
end
