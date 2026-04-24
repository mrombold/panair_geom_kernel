
function [Csec, info] = conicSurfaceSection(S, planeArg, arg2, varargin)
% geom.conicSurfaceSection
% Construct explicit section curve(s) from a geom.ConicSurface and a plane.
%
% Supported plane inputs:
%   [Csec,info] = geom.conicSurfaceSection(S, 'x', x0, ...)
%   [Csec,info] = geom.conicSurfaceSection(S, 'y', y0, ...)
%   [Csec,info] = geom.conicSurfaceSection(S, 'z', z0, ...)
%   [Csec,info] = geom.conicSurfaceSection(S, point, normal, ...)
%
% Name/value options
%   'Offset'            : extra offset along normal (default 0)
%   'SamplesU'          : number of sampled surface sections (default 61)
%   'RefinePasses'      : midpoint refinement in u (default 1)
%   'Tol'               : intersection tolerance (default 1e-9)
%   'MaxIter'           : max bisection iterations (default 60)
%   'Degree'            : fitted output degree (default min(3,nPts-1))
%   'Method'            : interpolation parameterization (default 'centripetal')
%   'MinBranchPoints'   : minimum points to keep a branch (default 3)
%   'BranchCountAction' : 'warn' | 'error' | 'split' (default 'warn')
%   'MultiBranchAction' : 'warn_primary' | 'error' | 'all' (default 'warn_primary')
%   'ExtendToGuides'    : if true, extend branch endpoints to upper/lower guide
%                         intersections when present (default true)
%   'GuideSamplesU'     : samples used to find guide-plane intersections
%                         (default max(201, 2*SamplesU+1))
%   'EndpointUTol'      : allowable u-gap when snapping to guide endpoints
%                         (default 3*mean sampled du)
%   'EndpointDistTol'   : allowable spatial gap when snapping to guides
%                         (default 2%% of surface bbox diagonal)
%
% Outputs
%   Csec : geom.NURBSCurve or cell array of curves when MultiBranchAction='all'
%   info : diagnostics struct

pa = inputParser;
addParameter(pa, 'Offset', 0.0);
addParameter(pa, 'SamplesU', 61);
addParameter(pa, 'RefinePasses', 1);
addParameter(pa, 'Tol', 1e-9);
addParameter(pa, 'MaxIter', 60);
addParameter(pa, 'Degree', []);
addParameter(pa, 'Method', 'centripetal');
addParameter(pa, 'MinBranchPoints', 3);
addParameter(pa, 'BranchCountAction', 'warn');
addParameter(pa, 'MultiBranchAction', 'warn_primary');
addParameter(pa, 'ExtendToGuides', true);
addParameter(pa, 'GuideSamplesU', []);
addParameter(pa, 'EndpointUTol', []);
addParameter(pa, 'EndpointDistTol', []);
parse(pa, varargin{:});
opt = pa.Results;

[planePoint, planeNormal, planeDesc] = localParsePlane(planeArg, arg2, opt.Offset);
planeNormal = planeNormal / norm(planeNormal);

[uVals, sections] = localSampleSections(S, opt.SamplesU, opt.RefinePasses);
surfaceScale = localSurfaceScale(S, sections, uVals);

[result, ok] = localBuildBranches(sections, uVals, planePoint, planeNormal, opt, surfaceScale);
if ~ok
    error('geom:conicSurfaceSection', ...
        'Could not build a usable section curve for the requested plane.');
end

if opt.ExtendToGuides
    guideOpt = opt;
    if isempty(guideOpt.GuideSamplesU)
        guideOpt.GuideSamplesU = max(201, 2*opt.SamplesU + 1);
    end
    if isempty(guideOpt.EndpointUTol)
        du = diff(uVals);
        if isempty(du)
            guideOpt.EndpointUTol = 1e-6;
        else
            guideOpt.EndpointUTol = 3 * mean(abs(du));
        end
    end
    if isempty(guideOpt.EndpointDistTol)
        guideOpt.EndpointDistTol = 0.02 * max(surfaceScale, 1.0);
    end

    guideHits = localFindGuideBoundaryHits(S, planePoint, planeNormal, guideOpt);
    [result, guideAttachInfo] = localExtendBranchesToGuides(result, guideHits, guideOpt);
else
    guideHits = struct('side', {}, 'u', {}, 't', {}, 'pt', {});
    guideAttachInfo = struct('numAttached', 0, 'numGuideHits', 0, 'attached', []);
end

[result.curves, result.fitRMS, result.numPtsEach] = localFitCurvesFromBranches( ...
    result.branchPoints, result.branchU, opt);

info = rmfield(result, 'curves');
info.planePoint = planePoint;
info.planeNormal = planeNormal;
info.planeDescription = planeDesc;
info.uSamples = uVals;
info.surfaceScale = surfaceScale;
info.guideHits = guideHits;
info.guideAttachment = guideAttachInfo;

nBranches = numel(result.curves);
info.numBranches = nBranches;

action = lower(string(opt.MultiBranchAction));
switch action
    case "warn_primary"
        if nBranches > 1
            warning('geom:conicSurfaceSection:MultipleBranches', ...
                'Plane produced %d section branches. Returning the primary branch only.', nBranches);
        end
        Csec = result.curves{1};
    case "error"
        if nBranches > 1
            error('geom:conicSurfaceSection:MultipleBranches', ...
                'Plane produced %d section branches. Use ''MultiBranchAction'',''all'' to retrieve all.', nBranches);
        end
        Csec = result.curves{1};
    case "all"
        if nBranches == 1
            Csec = result.curves{1};
        else
            Csec = result.curves;
        end
    otherwise
        error('geom:conicSurfaceSection', ...
            'Unknown MultiBranchAction "%s".', char(action));
end

end

function [uVals, sections] = localSampleSections(S, nU, refinePasses)
uVals = linspace(S.uDomain(1), S.uDomain(2), max(2, round(nU)));
for pass = 1:refinePasses
    mids = 0.5 * (uVals(1:end-1) + uVals(2:end));
    uVals = sort(unique([uVals, mids]));
end
uVals = uVals(:);

sections = cell(numel(uVals),1);
for i = 1:numel(uVals)
    sections{i} = S.sectionAt(uVals(i));
end
end

function scale = localSurfaceScale(S, sections, uVals)
pts = zeros(2*numel(uVals), 3);
row = 1;
for i = 1:numel(uVals)
    C = sections{i};
    dom = C.domain;
    pts(row,:)   = C.evaluate(dom(1));
    pts(row+1,:) = C.evaluate(dom(2));
    row = row + 2;
end
mins = min(pts, [], 1);
maxs = max(pts, [], 1);
scale = norm(maxs - mins);
if scale <= 0
    scale = 1.0;
end
end

function [result, ok] = localBuildBranches(sections, uVals, planePoint, planeNormal, opt, surfaceScale)
result = struct();
ok = false;

samples = struct('u', {}, 'pts', {}, 'params', {});
for k = 1:numel(uVals)
    C = sections{k};
    [ptsK, tK] = localIntersectCurvePlaneBisection(C, planePoint, planeNormal, ...
        'Tol', opt.Tol, 'MaxIter', opt.MaxIter);

    if isempty(tK)
        continue;
    end

    [tK, ord] = sort(tK(:));
    ptsK = ptsK(ord,:);

    samples(end+1).u = uVals(k); %#ok<AGROW>
    samples(end).pts = ptsK;
    samples(end).params = tK;
end

if isempty(samples)
    return;
end

counts = arrayfun(@(s) size(s.pts,1), samples);
countChangeIdx = find(diff(counts) ~= 0);

if ~isempty(countChangeIdx)
    msg = sprintf(['Intersection count changed across sampled conic sections at %d location(s). ' ...
        'Counts = [%s].'], numel(countChangeIdx), sprintf('%d ', counts));
    switch lower(string(opt.BranchCountAction))
        case "error"
            error('geom:conicSurfaceSection:ChangingIntersectionCount', '%s', msg);
        case "warn"
            warning('geom:conicSurfaceSection:ChangingIntersectionCount', '%s', msg);
        case "split"
            % split handled below
        otherwise
            error('geom:conicSurfaceSection', ...
                'Unknown BranchCountAction "%s".', char(opt.BranchCountAction));
    end
end

segments = localConstantCountSegments(counts);

branchesPts = {};
branchesU = {};
branchesT = {};

for s = 1:size(segments,1)
    i1 = segments(s,1);
    i2 = segments(s,2);
    c = counts(i1);

    if c <= 0
        continue;
    end

    for j = 1:c
        pts = zeros(i2-i1+1, 3);
        uu = zeros(i2-i1+1, 1);
        tt = zeros(i2-i1+1, 1);
        row = 1;
        for k = i1:i2
            pts(row,:) = samples(k).pts(j,:);
            uu(row,1) = samples(k).u;
            tt(row,1) = samples(k).params(j);
            row = row + 1;
        end

        [pts, uu, tt] = localDedupConsecutive(pts, uu, tt, 1e-9);

        if size(pts,1) >= opt.MinBranchPoints
            branchesPts{end+1,1} = pts; %#ok<AGROW>
            branchesU{end+1,1} = uu; %#ok<AGROW>
            branchesT{end+1,1} = tt; %#ok<AGROW>
        end
    end
end

if isempty(branchesPts)
    return;
end

[curves, fitRMS, numPtsEach, order] = localFitCurvesFromBranches(branchesPts, branchesU, opt);

result.curves = curves;
result.branchPoints = branchesPts(order);
result.branchU = branchesU(order);
result.branchSectionParams = branchesT(order);
result.numPts = sum(numPtsEach);
result.fitRMS = min(fitRMS);
result.numBranches = numel(curves);
result.intersectionCounts = counts;
result.countChangeIndices = countChangeIdx;
result.surfaceScale = surfaceScale;
ok = true;
end

function [pts, uu, tt] = localDedupConsecutive(pts, uu, tt, tol)
if size(pts,1) < 2
    return;
end
keep = true(size(pts,1),1);
for r = 2:size(pts,1)
    if norm(pts(r,:) - pts(r-1,:)) <= tol
        keep(r) = false;
    end
end
pts = pts(keep,:);
uu = uu(keep);
tt = tt(keep);
end

function segments = localConstantCountSegments(counts)
if isempty(counts)
    segments = zeros(0,2);
    return;
end
startIdx = 1;
segments = zeros(0,2);
for k = 2:numel(counts)
    if counts(k) ~= counts(k-1)
        segments(end+1,:) = [startIdx, k-1]; %#ok<AGROW>
        startIdx = k;
    end
end
segments(end+1,:) = [startIdx, numel(counts)];
end

function guideHits = localFindGuideBoundaryHits(S, planePoint, planeNormal, opt)
guideHits = struct('side', {}, 'u', {}, 't', {}, 'pt', {});
uVals = linspace(S.uDomain(1), S.uDomain(2), opt.GuideSamplesU).';

for sideCell = {'upper', 'lower'}
    side = sideCell{1};
    g = zeros(numel(uVals), 1);
    for k = 1:numel(uVals)
        [pt, ~] = localBoundaryGuidePoint(S, uVals(k), side);
        g(k) = dot(pt - planePoint, planeNormal);
    end

    rootU = localFindRoots1D(uVals, g, opt.Tol, opt.MaxIter);

    for i = 1:numel(rootU)
        [pt, tBnd] = localBoundaryGuidePoint(S, rootU(i), side);
        hit.side = side;
        hit.u = rootU(i);
        hit.t = tBnd;
        hit.pt = pt;
        guideHits(end+1,1) = hit; %#ok<AGROW>
    end
end

% Deduplicate coincident hits (e.g. tangency or both guides touching the plane)
if isempty(guideHits)
    return;
end

keep = true(numel(guideHits),1);
for i = 2:numel(guideHits)
    for j = 1:i-1
        if keep(j) && norm(guideHits(i).pt - guideHits(j).pt) <= max(1e-8, 10*opt.Tol)
            if abs(guideHits(i).u - guideHits(j).u) <= max(1e-8, 10*opt.Tol)
                keep(i) = false;
                break;
            end
        end
    end
end
guideHits = guideHits(keep);
end

function [pt, tBnd] = localBoundaryGuidePoint(S, u, side)
C = S.sectionAt(u);
dom = C.domain;
if strcmpi(side, 'upper')
    tBnd = dom(1);
else
    tBnd = dom(2);
end
pt = C.evaluate(tBnd);
end

function roots = localFindRoots1D(uVals, fVals, tol, maxIter)
roots = [];

for k = 1:numel(uVals)-1
    u0 = uVals(k);   u1 = uVals(k+1);
    f0 = fVals(k);   f1 = fVals(k+1);

    if abs(f0) <= tol
        roots(end+1,1) = u0; %#ok<AGROW>
    end

    if f0 == 0
        continue;
    end

    if f0 * f1 < 0
        a = u0; b = u1; fa = f0; fb = f1;
        for it = 1:maxIter
            m = 0.5 * (a + b);
            [fm] = interp1([a b],[fa fb],m,'linear');
            % More robust: direct evaluation from linearized interval is enough here
            % because the bracket is already tiny after the final iterations.
            % The caller recomputes the actual point at the final root.
            if abs(fm) <= tol || abs(b-a) <= tol
                a = m; b = m;
                break;
            end
            if fa * fm <= 0
                b = m; fb = fm;
            else
                a = m; fa = fm;
            end
        end
        roots(end+1,1) = 0.5 * (a + b); %#ok<AGROW>
    elseif abs(f1) <= tol
        roots(end+1,1) = u1; %#ok<AGROW>
    end
end

if isempty(roots)
    return;
end

roots = sort(roots(:));
keep = true(size(roots));
for k = 2:numel(roots)
    if abs(roots(k) - roots(k-1)) <= max(1e-8, tol)
        keep(k) = false;
    end
end
roots = roots(keep);
end

function [result, attachInfo] = localExtendBranchesToGuides(result, guideHits, opt)
attachInfo = struct('numAttached', 0, 'numGuideHits', numel(guideHits), 'attached', []);

if isempty(guideHits) || isempty(result.branchPoints)
    return;
end

attached = [];

for h = 1:numel(guideHits)
    hit = guideHits(h);
    bestBranch = [];
    bestEnd = '';
    bestScore = inf;

    for b = 1:numel(result.branchPoints)
        uu = result.branchU{b};
        pts = result.branchPoints{b};
        uStart = uu(1);
        uEnd   = uu(end);

        % Attach only to plausible endpoint based on u-order.
        if hit.u <= uStart + opt.EndpointUTol
            d = norm(hit.pt - pts(1,:));
            score = d + 0.1 * abs(hit.u - uStart);
            if score < bestScore
                bestScore = score;
                bestBranch = b;
                bestEnd = 'start';
            end
        end

        if hit.u >= uEnd - opt.EndpointUTol
            d = norm(hit.pt - pts(end,:));
            score = d + 0.1 * abs(hit.u - uEnd);
            if score < bestScore
                bestScore = score;
                bestBranch = b;
                bestEnd = 'end';
            end
        end
    end

    if isempty(bestBranch)
        continue;
    end

    if bestScore > opt.EndpointDistTol + 0.1 * opt.EndpointUTol
        continue;
    end

    pts = result.branchPoints{bestBranch};
    uu  = result.branchU{bestBranch};
    tt  = result.branchSectionParams{bestBranch};

    switch bestEnd
        case 'start'
            if norm(hit.pt - pts(1,:)) > 1e-9
                pts = [hit.pt; pts];
                uu  = [hit.u; uu];
                tt  = [hit.t; tt];
                [uu, idx] = sort(uu, 'ascend');
                pts = pts(idx,:);
                tt = tt(idx);
            end
        case 'end'
            if norm(hit.pt - pts(end,:)) > 1e-9
                pts = [pts; hit.pt];
                uu  = [uu; hit.u];
                tt  = [tt; hit.t];
                [uu, idx] = sort(uu, 'ascend');
                pts = pts(idx,:);
                tt = tt(idx);
            end
    end

    [pts, uu, tt] = localDedupConsecutive(pts, uu, tt, 1e-9);

    result.branchPoints{bestBranch} = pts;
    result.branchU{bestBranch} = uu;
    result.branchSectionParams{bestBranch} = tt;

    rec.branch = bestBranch;
    rec.endpoint = bestEnd;
    rec.side = hit.side;
    rec.u = hit.u;
    rec.point = hit.pt;
    attached = [attached; rec]; %#ok<AGROW>
end

attachInfo.numAttached = numel(attached);
attachInfo.attached = attached;
end

function [curves, fitRMS, numPtsEach, order] = localFitCurvesFromBranches(branchPoints, branchU, opt)
nB = numel(branchPoints);
curves = cell(nB,1);
fitRMS = inf(nB,1);
numPtsEach = zeros(nB,1);

for i = 1:nB
    pts = branchPoints{i};

    deg = opt.Degree;
    if isempty(deg)
        deg = min(3, size(pts,1)-1);
    end
    deg = max(1, min(deg, size(pts,1)-1));

    Cfit = geom.NURBSCurve.globalInterp(pts, deg, opt.Method);
    uData = geom.NURBSCurve.parameterizeData(pts, opt.Method, deg);
    fitPts = Cfit.evaluate(uData);
    fitRMS(i) = sqrt(mean(sum((fitPts - pts).^2, 2)));

    curves{i} = Cfit;
    numPtsEach(i) = size(pts,1);
end

[~, order] = sortrows([-numPtsEach(:), fitRMS(:)], [1 2]);

curves = curves(order);
fitRMS = fitRMS(order);
numPtsEach = numPtsEach(order);
end

function [pts, uvals, info] = localIntersectCurvePlaneBisection(C, planePoint, planeNormal, varargin)
pa = inputParser;
addParameter(pa, 'Samples', []);
addParameter(pa, 'Tol', 1e-9);
addParameter(pa, 'MaxIter', 60);
addParameter(pa, 'DedupTol', 1e-8);
parse(pa, varargin{:});
opts = pa.Results;

n = planeNormal(:).';
nn = norm(n);
if nn < eps
    error('geom:conicSurfaceSection', 'Plane normal must be nonzero.');
end
n = n / nn;
p0 = planePoint(:).';

if isempty(opts.Samples)
    uniqueKnots = unique(C.U);
    opts.Samples = max(41, 8 * max(numel(uniqueKnots)-1, 1) + 1);
end

uMin = C.domain(1);
uMax = C.domain(2);
us = linspace(uMin, uMax, opts.Samples).';
Xs = C.evaluate(us);
fs = (Xs - p0) * n(:);

brackets = zeros(0,2);
exactHits = [];

for k = 1:numel(us)-1
    f0 = fs(k);
    f1 = fs(k+1);

    if abs(f0) <= opts.Tol
        exactHits(end+1,1) = us(k); %#ok<AGROW>
    end

    if f0 == 0
        continue;
    end

    if f0 * f1 < 0
        brackets(end+1,:) = [us(k), us(k+1)]; %#ok<AGROW>
    elseif abs(f1) <= opts.Tol
        exactHits(end+1,1) = us(k+1); %#ok<AGROW>
    end
end

uRoots = exactHits(:);

for b = 1:size(brackets,1)
    a = brackets(b,1);
    c = brackets(b,2);
    fa = localSignedDistance(C, a, p0, n);
    fc = localSignedDistance(C, c, p0, n);

    if fa * fc > 0
        continue;
    end

    for it = 1:opts.MaxIter
        m = 0.5 * (a + c);
        fm = localSignedDistance(C, m, p0, n);

        if abs(fm) <= opts.Tol || abs(c-a) <= opts.DedupTol
            a = m;
            c = m;
            break;
        end

        if fa * fm <= 0
            c = m;
            fc = fm;
        else
            a = m;
            fa = fm;
        end
    end

    uRoots(end+1,1) = 0.5 * (a + c); %#ok<AGROW>
end

if isempty(uRoots)
    pts = zeros(0,3);
    uvals = zeros(0,1);
    info = struct('numBrackets', size(brackets,1), 'sampleU', us, 'sampleSignedDistance', fs);
    return;
end

uRoots = sort(uRoots);
keep = true(size(uRoots));
for k = 2:numel(uRoots)
    if abs(uRoots(k) - uRoots(k-1)) <= opts.DedupTol
        keep(k) = false;
    end
end
uvals = uRoots(keep);
pts = C.evaluate(uvals);

info = struct();
info.numBrackets = size(brackets,1);
info.sampleU = us;
info.sampleSignedDistance = fs;
info.brackets = brackets;
end

function f = localSignedDistance(C, u, p0, n)
X = C.evaluate(u);
f = dot(X - p0, n);
end

function [planePoint, planeNormal, desc] = localParsePlane(planeArg, arg2, offset)
if (ischar(planeArg) || isstring(planeArg)) && isscalar(planeArg)
    ax = lower(char(planeArg));
    val = arg2;
    switch ax
        case 'x'
            planeNormal = [1 0 0];
            planePoint = [val 0 0];
            desc = sprintf('X = %.6g', val);
        case 'y'
            planeNormal = [0 1 0];
            planePoint = [0 val 0];
            desc = sprintf('Y = %.6g', val);
        case 'z'
            planeNormal = [0 0 1];
            planePoint = [0 0 val];
            desc = sprintf('Z = %.6g', val);
        otherwise
            error('geom:conicSurfaceSection', 'Primary plane must be ''x'', ''y'', or ''z''.');
    end
    if abs(offset) > 0
        planePoint = planePoint + offset * planeNormal;
    end
else
    planePoint = planeArg(:).';
    planeNormal = arg2(:).';
    if numel(planePoint) ~= 3 || numel(planeNormal) ~= 3
        error('geom:conicSurfaceSection', ...
            'Arbitrary plane form requires point and normal as 1x3 vectors.');
    end
    nn = norm(planeNormal);
    if nn < eps
        error('geom:conicSurfaceSection', 'Plane normal must be nonzero.');
    end
    planeNormal = planeNormal / nn;
    planePoint = planePoint + offset * planeNormal;
    desc = sprintf('point [%g %g %g], normal [%g %g %g], offset %g', ...
        planePoint(1), planePoint(2), planePoint(3), ...
        planeNormal(1), planeNormal(2), planeNormal(3), offset);
end
end
