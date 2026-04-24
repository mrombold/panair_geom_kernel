function [Csec, info] = surfacePlaneSection_multibranch(S, planeArg, arg2, varargin)
% geom.surfacePlaneSection_multibranch
% Construct explicit NURBS section curve(s) from a NURBS surface and a plane.
%
% This version uses ORDERED intersection tracking instead of nearest-neighbor
% spatial stitching, which avoids false bridge curves when a plane intersects
% the sampled iso-curves multiple times.
%
% Supported plane inputs:
%   [Csec,info] = geom.surfacePlaneSection_multibranch(S, 'x', x0, ...)
%   [Csec,info] = geom.surfacePlaneSection_multibranch(S, 'y', y0, ...)
%   [Csec,info] = geom.surfacePlaneSection_multibranch(S, 'z', z0, ...)
%   [Csec,info] = geom.surfacePlaneSection_multibranch(S, point, normal, ...)
%
% Name-value options
%   'Offset'            : extra offset along normal (default 0)
%   'Direction'         : 'auto' | 'u' | 'v' (default 'auto')
%   'FixedValues'       : 'knot' | 'uniform' | numeric vector (default 'knot')
%   'RefinePasses'      : midpoint refinement count (default 1)
%   'Tol'               : plane/curve bisection tolerance (default 1e-9)
%   'MaxIter'           : max bisection iterations (default 60)
%   'Degree'            : output curve degree (default min(3,nPts-1))
%   'Method'            : interpolation parameterization (default 'centripetal')
%   'MinBranchPoints'   : minimum points to keep a branch (default 3)
%   'BranchCountAction' : 'warn' | 'error' | 'split'
%                         behavior when the number of intersections changes
%                         from one iso-curve sample to the next (default 'warn')
%   'MultiBranchAction' : 'warn_primary' | 'error' | 'all'
%                         (default 'warn_primary')
%
% Outputs
%   Csec : geom.NURBSCurve or cell array of curves when MultiBranchAction='all'
%   info : diagnostics struct

pa = inputParser;
addParameter(pa, 'Offset', 0.0);
addParameter(pa, 'Direction', 'auto');
addParameter(pa, 'FixedValues', 'knot');
addParameter(pa, 'RefinePasses', 1);
addParameter(pa, 'Tol', 1e-9);
addParameter(pa, 'MaxIter', 60);
addParameter(pa, 'Degree', []);
addParameter(pa, 'Method', 'centripetal');
addParameter(pa, 'MinBranchPoints', 3);
addParameter(pa, 'BranchCountAction', 'warn');
addParameter(pa, 'MultiBranchAction', 'warn_primary');
parse(pa, varargin{:});
opts = pa.Results;

[planePoint, planeNormal, planeDesc] = localParsePlane(planeArg, arg2, opts.Offset);
planeNormal = planeNormal / norm(planeNormal);

dirMode = lower(string(opts.Direction));

if dirMode == "auto"
    candidates = localTryDirectionSet(S, planePoint, planeNormal, 'u', opts);
    candidates = [candidates; localTryDirectionSet(S, planePoint, planeNormal, 'v', opts)]; %#ok<AGROW>
else
    candidates = localTryDirectionSet(S, planePoint, planeNormal, char(dirMode), opts);
    if isempty(candidates)
        if char(dirMode) == 'u'
            candidates = localTryDirectionSet(S, planePoint, planeNormal, 'v', opts);
        else
            candidates = localTryDirectionSet(S, planePoint, planeNormal, 'u', opts);
        end
    end
end

if isempty(candidates)
    error('geom:surfacePlaneSection_multibranch', ...
        'Could not build a usable section curve in requested direction "%s".', char(dirMode));
end

pick = localPickBest(candidates);

info = rmfield(pick, 'curves');
info.planePoint = planePoint;
info.planeNormal = planeNormal;
info.planeDescription = planeDesc;

nBranches = numel(pick.curves);
info.numBranches = nBranches;

action = lower(string(opts.MultiBranchAction));
switch action
    case "warn_primary"
        if nBranches > 1
            warning('geom:surfacePlaneSection_multibranch:MultipleBranches', ...
                'Plane produced %d section branches. Returning the primary branch only.', nBranches);
        end
        Csec = pick.curves{1};
    case "error"
        if nBranches > 1
            error('geom:surfacePlaneSection_multibranch:MultipleBranches', ...
                'Plane produced %d section branches. Use ''MultiBranchAction'',''all'' to retrieve all.', nBranches);
        end
        Csec = pick.curves{1};
    case "all"
        if nBranches == 1
            Csec = pick.curves{1};
        else
            Csec = pick.curves;
        end
    otherwise
        error('geom:surfacePlaneSection_multibranch', ...
            'Unknown MultiBranchAction "%s".', char(action));
end

end

function candidates = localTryDirectionSet(S, planePoint, planeNormal, dirFlag, opts)
candidates = struct('curves', {}, 'direction', {}, 'fixedValuesRequested', {}, ...
    'fixedValuesUsed', {}, 'branchPoints', {}, 'branchIsoParams', {}, ...
    'branchFixedValues', {}, 'numPts', {}, 'fitRMS', {}, 'strategy', {}, ...
    'numBranches', {}, 'intersectionCounts', {}, 'countChangeIndices', {});

tries = { ...
    struct('FixedValues', opts.FixedValues, 'RefinePasses', opts.RefinePasses, 'strategy', 'requested'), ...
    struct('FixedValues', 'knot',    'RefinePasses', max(opts.RefinePasses, 2), 'strategy', 'knot_refined'), ...
    struct('FixedValues', 'uniform', 'RefinePasses', max(opts.RefinePasses, 2), 'strategy', 'uniform_refined'), ...
    struct('FixedValues', localDenseValues(S, dirFlag, 31), 'RefinePasses', 0, 'strategy', 'dense_uniform_31'), ...
    struct('FixedValues', localDenseValues(S, dirFlag, 61), 'RefinePasses', 0, 'strategy', 'dense_uniform_61') ...
    };

for t = 1:numel(tries)
    optsTry = opts;
    optsTry.FixedValues = tries{t}.FixedValues;
    optsTry.RefinePasses = tries{t}.RefinePasses;

    [result, ok] = localBuildSectionBranches(S, planePoint, planeNormal, dirFlag, optsTry);
    if ok
        result.strategy = tries{t}.strategy;
        candidates(end+1,1) = result; %#ok<AGROW>
    end
end
end

function pick = localPickBest(candidates)
pick = candidates(1);
for k = 2:numel(candidates)
    a = pick;
    b = candidates(k);

    % Prefer candidates with fewer count changes, then more branches,
    % then more points, then lower fit RMS.
    aChanges = numel(a.countChangeIndices);
    bChanges = numel(b.countChangeIndices);

    if bChanges < aChanges
        pick = b;
    elseif bChanges == aChanges
        if b.numBranches > a.numBranches
            pick = b;
        elseif b.numBranches == a.numBranches
            if b.numPts > a.numPts
                pick = b;
            elseif b.numPts == a.numPts && b.fitRMS < a.fitRMS
                pick = b;
            end
        end
    end
end
end

function vals = localDenseValues(S, dirFlag, n)
if dirFlag == 'u'
    vals = linspace(S.domainV(1), S.domainV(2), n);
else
    vals = linspace(S.domainU(1), S.domainU(2), n);
end
end

function [result, ok] = localBuildSectionBranches(S, planePoint, planeNormal, dirFlag, opts)
result = struct();
ok = false;

fixedVals = localFixedValues(S, dirFlag, opts.FixedValues, opts.RefinePasses);
fixedVals = fixedVals(:);

samples = struct('fixed', {}, 'pts', {}, 'params', {});
for k = 1:numel(fixedVals)
    fk = fixedVals(k);

    if dirFlag == 'u'
        Ciso = S.isoCurveU(fk);   % fix v=fk, vary u
    else
        Ciso = S.isoCurveV(fk);   % fix u=fk, vary v
    end

    [ptsK, uK] = geom.intersectCurvePlaneBisection(Ciso, planePoint, planeNormal, ...
        'Tol', opts.Tol, 'MaxIter', opts.MaxIter);

    if isempty(uK)
        continue;
    end

    [uK, ord] = sort(uK(:));
    ptsK = ptsK(ord,:);

    samples(end+1).fixed = fk; %#ok<AGROW>
    samples(end).pts = ptsK;
    samples(end).params = uK;
end

if isempty(samples)
    return;
end

counts = arrayfun(@(s) size(s.pts,1), samples);
countChangeIdx = find(diff(counts) ~= 0);

if ~isempty(countChangeIdx)
    msg = sprintf(['Intersection count changed across sampled iso-curves at %d location(s). ' ...
        'Counts = [%s].'], numel(countChangeIdx), sprintf('%d ', counts));
    switch lower(string(opts.BranchCountAction))
        case "error"
            error('geom:surfacePlaneSection_multibranch:ChangingIntersectionCount', '%s', msg);
        case "warn"
            warning('geom:surfacePlaneSection_multibranch:ChangingIntersectionCount', '%s', msg);
        case "split"
            % no-op here; handled by segmenting below
        otherwise
            error('geom:surfacePlaneSection_multibranch', ...
                'Unknown BranchCountAction "%s".', char(opts.BranchCountAction));
    end
end

segments = localConstantCountSegments(counts);

branchesPts = {};
branchesIso = {};
branchesFix = {};

for s = 1:size(segments,1)
    i1 = segments(s,1);
    i2 = segments(s,2);
    c = counts(i1);

    if c <= 0
        continue;
    end

    for j = 1:c
        pts = zeros(i2-i1+1, 3);
        iso = zeros(i2-i1+1, 1);
        fixv = zeros(i2-i1+1, 1);
        row = 1;
        for k = i1:i2
            pts(row,:) = samples(k).pts(j,:);
            iso(row,1) = samples(k).params(j);
            fixv(row,1) = samples(k).fixed;
            row = row + 1;
        end

        % dedupe successive identical points
        if size(pts,1) >= 2
            keep = true(size(pts,1),1);
            for r = 2:size(pts,1)
                if norm(pts(r,:) - pts(r-1,:)) <= 1e-9
                    keep(r) = false;
                end
            end
            pts = pts(keep,:);
            iso = iso(keep);
            fixv = fixv(keep);
        end

        if size(pts,1) >= opts.MinBranchPoints
            branchesPts{end+1,1} = pts; %#ok<AGROW>
            branchesIso{end+1,1} = iso; %#ok<AGROW>
            branchesFix{end+1,1} = fixv; %#ok<AGROW>
        end
    end
end

if isempty(branchesPts)
    return;
end

curves = cell(numel(branchesPts),1);
fitRMS = inf(numel(branchesPts),1);
numPtsEach = zeros(numel(branchesPts),1);

for i = 1:numel(branchesPts)
    pts = branchesPts{i};

    deg = opts.Degree;
    if isempty(deg)
        deg = min(3, size(pts,1)-1);
    end
    deg = max(1, min(deg, size(pts,1)-1));

    Cfit = geom.NURBSCurve.globalInterp(pts, deg, opts.Method);
    uData = geom.NURBSCurve.parameterizeData(pts, opts.Method, deg);
    fitPts = Cfit.evaluate(uData);
    fitRMS(i) = sqrt(mean(sum((fitPts - pts).^2, 2)));

    curves{i} = Cfit;
    numPtsEach(i) = size(pts,1);
end

% sort by descending point count, then fit RMS
[~, order] = sortrows([-numPtsEach(:), fitRMS(:)], [1 2]);
curves = curves(order);
branchesPts = branchesPts(order);
branchesIso = branchesIso(order);
branchesFix = branchesFix(order);
fitRMS = fitRMS(order);
numPtsEach = numPtsEach(order);

result.curves = curves;
result.direction = dirFlag;
result.fixedValuesRequested = fixedVals;
result.fixedValuesUsed = unique(vertcat(branchesFix{:}));
result.branchPoints = branchesPts;
result.branchIsoParams = branchesIso;
result.branchFixedValues = branchesFix;
result.numPts = sum(numPtsEach);
result.fitRMS = min(fitRMS);
result.numBranches = numel(curves);
result.intersectionCounts = counts;
result.countChangeIndices = countChangeIdx;
ok = true;
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

function vals = localFixedValues(S, dirFlag, spec, refinePasses)
if isnumeric(spec)
    vals = spec(:).';
else
    switch lower(char(spec))
        case 'knot'
            if dirFlag == 'u'
                vals = unique(S.V);
            else
                vals = unique(S.U);
            end
        case 'uniform'
            if dirFlag == 'u'
                vals = linspace(S.domainV(1), S.domainV(2), max(8, numel(unique(S.V))));
            else
                vals = linspace(S.domainU(1), S.domainU(2), max(8, numel(unique(S.U))));
            end
        otherwise
            error('geom:surfacePlaneSection_multibranch', 'Unknown FixedValues option "%s".', char(spec));
    end
end

if dirFlag == 'u'
    a = S.domainV(1); b = S.domainV(2);
else
    a = S.domainU(1); b = S.domainU(2);
end
vals = vals(vals >= a - 1e-12 & vals <= b + 1e-12);
vals = unique(max(a, min(b, vals)));

for pass = 1:refinePasses
    if numel(vals) < 2, break; end
    mids = 0.5 * (vals(1:end-1) + vals(2:end));
    vals = sort(unique([vals, mids]));
end
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
            error('geom:surfacePlaneSection_multibranch', 'Primary plane must be ''x'', ''y'', or ''z''.');
    end
    if abs(offset) > 0
        planePoint = planePoint + offset * planeNormal;
    end
else
    planePoint = planeArg(:).';
    planeNormal = arg2(:).';
    if numel(planePoint) ~= 3 || numel(planeNormal) ~= 3
        error('geom:surfacePlaneSection_multibranch', ...
            'Arbitrary plane form requires point and normal as 1x3 vectors.');
    end
    nn = norm(planeNormal);
    if nn < eps
        error('geom:surfacePlaneSection_multibranch', 'Plane normal must be nonzero.');
    end
    planeNormal = planeNormal / nn;
    planePoint = planePoint + offset * planeNormal;
    desc = sprintf('point [%g %g %g], normal [%g %g %g], offset %g', ...
        planePoint(1), planePoint(2), planePoint(3), ...
        planeNormal(1), planeNormal(2), planeNormal(3), offset);
end
end
