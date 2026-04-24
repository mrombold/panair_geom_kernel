function [Csec, info] = surfacePlaneSection(S, planeArg, arg2, varargin)
% geom.surfacePlaneSection
% Construct an explicit NURBS section curve from a NURBS surface and a plane.
%
% Supported plane inputs:
%   [Csec,info] = geom.surfacePlaneSection(S, 'x', x0, ...)
%   [Csec,info] = geom.surfacePlaneSection(S, 'y', y0, ...)
%   [Csec,info] = geom.surfacePlaneSection(S, 'z', z0, ...)
%   [Csec,info] = geom.surfacePlaneSection(S, point, normal, ...)
%
% Name-value options
%   'Offset'        : extra offset along the supplied normal (default 0)
%   'Direction'     : 'auto' | 'u' | 'v' (default 'auto')
%   'FixedValues'   : 'knot' | 'uniform' | numeric vector (default 'knot')
%   'RefinePasses'  : integer midpoint refinement of fixed values (default 1)
%   'Tol'           : plane/curve bisection tolerance (default 1e-9)
%   'MaxIter'       : max bisection iterations (default 60)
%   'Degree'        : output curve degree (default min(3,nPts-1))
%   'Method'        : interpolation parameterization for output curve
%                     'uniform'|'chord'|'centripetal' (default 'centripetal')
%   'PreferNearest' : if multiple intersections on one iso-curve, follow the
%                     branch nearest the previous selected point (default true)

pa = inputParser;
addParameter(pa, 'Offset', 0.0);
addParameter(pa, 'Direction', 'auto');
addParameter(pa, 'FixedValues', 'knot');
addParameter(pa, 'RefinePasses', 1);
addParameter(pa, 'Tol', 1e-9);
addParameter(pa, 'MaxIter', 60);
addParameter(pa, 'Degree', []);
addParameter(pa, 'Method', 'centripetal');
addParameter(pa, 'PreferNearest', true);
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

    % robust fallback: if requested direction fails, also try the opposite
    if isempty(candidates)
        if char(dirMode) == 'u'
            candidates = localTryDirectionSet(S, planePoint, planeNormal, 'v', opts);
        else
            candidates = localTryDirectionSet(S, planePoint, planeNormal, 'u', opts);
        end
    end
end

if isempty(candidates)
    error('geom:surfacePlaneSection', ...
        'Could not build a usable section curve in requested direction "%s".', char(dirMode));
end

pick = localPickBest(candidates);

Csec = pick.C;
info = rmfield(pick, 'C');
info.planePoint = planePoint;
info.planeNormal = planeNormal;
info.planeDescription = planeDesc;
end

function candidates = localTryDirectionSet(S, planePoint, planeNormal, dirFlag, opts)
candidates = struct('C', {}, 'direction', {}, 'fixedValuesRequested', {}, ...
    'fixedValuesUsed', {}, 'points', {}, 'isoCurveParams', {}, ...
    'numPts', {}, 'fitRMS', {}, 'strategy', {});

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

    [result, ok] = localBuildSection(S, planePoint, planeNormal, dirFlag, optsTry);
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

    if b.numPts > a.numPts
        pick = b;
    elseif b.numPts == a.numPts
        if b.fitRMS < a.fitRMS
            pick = b;
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

function [result, ok] = localBuildSection(S, planePoint, planeNormal, dirFlag, opts)
result = struct();
ok = false;

fixedVals = localFixedValues(S, dirFlag, opts.FixedValues, opts.RefinePasses);
fixedVals = fixedVals(:);

pts = zeros(0,3);
curveParams = zeros(0,1);
fixedKept = zeros(0,1);

prevPoint = [];
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

    if numel(uK) == 1 || isempty(prevPoint) || ~opts.PreferNearest
        idx = 1;
    else
        d = vecnorm(ptsK - prevPoint, 2, 2);
        [~, idx] = min(d);
    end

    pts(end+1,:) = ptsK(idx,:); %#ok<AGROW>
    curveParams(end+1,1) = uK(idx); %#ok<AGROW>
    fixedKept(end+1,1) = fk; %#ok<AGROW>
    prevPoint = ptsK(idx,:);
end

if size(pts,1) >= 2
    keep = true(size(pts,1),1);
    for i = 2:size(pts,1)
        if norm(pts(i,:) - pts(i-1,:)) <= 1e-9
            keep(i) = false;
        end
    end
    pts = pts(keep,:);
    curveParams = curveParams(keep);
    fixedKept = fixedKept(keep);
end

if size(pts,1) < 2
    return;
end

deg = opts.Degree;
if isempty(deg)
    deg = min(3, size(pts,1)-1);
end
deg = max(1, min(deg, size(pts,1)-1));

Cfit = geom.NURBSCurve.globalInterp(pts, deg, opts.Method);

uData = geom.NURBSCurve.parameterizeData(pts, opts.Method, deg);
fitPts = Cfit.evaluate(uData);
fitRMS = sqrt(mean(sum((fitPts - pts).^2, 2)));

result.C = Cfit;
result.direction = dirFlag;
result.fixedValuesRequested = fixedVals;
result.fixedValuesUsed = fixedKept;
result.points = pts;
result.isoCurveParams = curveParams;
result.numPts = size(pts,1);
result.fitRMS = fitRMS;
ok = true;
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
            error('geom:surfacePlaneSection', 'Unknown FixedValues option "%s".', char(spec));
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
            error('geom:surfacePlaneSection', 'Primary plane must be ''x'', ''y'', or ''z''.');
    end
    if abs(offset) > 0
        planePoint = planePoint + offset * planeNormal;
    end
else
    planePoint = planeArg(:).';
    planeNormal = arg2(:).';
    if numel(planePoint) ~= 3 || numel(planeNormal) ~= 3
        error('geom:surfacePlaneSection', ...
            'Arbitrary plane form requires point and normal as 1x3 vectors.');
    end
    nn = norm(planeNormal);
    if nn < eps
        error('geom:surfacePlaneSection', 'Plane normal must be nonzero.');
    end
    planeNormal = planeNormal / nn;
    planePoint = planePoint + offset * planeNormal;
    desc = sprintf('point [%g %g %g], normal [%g %g %g], offset %g', ...
        planePoint(1), planePoint(2), planePoint(3), ...
        planeNormal(1), planeNormal(2), planeNormal(3), offset);
end
end
