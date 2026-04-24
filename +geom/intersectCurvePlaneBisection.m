function [pts, uvals, info] = intersectCurvePlaneBisection(C, planePoint, planeNormal, varargin)
% geom.intersectCurvePlaneBisection
% Find curve/plane intersections by sign-bracketing and bisection.

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
    error('geom:intersectCurvePlaneBisection', 'Plane normal must be nonzero.');
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
