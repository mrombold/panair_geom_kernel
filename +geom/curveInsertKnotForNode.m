function [C2, info] = curveInsertKnotForNode(C, ubar, k)
% geom.curveInsertKnotForNode
% Exact implementation of The NURBS Book Sec. 11.2 Eq. (11.6)

p = C.p;
U = C.U(:).';
nodes = geom.curveControlNodes(C);
n = size(C.P,1);

if nargin < 3 || isempty(k)
    k = find(nodes(1:end-1) <= ubar & ubar <= nodes(2:end), 1, 'first');
    if isempty(k)
        [~, k] = min(abs(nodes(1:end-1) - ubar));
    end
end

if k < 1 || k >= n
    error('geom:curveInsertKnotForNode:index', ...
        'Need 1 <= k < number of control points.');
end

if p < 2
    error('geom:curveInsertKnotForNode:degree', ...
        'Eq. (11.6) is only meaningful for degree p >= 2.');
end
uhat = p*ubar - sum(U((k+2):(k+p)));

[C2, span, mult] = localInsertOneKnot(C, uhat);

info = struct();
info.ubar = ubar;
info.k = k;
info.uhat = uhat;
info.insertSpan = span;
info.previousMultiplicity = mult;
info.qIndex = k + 1;
info.equation = 'Implements Sec. 11.2 Eq. (11.6).';
end

function [C2, kSpan, s] = localInsertOneKnot(C, u)
p = C.p;
U = C.U(:).';
P = C.P;
W = C.W(:);

Pw = [W .* P, W];

n = size(Pw,1) - 1;
kSpan = localFindSpan(n, p, u, U);
s = sum(abs(U - u) < 1e-12) - 1;
if s > p
    error('geom:curveInsertKnotForNode:multiplicity', ...
        'Knot multiplicity exceeds degree.');
end

Qw = zeros(n+2, size(Pw,2));
UQ = [U(1:kSpan+1), u, U(kSpan+2:end)];

Qw(1:(kSpan-p+1), :) = Pw(1:(kSpan-p+1), :);
Qw((kSpan-s+2):(n+2), :) = Pw((kSpan-s+1):(n+1), :);

for i0 = (kSpan-p+1):(kSpan-s)
    alpha = (u - U(i0+1)) / (U(i0+p+1) - U(i0+1));
    Qw(i0+1,:) = alpha * Pw(i0+1,:) + (1-alpha) * Pw(i0,:);
end

WQ = Qw(:,4);
PQ = Qw(:,1:3) ./ WQ;

C2 = geom.NURBSCurve(PQ, p, UQ, WQ);
end

function k = localFindSpan(n, p, u, U)
if u >= U(n+2)
    k = n;
    return;
end
if u <= U(p+1)
    k = p;
    return;
end

low = p;
high = n + 1;
mid = floor((low + high)/2);

while u < U(mid+1) || u >= U(mid+2)
    if u < U(mid+1)
        high = mid;
    else
        low = mid;
    end
    mid = floor((low + high)/2);
end
k = mid;
end
