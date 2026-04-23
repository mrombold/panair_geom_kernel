function [C2, info] = curveRepositionOneControl(C, ubar, deltaP, k)
% geom.curveRepositionOneControl
% Exact implementation of The NURBS Book Sec. 11.2, Eqs. (11.1)-(11.5)

if nargin < 4 || isempty(k)
    nodes = geom.curveControlNodes(C);
    [~, k] = min(abs(nodes - ubar));
else
    nodes = geom.curveControlNodes(C);
end

n = size(C.P,1);
if k < 1 || k > n
    error('geom:curveRepositionOneControl:index', 'Control index k out of range.');
end

R = localRationalBasis(C, ubar);
Rk = R(k);

if abs(Rk) < 1e-14
    error('geom:curveRepositionOneControl:basis', ...
        'R_{k,p}(ubar) is zero or too small; choose another k.');
end

P2 = C.P;
P2(k,:) = P2(k,:) + deltaP ./ Rk;

C2 = geom.NURBSCurve(P2, C.p, C.U, C.W);

info = struct();
info.ubar = ubar;
info.nodes = nodes;
info.k = k;
info.R = R;
info.Rk = Rk;
info.deltaP = deltaP;
info.controlTranslation = deltaP ./ Rk;
info.equation = 'Implements Sec. 11.2 Eqs. (11.1)-(11.5).';
end

function R = localRationalBasis(C, u)
p = C.p;
U = C.U(:).';
W = C.W(:).';
n = numel(W) - 1;

N = zeros(n+1,1);
for i0 = 0:n
    N(i0+1) = localBasisFun(i0, p, U, u);
end

den = sum(N .* W);
if abs(den) < 1e-14
    error('geom:curveRepositionOneControl:den', ...
        'Rational basis denominator is zero.');
end
R = (N .* W) / den;
end

function Nip = localBasisFun(i0, p, U, u)
if p == 0
    if (U(i0+1) <= u && u < U(i0+2)) || ...
       (u == U(end) && U(i0+1) <= u && u <= U(i0+2) && U(i0+2) == U(end))
        Nip = 1;
    else
        Nip = 0;
    end
    return;
end

leftDen  = U(i0+p+1) - U(i0+1);
rightDen = U(i0+p+2) - U(i0+2);

term1 = 0;
term2 = 0;
if leftDen ~= 0
    term1 = ((u - U(i0+1)) / leftDen) * localBasisFun(i0, p-1, U, u);
end
if rightDen ~= 0
    term2 = ((U(i0+p+2) - u) / rightDen) * localBasisFun(i0+1, p-1, U, u);
end
Nip = term1 + term2;
end
