function [C2, info] = curveRepositionTwoControls(C, ubar, deltaP, k, gamma)
% geom.curveRepositionTwoControls
% Exact implementation of The NURBS Book Sec. 11.2, Eqs. (11.7)-(11.9)

nodes = geom.curveControlNodes(C);
n = size(C.P,1);

if nargin < 4 || isempty(k)
    k = find(nodes(1:end-1) <= ubar & ubar <= nodes(2:end), 1, 'first');
    if isempty(k)
        [~, k] = min(abs(nodes(1:end-1) - ubar));
    end
end

if k < 1 || k >= n
    error('geom:curveRepositionTwoControls:index', ...
        'Need 1 <= k < number of control points.');
end

if nargin < 5 || isempty(gamma)
    tk  = nodes(k);
    tk1 = nodes(k+1);
    if abs(tk1 - tk) < 1e-14
        gamma = 0.5;
    else
        gamma = (ubar - tk) / (tk1 - tk);
    end
end

R = localRationalBasis(C, ubar);
Rk = R(k);
Rk1 = R(k+1);

w = (1-gamma)*Rk + gamma*Rk1;
if abs(w) < 1e-14
    error('geom:curveRepositionTwoControls:basis', ...
        'Weighted basis value is zero or too small.');
end

alphaV = deltaP ./ w;
dPk  = (1-gamma) * alphaV;
dPk1 = gamma      * alphaV;

P2 = C.P;
P2(k,:)   = P2(k,:)   + dPk;
P2(k+1,:) = P2(k+1,:) + dPk1;

C2 = geom.NURBSCurve(P2, C.p, C.U, C.W);

info = struct();
info.ubar = ubar;
info.nodes = nodes;
info.k = k;
info.gamma = gamma;
info.R = R;
info.Rk = Rk;
info.Rk1 = Rk1;
info.weightedBasis = w;
info.deltaP = deltaP;
info.dPk = dPk;
info.dPk1 = dPk1;
info.equation = 'Implements Sec. 11.2 Eqs. (11.7)-(11.9).';
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
    error('geom:curveRepositionTwoControls:den', ...
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
