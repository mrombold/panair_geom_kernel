function t = curveControlNodes(C)
% geom.curveControlNodes
% Compute the control-point nodes t_i from The NURBS Book Eq. (11.4):
%   t_i = (1/p) * sum_{j=1}^p u_{i+j},   i = 0,...,n

p = C.p;
n = size(C.P,1) - 1;
U = C.U(:).';

t = zeros(n+1,1);
for i0 = 0:n
    t(i0+1) = mean(U((i0+1)+1 : (i0+1)+p));
end
end
