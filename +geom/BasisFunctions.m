classdef BasisFunctions
% BASISFUNCTIONS  Core NURBS/B-spline basis function algorithms.
%
%   Static utility class implementing the fundamental algorithms from
%   Piegl & Tiller, "The NURBS Book" (2nd ed.), Springer 1997.
%   Algorithm numbers reference that text.
%
%   All methods are static - no instantiation needed.
%   Usage:  geom.BasisFunctions.FindSpan(n, p, u, U)
%
% Author:  Geometry Kernel v0.1

    methods (Static)

        function i = FindSpan(n, p, u, U)
        % FINDSPAN  Find the knot span index.  Algorithm A2.1.
        %
        %   i = FindSpan(n, p, u, U)
        %
        %   n  - last basis function index (num control points - 1)
        %   p  - degree
        %   u  - parameter value  [U(1), U(end)]
        %   U  - knot vector (row), length n+p+2
        %
        %   Returns i such that U(i+1) <= u < U(i+2)  (1-based index)

            % Special case: u at end of domain
            if u >= U(n+2)
                i = n + 1;  % 1-based: last internal span (n in 0-based -> n+1 in 1-based)
                return;
            end
            if u <= U(1)
                i = p + 1;  % first non-zero span (1-based)
                return;
            end

            low  = p + 1;     % 1-based
            high = n + 2;     % 1-based
            mid  = floor((low + high) / 2);

            while u < U(mid) || u >= U(mid+1)
                if u < U(mid)
                    high = mid;
                else
                    low  = mid;
                end
                mid = floor((low + high) / 2);
            end
            i = mid;  % 1-based span index
        end

        function N = BasisFuns(i, u, p, U)
        % BASISFUNS  Compute non-zero B-spline basis functions.  Algorithm A2.2.
        %
        %   N = BasisFuns(i, u, p, U)
        %
        %   i  - knot span index (1-based, from FindSpan)
        %   u  - parameter value
        %   p  - degree
        %   U  - knot vector
        %
        %   Returns N(1..p+1): the p+1 non-zero basis functions N_{i-p,p}..N_{i,p}

            N    = zeros(1, p+1);
            left = zeros(1, p+1);
            right= zeros(1, p+1);

            N(1) = 1.0;
            for j = 1:p
                left(j+1)  = u - U(i+1-j);
                right(j+1) = U(i+j) - u;
                saved = 0.0;
                for r = 0:j-1
                    temp     = N(r+1) / (right(r+2) + left(j-r+1));
                    N(r+1)   = saved + right(r+2) * temp;
                    saved    = left(j-r+1) * temp;
                end
                N(j+1) = saved;
            end
        end

        function ders = DersBasisFuns(i, u, p, n_ders, U)
        % DERSBASISFUNS  Compute basis functions and derivatives.  Algorithm A2.3.
        %
        %   ders = DersBasisFuns(i, u, p, n_ders, U)
        %
        %   i      - knot span index (1-based)
        %   u      - parameter value
        %   p      - degree
        %   n_ders - number of derivatives to compute (0 = function only)
        %   U      - knot vector
        %
        %   Returns ders(k+1, j+1): k-th derivative of N_{i-p+j, p}
        %   ders(1,:) = function values, ders(2,:) = first derivatives, etc.

            ndu  = zeros(p+1, p+1);
            a    = zeros(2,   p+1);
            ders = zeros(n_ders+1, p+1);
            left = zeros(p+1, 1);
            right= zeros(p+1, 1);

            ndu(1,1) = 1.0;
            for j = 1:p
                left(j+1)  = u - U(i+1-j);
                right(j+1) = U(i+j) - u;
                saved = 0.0;
                for r = 0:j-1
                    ndu(j+1, r+1) = right(r+2) + left(j-r+1);
                    temp = ndu(r+1,j) / ndu(j+1,r+1);
                    ndu(r+1, j+1) = saved + right(r+2)*temp;
                    saved = left(j-r+1)*temp;
                end
                ndu(j+1, j+1) = saved;
            end
            ders(1,:) = ndu(:,p+1)';

            for r = 0:p
                s1 = 0; s2 = 1;
                a(1,1) = 1.0;
                for k = 1:n_ders
                    d  = 0.0;
                    rk = r-k; pk = p-k;
                    if r >= k
                        a(s2+1,1) = a(s1+1,1) / ndu(pk+2, rk+1);
                        d = a(s2+1,1) * ndu(rk+1, pk+1);
                    end
                    if rk >= -1
                        j1 = 1;
                    else
                        j1 = -rk;
                    end
                    if r-1 <= pk
                        j2 = k-1;
                    else
                        j2 = p-r;
                    end
                    for j = j1:j2
                        a(s2+1, j+1) = (a(s1+1,j+1) - a(s1+1,j)) / ndu(pk+2, rk+j+1);
                        d = d + a(s2+1,j+1)*ndu(rk+j+1, pk+1);
                    end
                    if r <= pk
                        a(s2+1, k+1) = -a(s1+1,k) / ndu(pk+2, r+1);
                        d = d + a(s2+1,k+1)*ndu(r+1,pk+1);
                    end
                    ders(k+1,r+1) = d;
                    j = s1; s1 = s2; s2 = j;  % swap
                end
            end

            r = p;
            for k = 1:n_ders
                for j = 0:p
                    ders(k+1,j+1) = ders(k+1,j+1)*r;
                end
                r = r*(p-k);
            end
        end

        function U = MakeUniformKnotVector(n, p)
        % MAKEUNIFORMKNOTVECTOR  Clamped uniform knot vector.
        %
        %   U = MakeUniformKnotVector(n, p)
        %
        %   n  - last control point index (num_pts - 1)
        %   p  - degree
        %
        %   Returns clamped uniform knot vector of length n+p+2.

            m = n + p + 1;  % last knot index (0-based) => length m+1
            U = zeros(1, m+1);
            % First p+1 = 0, last p+1 = 1
            for j = 1:m+1
                if j <= p+1
                    U(j) = 0;
                elseif j >= m-p+1
                    U(j) = 1;
                else
                    U(j) = (j - p - 1) / (m - 2*p);
                end
            end
        end

        function U = ChordLengthKnotVector(pts, p)
        % CHORDLENGTHKNOTVECTOR  Chord-length parameterized clamped knot vector.
        %
        %   U = ChordLengthKnotVector(pts, p)
        %
        %   pts  - [n+1 x d] control (or data) points
        %   p    - degree
        %
        %   Returns clamped knot vector with interior knots placed at
        %   chord-length averaged parameter values.

            n   = size(pts,1) - 1;
            d   = size(pts,1);  % total points
            % Chord-length parameterization
            dk  = sqrt(sum(diff(pts).^2, 2));  % segment lengths
            tot = sum(dk);
            if tot < eps
                t = linspace(0,1,d)';
            else
                t = [0; cumsum(dk)/tot];
            end

            m   = n + p + 1;
            U   = zeros(1, m+1);
            % Interior knots (averaging)
            for j = 1:n-p
                U(p+1+j) = mean(t(j+1:j+p));
            end
            U(end-p:end) = 1;
        end

    end  % methods (Static)

end  % classdef
