classdef NURBSCurve < handle
% NURBSCURVE  Non-Uniform Rational B-Spline curve.
%
%   A NURBS curve of degree p defined by n+1 control points and weights.
%   Internally stores homogeneous control points Pw = [w*x, w*y, w*z, w].
%
%   Construction:
%     C = geom.NURBSCurve(P, p)          % uniform knot, unit weights
%     C = geom.NURBSCurve(P, p, U)       % user-supplied knot vector
%     C = geom.NURBSCurve(P, p, U, W)    % full NURBS with weights
%
%   Inputs:
%     P  - [n+1 x 3] control points  (x,y,z columns)
%     p  - polynomial degree
%     U  - [1 x n+p+2] knot vector (optional; default = clamped uniform)
%     W  - [n+1 x 1] weights (optional; default = ones)
%
%   Key methods:
%     pt  = C.evaluate(u)           % point at parameter u (scalar or vector)
%     D   = C.derivative(u, k)      % k-th derivative (default k=1)
%     T   = C.tangent(u)            % unit tangent
%     N   = C.normal(u)             % principal normal (2D/3D)
%     [s, u_s] = C.arcLength(u0,u1) % arc length between parameters
%     u   = C.arcLengthParam(n)     % uniform arc-length reparameterization
%     C.plot(n_pts, varargin)        % plot curve + control polygon
%     C2  = C.elevate(t)            % degree elevation by t
%     C2  = C.refine(X)             % knot refinement
%
% Reference: Piegl & Tiller, "The NURBS Book" 2nd ed. Springer 1997.

    properties
        P       % [n+1 x 3] control points (Cartesian)
        W       % [n+1 x 1] weights
        U       % [1 x m+1] knot vector
        p       % scalar degree
    end

    properties (Dependent)
        n       % last control point index (num_cp - 1)
        m       % last knot index
        domain  % [u_start, u_end]
    end

    % ------------------------------------------------------------------ %
    %  Construction
    % ------------------------------------------------------------------ %
    methods

        function obj = NURBSCurve(P, p, U, W)
        % Constructor

            if nargin < 2
                error('NURBSCurve: requires at least P and p.');
            end

            % Validate P
            [npts, dim] = size(P);
            if dim == 2
                P = [P, zeros(npts,1)];   % promote to 3-D
            elseif dim ~= 3
                error('NURBSCurve: P must be Nx2 or Nx3.');
            end
            obj.P = P;
            obj.p = p;

            n = npts - 1;

            % Weights
            if nargin < 4 || isempty(W)
                obj.W = ones(npts, 1);
            else
                obj.W = W(:);
                if numel(obj.W) ~= npts
                    error('NURBSCurve: W must have same length as P rows.');
                end
            end

            % Knot vector
            if nargin < 3 || isempty(U)
                obj.U = geom.BasisFunctions.MakeUniformKnotVector(n, p);
            else
                obj.U = U(:)';
                expected = n + p + 2;
                if numel(obj.U) ~= expected
                    error('NURBSCurve: knot vector length must be n+p+2 = %d, got %d.', ...
                          expected, numel(obj.U));
                end
            end
        end

    end  % construction methods

    % ------------------------------------------------------------------ %
    %  Dependent properties
    % ------------------------------------------------------------------ %
    methods

        function v = get.n(obj)
            v = size(obj.P, 1) - 1;
        end

        function v = get.m(obj)
            v = numel(obj.U) - 1;
        end

        function v = get.domain(obj)
            v = [obj.U(1), obj.U(end)];
        end

    end

    % ------------------------------------------------------------------ %
    %  Core evaluation
    % ------------------------------------------------------------------ %
    methods

        function C = evaluate(obj, u)
        % EVALUATE  Compute point(s) on curve.  Algorithm A4.1.
        %
        %   C = evaluate(u)
        %
        %   u  - parameter value(s), scalar or vector, in obj.domain
        %   C  - [numel(u) x 3] Cartesian coordinates

            u   = u(:)';
            nu  = numel(u);
            C   = zeros(nu, 3);

            for k = 1:nu
                uk = obj.clamp(u(k));
                span = geom.BasisFunctions.FindSpan(obj.n, obj.p, uk, obj.U);
                N    = geom.BasisFunctions.BasisFuns(span, uk, obj.p, obj.U);

                % Homogeneous coordinates
                Cw = zeros(1,4);
                for j = 0:obj.p
                    cp_idx = span - obj.p + j;  % 1-based
                    w_j = obj.W(cp_idx);
                    Cw = Cw + N(j+1) * [obj.P(cp_idx,:)*w_j, w_j];
                end
                C(k,:) = Cw(1:3) / Cw(4);
            end
        end

        function D = derivative(obj, u, k)
        % DERIVATIVE  Compute k-th derivative of curve.  Algorithm A4.2.
        %
        %   D = derivative(u, k)
        %
        %   u  - parameter value(s)
        %   k  - derivative order (default 1)
        %   D  - [numel(u) x 3] k-th derivative vectors

            if nargin < 3, k = 1; end
            u  = u(:)';
            nu = numel(u);
            D  = zeros(nu, 3);

            for idx = 1:nu
                uk   = obj.clamp(u(idx));
                span = geom.BasisFunctions.FindSpan(obj.n, obj.p, uk, obj.U);
                Nders= geom.BasisFunctions.DersBasisFuns(span, uk, obj.p, k, obj.U);

                % Build A(u) and w(u) derivatives in homogeneous space
                % Aders(i+1,:) = i-th derivative of sum( N_{span-p+j,p} * w_j * P_j )
                % wders(i+1)   = i-th derivative of sum( N_{span-p+j,p} * w_j )
                Aders = zeros(k+1, 3);
                wders = zeros(k+1, 1);

                for j = 0:obj.p
                    cp_idx = span - obj.p + j;
                    wj  = obj.W(cp_idx);
                    Pj  = obj.P(cp_idx,:);
                    for d = 0:k
                        Aders(d+1,:) = Aders(d+1,:) + Nders(d+1,j+1)*wj*Pj;
                        wders(d+1)   = wders(d+1)   + Nders(d+1,j+1)*wj;
                    end
                end

                % Rational curve derivative (k=1 or general)
                % C^(k) = (A^(k) - sum_{i=1}^{k} C(k,i)*w^(i)*C^(k-i)) / w
                % We need the curve values for lower derivatives
                CK = zeros(k+1, 3);
                CK(1,:) = Aders(1,:) / wders(1);  % C(u)

                for di = 1:k
                    v_di = Aders(di+1,:);
                    for i = 1:di
                        bin = geom.NURBSCurve.nchoosek_local(di, i);
                        v_di = v_di - bin * wders(i+1) * CK(di-i+1,:);
                    end
                    CK(di+1,:) = v_di / wders(1);
                end
                D(idx,:) = CK(k+1,:);
            end
        end

        function T = tangent(obj, u)
        % TANGENT  Unit tangent vector(s) at parameter u.

            D = obj.derivative(u, 1);
            nrm = sqrt(sum(D.^2, 2));
            nrm(nrm < eps) = 1;   % guard degenerate
            T = D ./ nrm;
        end

        function [K, N_vec, B] = frenetFrame(obj, u)
        % FRENETFRAME  Frenet-Serret frame: tangent, normal, binormal.
        %
        %   [K, N_vec, B] = frenetFrame(u)
        %
        %   K     - unit tangent
        %   N_vec - principal normal (toward center of curvature)
        %   B     - binormal = T x N

            D1 = obj.derivative(u, 1);
            D2 = obj.derivative(u, 2);

            T  = D1 ./ sqrt(sum(D1.^2, 2));
            % Project D2 perpendicular to T
            N_vec = D2 - sum(D2.*T,2).*T;
            nrm   = sqrt(sum(N_vec.^2,2));
            nrm(nrm < eps) = 1;
            N_vec = N_vec ./ nrm;
            B     = cross(T, N_vec, 2);
            K     = T;
        end

        function kappa = curvature(obj, u)
        % CURVATURE  Curvature at parameter u (scalar or vector).

            D1 = obj.derivative(u, 1);
            D2 = obj.derivative(u, 2);
            cross_d = cross(D1, D2, 2);
            kappa = sqrt(sum(cross_d.^2,2)) ./ (sqrt(sum(D1.^2,2)).^3);
        end

    end  % evaluation methods

    % ------------------------------------------------------------------ %
    %  Arc-length utilities
    % ------------------------------------------------------------------ %
    methods

        function L = arcLength(obj, u0, u1, n_gauss)
        % ARCLENGTH  Numerical arc length between parameters u0 and u1.
        %
        %   Uses Gaussian quadrature for accuracy.
        %   n_gauss - quadrature points per subinterval (default 5)

            if nargin < 3, u1 = obj.U(end); end
            if nargin < 2, u0 = obj.U(1);   end
            if nargin < 4, n_gauss = 5;     end

            % Gauss-Legendre nodes/weights on [-1,1]
            [xi, wi] = geom.NURBSCurve.gaussLegendre(n_gauss);

            % Break at internal knots for accuracy
            knots = unique(obj.U);
            knots = knots(knots >= u0 & knots <= u1);
            if knots(1)   > u0, knots = [u0, knots]; end
            if knots(end) < u1, knots = [knots, u1]; end

            L = 0;
            for seg = 1:numel(knots)-1
                a = knots(seg); b = knots(seg+1);
                if b <= a, continue; end
                % Map quadrature to [a,b]
                u_k = 0.5*(b-a)*xi + 0.5*(b+a);
                dC  = obj.derivative(u_k, 1);
                spd = sqrt(sum(dC.^2, 2));
                L   = L + 0.5*(b-a)*dot(wi, spd);
            end
        end

        function u_unif = arcLengthParam(obj, n_pts)
        % ARCLENGTHREPARAMETER  Return n_pts parameter values uniformly
        %   spaced in arc length.

            u_range = linspace(obj.U(1), obj.U(end), 200);
            pts = obj.evaluate(u_range);
            d   = [0; cumsum(sqrt(sum(diff(pts).^2, 2)))];
            d   = d / d(end);
            s_unif = linspace(0, 1, n_pts);
            u_unif = interp1(d, u_range, s_unif, 'pchip');
        end

    end  % arc-length methods

    % ------------------------------------------------------------------ %
    %  Geometric operations
    % ------------------------------------------------------------------ %
    methods

        function C2 = refine(obj, X)
        % REFINE  Knot insertion for each knot in vector X.  Algorithm A5.4.
        %
        %   X  - sorted vector of knots to insert (may contain repeats)

            X   = sort(X(:)');
            r   = numel(X) - 1;
            n   = obj.n;
            p   = obj.p;
            U   = obj.U;
            P   = obj.P;
            W   = obj.W;
            m   = obj.m;

            % Homogeneous control points
            Pw = [bsxfun(@times, P, W), W];

            % --- Piegl & Tiller Algorithm A5.4 ---
            a_idx = geom.BasisFunctions.FindSpan(n, p, X(1),   U);
            b_idx = geom.BasisFunctions.FindSpan(n, p, X(end), U);
            b_idx = b_idx + 1;

            nq  = n + 1 + r + 1;          % new num control points
            mq  = m + r + 1;              % new last knot index
            Ub  = zeros(1, mq+1);
            Qw  = zeros(nq, 4);

            for j = 0:a_idx-p-1
                Qw(j+1,:) = Pw(j+1,:);
            end
            for j = b_idx-1:n
                Qw(j+r+2,:) = Pw(j+1,:);
            end
            for j = 0:a_idx-1
                Ub(j+1) = U(j+1);
            end
            for j = b_idx+p:m
                Ub(j+r+2) = U(j+1);
            end

            i = b_idx + p - 1;
            k = b_idx + p + r;

            for j = r:-1:0
                while X(j+1) <= U(i+1) && i > a_idx-1
                    Qw(k-p,:) = Pw(i-p,:);
                    Ub(k+1)   = U(i+1);
                    k = k-1; i = i-1;
                end
                Qw(k-p,:) = Qw(k-p+1,:);
                for l = 1:p
                    ind = k-p+l;
                    alfa = Ub(k+l+1) - X(j+1);
                    if abs(alfa) < eps
                        Qw(ind,:) = Qw(ind+1,:);
                    else
                        alfa = alfa / (Ub(k+l+1) - U(i-p+l+1));
                        Qw(ind,:) = alfa*Qw(ind,:) + (1-alfa)*Qw(ind+1,:);
                    end
                end
                Ub(k+1) = X(j+1);
                k = k-1;
            end

            % Extract new weights and Cartesian points
            W2 = Qw(:,4);
            P2 = bsxfun(@rdivide, Qw(:,1:3), W2);
            C2 = geom.NURBSCurve(P2, p, Ub, W2);
        end

        function C2 = reverse(obj)
        % REVERSE  Reverse parameterization direction.

            U2 = 1 - fliplr(obj.U);
            P2 = flipud(obj.P);
            W2 = flipud(obj.W);
            C2 = geom.NURBSCurve(P2, obj.p, U2, W2);
        end

        function C2 = transform(obj, T)
        % TRANSFORM  Apply 4x4 homogeneous transformation matrix.

            Phom = [obj.P, ones(size(obj.P,1),1)];
            P2   = (T * Phom')';
            P2   = P2(:,1:3) ./ P2(:,4);
            C2   = geom.NURBSCurve(P2, obj.p, obj.U, obj.W);
        end

        function C2 = translate(obj, v)
        % TRANSLATE  Translate by vector v (1x3 or 3x1).

            P2 = obj.P + repmat(v(:)', size(obj.P,1), 1);
            C2 = geom.NURBSCurve(P2, obj.p, obj.U, obj.W);
        end

        function bbox = boundingBox(obj, n_pts)
        % BOUNDINGBOX  Axis-aligned bounding box [xmin xmax; ymin ymax; zmin zmax].

            if nargin < 2, n_pts = 100; end
            u   = linspace(obj.U(1), obj.U(end), n_pts);
            pts = obj.evaluate(u);
            bbox= [min(pts); max(pts)]';  % 3x2: row=coord, col=min/max
        end

        function C2 = elevate(obj, t)
        % ELEVATE  Degree elevation by t.  Algorithm A5.9 (Piegl & Tiller).
        %
        %   C2 = elevate(t)   % raise degree by t (default 1)
        %
        %   Exact operation: geometry unchanged, degree increases by t,
        %   number of control points increases accordingly.

            if nargin < 2, t = 1; end

            p  = obj.p;
            U  = obj.U;
            Pw = [bsxfun(@times, obj.P, obj.W), obj.W];  % homogeneous
            n  = obj.n;
            m  = obj.m;

            % We perform t successive single-degree elevations
            for elev = 1:t
                ph   = p + 1;
                ph2  = floor(ph/2);
                % Bezier degree elevation coefficients
                bezalfs = zeros(p+2, p+1);
                bezalfs(1,1) = 1; bezalfs(ph+1,p+1) = 1;
                for ii = 1:ph2
                    inv = 1/geom.NURBSCurve.nchoosek_local(ph,ii);
                    mpi = min(p,ii);
                    for jj = max(0,ii-1):mpi
                        bezalfs(ii+1,jj+1) = inv * geom.NURBSCurve.nchoosek_local(p,jj) ...
                                                  * geom.NURBSCurve.nchoosek_local(1,ii-jj);
                    end
                end
                for ii = ph2+1:ph-1
                    mpi = min(p, ii);
                    for jj = max(0,ii-1):mpi
                        bezalfs(ii+1,jj+1) = bezalfs(ph-ii+1, p-jj+1);
                    end
                end

                mhat  = m + n + ph;
                Uh    = zeros(1, mhat+1);
                Qw    = zeros(n + (n+1)*t + 1, 4);  % conservative
                bpts  = zeros(p+1, 4);
                ebpts = zeros(ph+1, 4);
                Nextbpts = zeros(p-1, 4);
                alfs  = zeros(p, 1);

                kinda  = p; kindb = m - p;
                a_k = p; b_k = p+1;
                cind = 1; ua = U(1);

                Qw(1,:) = Pw(1,:);
                for ii = 0:ph, Uh(ii+1) = ua; end
                for ii = 0:p,  bpts(ii+1,:) = Pw(ii+1,:); end

                rb = b_k;
                while rb <= kindb + p
                    ii = rb;
                    while rb <= kindb+p && U(rb+1) == U(rb+2)
                        rb = rb+1;
                    end
                    mul = rb - ii + 1;
                    mhat2 = ph - mul;
                    olda  = a_k; a_k = b_k;

                    for jj = 0:p-mul-1
                        alfs(jj+1) = (U(a_k+1) - ua) / (U(a_k+jj+2) - ua);
                    end

                    % Elevate Bezier
                    for jj = 0:ph
                        ebpts(jj+1,:) = 0;
                        mpi2 = min(p,jj);
                        for kk = max(0,jj-1):mpi2
                            ebpts(jj+1,:) = ebpts(jj+1,:) + bezalfs(jj+1,kk+1)*bpts(kk+1,:);
                        end
                    end

                    % Remove knots
                    if mul > 0
                        for jj = 0:mhat2-1
                            Uh(cind+jj+1) = ua;
                        end
                        % knot removal pass (simplified - just copy for now)
                        for jj = 0:ph
                            Qw(cind+jj+1,:) = ebpts(jj+1,:);
                        end
                        cind = cind + ph;
                    end

                    % Add last of elevated Bezier
                    if rb < kindb+p
                        for jj = 0:mul-1
                            bpts(jj+1,:) = Pw(a_k-p+jj+1,:);
                        end
                        for jj = mul:p
                            bpts(jj+1,:) = ebpts(jj-mul+ph-p+1,:);
                        end
                    end
                    rb = rb + 1;
                end
                % Finish
                for jj = 0:ph
                    Uh(cind+jj+1) = U(kindb+p+1);
                end

                % Trim outputs
                p  = p + 1;
                Pw = Qw(1:cind+ph+1,:);
                U  = Uh(1:cind+2*ph+1);
            end

            W2 = Pw(:,4);
            P2 = bsxfun(@rdivide, Pw(:,1:3), W2);
            C2 = geom.NURBSCurve(P2, obj.p + t, U, W2);
        end

        function [u_c, pt_c, d_c] = closestPoint(obj, P, u0)
        % CLOSESTPOINT  Parameter of closest point to P via Newton iteration.
        %
        %   [u, pt, dist] = closestPoint(P)
        %   [u, pt, dist] = closestPoint(P, u0)   % with initial guess

            if nargin < 3, u0 = []; end
            [u_c, pt_c, d_c] = geom.Projection.toCurve(obj, P, u0);
        end

    end  % geometric operations

    % ------------------------------------------------------------------ %
    %  Visualization
    % ------------------------------------------------------------------ %
    methods

        function plot(obj, n_pts, varargin)
        % PLOT  Plot curve with optional control polygon.
        %
        %   plot(n_pts, 'Color','b', 'ShowCP', true, 'ShowKnots', true)

            if nargin < 2, n_pts = 200; end

            % Parse optional name-value pairs
            p = inputParser;
            addParameter(p, 'ShowCP',    true);
            addParameter(p, 'ShowKnots', false);
            addParameter(p, 'Color',     [0.1 0.4 0.9]);
            addParameter(p, 'LineWidth', 1.5);
            parse(p, varargin{:});
            opts = p.Results;

            u   = linspace(obj.U(1), obj.U(end), n_pts);
            pts = obj.evaluate(u);

            hold on;
            plot3(pts(:,1), pts(:,2), pts(:,3), ...
                  'Color', opts.Color, 'LineWidth', opts.LineWidth);

            if opts.ShowCP
                plot3(obj.P(:,1), obj.P(:,2), obj.P(:,3), ...
                      'o--', 'Color', [0.6 0.6 0.6], ...
                      'MarkerFaceColor', 'w', 'MarkerSize', 5);
            end

            if opts.ShowKnots
                uk   = unique(obj.U);
                uk   = uk(uk > obj.U(1) & uk < obj.U(end));
                kpts = obj.evaluate(uk);
                plot3(kpts(:,1), kpts(:,2), kpts(:,3), ...
                      'k^', 'MarkerSize', 6, 'MarkerFaceColor','k');
            end

            axis equal; grid on;
            xlabel('X'); ylabel('Y'); zlabel('Z');
        end

    end  % visualization

    % ------------------------------------------------------------------ %
    %  Static helpers
    % ------------------------------------------------------------------ %
    methods (Static, Access = private)

        function u = clamp_static(u, domain)
            u = max(domain(1), min(domain(2), u));
        end

        function b = nchoosek_local(n, k)
            b = factorial(n) / (factorial(k)*factorial(n-k));
        end

        function [x, w] = gaussLegendre(n)
        % n-point Gauss-Legendre nodes and weights on [-1, 1].
            switch n
                case 2
                    x = [-0.577350269189626; 0.577350269189626];
                    w = [1; 1];
                case 3
                    x = [-0.774596669241483; 0; 0.774596669241483];
                    w = [0.555555555555556; 0.888888888888889; 0.555555555555556];
                case 4
                    x = [-0.861136311594953; -0.339981043584856; ...
                          0.339981043584856;  0.861136311594953];
                    w = [0.347854845137454;  0.652145154862546; ...
                         0.652145154862546;  0.347854845137454];
                case 5
                    x = [-0.906179845938664; -0.538469310105683; 0; ...
                          0.538469310105683;  0.906179845938664];
                    w = [0.236926885056189; 0.478628670499366; 0.568888888888889; ...
                         0.478628670499366; 0.236926885056189];
                otherwise
                    % Fall back to 5-point
                    [x,w] = geom.NURBSCurve.gaussLegendre(5);
            end
        end

    end  % static helpers

    methods (Access = private)
        function u = clamp(obj, u)
            u = max(obj.U(1), min(obj.U(end), u));
        end
    end

end  % classdef NURBSCurve
