classdef NURBSCurve < handle
% NURBSCURVE  Non-Uniform Rational B-Spline curve.
%
%   Rational B-spline curve using Cartesian control points P, weights W,
%   clamped or user-supplied knot vector U, and degree p.
%
%   This cleaned version keeps the exact rational evaluation/derivative path
%   and removes later duplicate methods that referenced inconsistent APIs.

    properties
        P       % [n+1 x 3] Cartesian control points
        W       % [n+1 x 1] weights
        U       % [1 x m+1] knot vector
        p       % degree
    end

    properties (Dependent)
        n       % last control point index
        m       % last knot index
        domain  % active parameter domain [U(p+1), U(m-p+1)]
    end

    methods
        function obj = NURBSCurve(P, p, U, W)
            if nargin < 2
                error('NURBSCurve: requires at least P and p.');
            end

            [npts, dim] = size(P);
            if dim == 2
                P = [P, zeros(npts,1)];
            elseif dim ~= 3
                error('NURBSCurve: P must be Nx2 or Nx3.');
            end
            if ~isscalar(p) || p < 0 || p ~= floor(p)
                error('NURBSCurve: p must be a nonnegative integer.');
            end
            n = npts - 1;
            if p > n
                error('NURBSCurve: degree p=%d exceeds n=%d.', p, n);
            end

            if nargin < 4 || isempty(W)
                W = ones(npts, 1);
            else
                W = W(:);
                if numel(W) ~= npts
                    error('NURBSCurve: W must have same length as P rows.');
                end
                if any(W <= 0)
                    error('NURBSCurve: weights must be strictly positive.');
                end
            end

            if nargin < 3 || isempty(U)
                U = geom.BasisFunctions.MakeUniformKnotVector(n, p);
            else
                U = U(:)';
            end
            expected = n + p + 2;
            if numel(U) ~= expected
                error('NURBSCurve: knot vector length must be n+p+2 = %d, got %d.', ...
                      expected, numel(U));
            end
            if any(diff(U) < 0)
                error('NURBSCurve: knot vector must be nondecreasing.');
            end

            obj.P = P;
            obj.W = W;
            obj.U = U;
            obj.p = p;
        end

        function v = get.n(obj)
            v = size(obj.P, 1) - 1;
        end

        function v = get.m(obj)
            v = numel(obj.U) - 1;
        end

        function v = get.domain(obj)
            % Active NURBS domain, not the full knot-array extent.
            v = [obj.U(obj.p+1), obj.U(end-obj.p)];
        end
    end

    methods
        function C = evaluate(obj, u)
            u = u(:)';
            nu = numel(u);
            C = zeros(nu, 3);

            for k = 1:nu
                uk = obj.clamp(u(k));
                span = geom.BasisFunctions.FindSpan(obj.n, obj.p, uk, obj.U);
                N = geom.BasisFunctions.BasisFuns(span, uk, obj.p, obj.U);

                Cw = zeros(1, 4);
                for j = 0:obj.p
                    cp_idx = span - obj.p + j;
                    wj = obj.W(cp_idx);
                    Cw = Cw + N(j+1) * [obj.P(cp_idx,:) * wj, wj];
                end
                C(k,:) = Cw(1:3) / Cw(4);
            end
        end

        function CK = derivatives(obj, u, d)
        % DERIVATIVES  Exact rational curve derivatives up to order d.
        % Returns CK{k+1} for the k-th derivative, or [d+1 x 3] for scalar u.
            if nargin < 3 || isempty(d)
                d = 1;
            end
            if numel(u) ~= 1
                error('NURBSCurve:derivatives', 'Use scalar u for derivatives().');
            end

            u = obj.clamp(u);
            d = min(d, obj.p);
            span = geom.BasisFunctions.FindSpan(obj.n, obj.p, u, obj.U);
            Nders = geom.BasisFunctions.DersBasisFuns(span, u, obj.p, d, obj.U);

            Aders = zeros(d+1, 3);
            wders = zeros(d+1, 1);

            for j = 0:obj.p
                cp_idx = span - obj.p + j;
                wj = obj.W(cp_idx);
                Pj = obj.P(cp_idx,:);
                for k = 0:d
                    Aders(k+1,:) = Aders(k+1,:) + Nders(k+1,j+1) * wj * Pj;
                    wders(k+1)   = wders(k+1)   + Nders(k+1,j+1) * wj;
                end
            end

            CK = zeros(d+1, 3);
            CK(1,:) = Aders(1,:) / wders(1);
            for k = 1:d
                v = Aders(k+1,:);
                for i = 1:k
                    v = v - nchoosek(k,i) * wders(i+1) * CK(k-i+1,:);
                end
                CK(k+1,:) = v / wders(1);
            end
        end

        function D = derivative(obj, u, k)
            if nargin < 3, k = 1; end
            u = u(:)';
            D = zeros(numel(u), 3);
            for idx = 1:numel(u)
                CK = obj.derivatives(u(idx), k);
                D(idx,:) = CK(k+1,:);
            end
        end

        function T = tangent(obj, u)
            D = obj.derivative(u, 1);
            nrm = sqrt(sum(D.^2, 2));
            nrm(nrm < eps) = 1;
            T = D ./ nrm;
        end

        function [T, Nvec, B] = frenetFrame(obj, u)
            D1 = obj.derivative(u, 1);
            D2 = obj.derivative(u, 2);

            T = D1 ./ sqrt(sum(D1.^2, 2));
            Nvec = D2 - sum(D2 .* T, 2) .* T;
            nrm = sqrt(sum(Nvec.^2, 2));
            nrm(nrm < eps) = 1;
            Nvec = Nvec ./ nrm;
            B = cross(T, Nvec, 2);
        end

        function kappa = curvature(obj, u)
            D1 = obj.derivative(u, 1);
            D2 = obj.derivative(u, 2);
            cross_d = cross(D1, D2, 2);
            kappa = sqrt(sum(cross_d.^2,2)) ./ max(sqrt(sum(D1.^2,2)).^3, eps);
        end
    end

    methods
        function L = arcLength(obj, u0, u1, n_gauss)
            if nargin < 4, n_gauss = 5; end
            if nargin < 3 || isempty(u1), u1 = obj.domain(2); end
            if nargin < 2 || isempty(u0), u0 = obj.domain(1); end
            u0 = obj.clamp(u0);
            u1 = obj.clamp(u1);
            if u1 < u0
                tmp = u0; u0 = u1; u1 = tmp;
            end

            [xi, wi] = geom.NURBSCurve.gaussLegendre(n_gauss);

            knots = unique(obj.U);
            knots = knots(knots >= u0 & knots <= u1);
            if isempty(knots) || knots(1) > u0, knots = [u0, knots]; end
            if knots(end) < u1, knots = [knots, u1]; end

            L = 0;
            for seg = 1:numel(knots)-1
                a = knots(seg); b = knots(seg+1);
                if b <= a, continue; end
                u_k = 0.5 * (b-a) * xi + 0.5 * (b+a);
                dC = obj.derivative(u_k, 1);
                spd = sqrt(sum(dC.^2, 2));
                L = L + 0.5 * (b-a) * dot(wi, spd);
            end
        end

        function u_unif = arcLengthParam(obj, n_pts)
            u_range = linspace(obj.domain(1), obj.domain(2), 200);
            pts = obj.evaluate(u_range);
            d = [0; cumsum(sqrt(sum(diff(pts).^2, 2)))];
            d = d / d(end);
            s_unif = linspace(0, 1, n_pts);
            u_unif = interp1(d, u_range, s_unif, 'pchip');
        end
    end

    methods
        function C2 = refine(obj, X)
            X = sort(X(:)');
            if isempty(X)
                C2 = geom.NURBSCurve(obj.P, obj.p, obj.U, obj.W);
                return;
            end

            r = numel(X) - 1;
            n = obj.n;
            p = obj.p;
            U = obj.U;
            P = obj.P;
            W = obj.W;
            m = obj.m;

            Pw = [P .* W, W];

            a_idx = geom.BasisFunctions.FindSpan(n, p, X(1), U);
            b_idx = geom.BasisFunctions.FindSpan(n, p, X(end), U) + 1;

            nq = n + r + 1;
            mq = m + r + 1;
            Ub = zeros(1, mq+1);
            Qw = zeros(nq+1, 4);

            for j = 0:a_idx-p
                Qw(j+1,:) = Pw(j+1,:);
            end
            for j = b_idx-1:n
                Qw(j+r+2,:) = Pw(j+1,:);
            end
            for j = 0:a_idx
                Ub(j+1) = U(j+1);
            end
            for j = b_idx+p:m
                Ub(j+r+2) = U(j+1);
            end

            i = b_idx + p - 1;
            k = b_idx + p + r;

            for j = r:-1:0
                while X(j+1) <= U(i+1) && i > a_idx
                    Qw(k-p,:) = Pw(i-p,:);
                    Ub(k+1) = U(i+1);
                    k = k - 1;
                    i = i - 1;
                end
                Qw(k-p,:) = Qw(k-p+1,:);
                for l = 1:p
                    ind = k - p + l;
                    alpha_num = Ub(k+l+1) - X(j+1);
                    alpha_den = Ub(k+l+1) - U(i-p+l+1);
                    if abs(alpha_num) < eps || abs(alpha_den) < eps
                        Qw(ind,:) = Qw(ind+1,:);
                    else
                        alpha = alpha_num / alpha_den;
                        Qw(ind,:) = alpha * Qw(ind,:) + (1 - alpha) * Qw(ind+1,:);
                    end
                end
                Ub(k+1) = X(j+1);
                k = k - 1;
            end

            W2 = Qw(:,4);
            P2 = Qw(:,1:3) ./ W2;
            C2 = geom.NURBSCurve(P2, p, Ub, W2);
        end

        function C2 = reverse(obj)
            U2 = (obj.domain(1) + obj.domain(2)) - fliplr(obj.U);
            P2 = flipud(obj.P);
            W2 = flipud(obj.W);
            C2 = geom.NURBSCurve(P2, obj.p, U2, W2);
        end

        function C2 = transform(obj, T)
            Ph = [obj.P, ones(size(obj.P,1),1)];
            Qh = (T * Ph')';
            P2 = Qh(:,1:3) ./ Qh(:,4);
            C2 = geom.NURBSCurve(P2, obj.p, obj.U, obj.W);
        end

        function C2 = translate(obj, v)
            P2 = obj.P + repmat(v(:)', size(obj.P,1), 1);
            C2 = geom.NURBSCurve(P2, obj.p, obj.U, obj.W);
        end

        function bbox = boundingBox(obj, n_pts)
            if nargin < 2, n_pts = 100; end
            u = linspace(obj.domain(1), obj.domain(2), n_pts);
            pts = obj.evaluate(u);
            bbox = [min(pts); max(pts)]';
        end

        function C2 = elevate(obj, t)
            %#ok<INUSD>
            error(['NURBSCurve:elevateNotImplementedExact', ...
                   ' Exact degree elevation has not been carried over yet. ', ...
                   'The previous implementation was incomplete and was removed ', ...
                   'to avoid a false-exact method in the cleaned kernel.']);
        end

        function [u_c, pt_c, d_c] = closestPoint(obj, P, u0)
            P = P(:)';
            dom = obj.domain;

            if nargin < 3 || isempty(u0)
                u0 = obj.coarseSearchPoint(P, 32);
            end
            u = geom.NURBSCurve.clamp_static(u0, dom);

            for iter = 1:50
                CK = obj.derivatives(u, 2);
                Ck = CK(1,:);
                D1 = CK(2,:);
                D2 = CK(3,:);

                r = Ck - P;
                f1 = dot(r, D1);
                f2 = dot(D1, D1) + dot(r, D2);

                if abs(f2) < eps
                    break;
                end

                du = -f1 / f2;
                u_new = geom.NURBSCurve.clamp_static(u + du, dom);
                if abs(u_new - u) < 1e-12
                    u = u_new;
                    break;
                end
                u = u_new;
            end

            pt_c = obj.evaluate(u);
            d_c = norm(pt_c - P);
            u_c = u;
        end

        function [u_all, pt_all, d_all] = closestPointBatch(obj, pts, n_coarse)
            if nargin < 3 || isempty(n_coarse)
                n_coarse = 64;
            end

            N = size(pts, 1);
            u_all = zeros(N, 1);
            pt_all = zeros(N, 3);
            d_all = zeros(N, 1);

            u_seed = linspace(obj.domain(1), obj.domain(2), n_coarse);
            pts_s = obj.evaluate(u_seed);

            for k = 1:N
                P = pts(k, :);
                d2 = sum((pts_s - P).^2, 2);
                [~, idx] = min(d2);
                [u_all(k), pt_all(k,:), d_all(k)] = obj.closestPoint(P, u_seed(idx));
            end
        end

        function [u_inv, pt_inv, err] = invertPoint(obj, P, u0)
            if nargin < 3, u0 = []; end
            [u_inv, pt_inv, err] = obj.closestPoint(P, u0);
        end
    end

    methods
        function plot(obj, n_pts, varargin)
            if nargin < 2, n_pts = 200; end

            p = inputParser;
            addParameter(p, 'ShowCP', true);
            addParameter(p, 'ShowKnots', false);
            addParameter(p, 'Color', [0.1 0.4 0.9]);
            addParameter(p, 'LineWidth', 1.5);
            parse(p, varargin{:});
            opts = p.Results;

            u = linspace(obj.domain(1), obj.domain(2), n_pts);
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
                uk = unique(obj.U);
                uk = uk(uk > obj.domain(1) & uk < obj.domain(2));
                if ~isempty(uk)
                    kpts = obj.evaluate(uk);
                    plot3(kpts(:,1), kpts(:,2), kpts(:,3), ...
                          'k^', 'MarkerSize', 6, 'MarkerFaceColor','k');
                end
            end

            axis equal; grid on;
            xlabel('X'); ylabel('Y'); zlabel('Z');
        end
    end

    methods (Static, Access = private)
        function u = clamp_static(u, domain)
            u = max(domain(1), min(domain(2), u));
        end

        function [x, w] = gaussLegendre(n)
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
                    [x, w] = geom.NURBSCurve.gaussLegendre(5);
            end
        end
    end

    methods (Access = private)
        function u = clamp(obj, u)
            u = max(obj.domain(1), min(obj.domain(2), u));
        end

        function u0 = coarseSearchPoint(obj, P, n)
            u_c = linspace(obj.domain(1), obj.domain(2), n);
            pts = obj.evaluate(u_c);
            d2 = sum((pts - P).^2, 2);
            [~, idx] = min(d2);
            u0 = u_c(idx);
        end
    end
end