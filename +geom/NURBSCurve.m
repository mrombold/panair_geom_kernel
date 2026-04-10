classdef NURBSCurve < handle
% NURBSCURVE  Rational B-spline curve with exact evaluation/derivative core.
%
%   This version is cleaned up around the standard Piegl & Tiller curve
%   workflow: evaluation, exact rational derivatives, knot insertion,
%   Bezier decomposition, exact degree elevation, splitting, interpolation,
%   and projection.

    properties
        P   % [n+1 x 3] Cartesian control points
        W   % [n+1 x 1] positive weights
        U   % [1 x m+1] nondecreasing knot vector
        p   % degree
    end

    properties (Dependent)
        n       % last control-point index
        m       % last knot index
        domain  % active parameter domain [U(p+1), U(end-p)]
        Pw      % homogeneous control points [wx wy wz w]
    end

    methods
        function obj = NURBSCurve(P, p, U, W)
            if nargin < 2
                error('NURBSCurve:Constructor', 'Requires at least P and p.');
            end

            [npts, dim] = size(P);
            if dim == 2
                P = [P, zeros(npts,1)];
            elseif dim ~= 3
                error('NURBSCurve:Constructor', 'P must be Nx2 or Nx3.');
            end
            if ~isscalar(p) || p < 0 || p ~= floor(p)
                error('NURBSCurve:Constructor', 'Degree p must be a nonnegative integer.');
            end
            n = npts - 1;
            if p > n
                error('NURBSCurve:Constructor', 'Degree p=%d exceeds n=%d.', p, n);
            end

            if nargin < 4 || isempty(W)
                W = ones(npts,1);
            else
                W = W(:);
                if numel(W) ~= npts
                    error('NURBSCurve:Constructor', 'W must match the number of control points.');
                end
                if any(W <= 0)
                    error('NURBSCurve:Constructor', 'Weights must be strictly positive.');
                end
            end

            if nargin < 3 || isempty(U)
                U = geom.BasisFunctions.MakeUniformKnotVector(n, p);
            else
                U = U(:)';
            end
            if numel(U) ~= n + p + 2
                error('NURBSCurve:Constructor', 'Knot-vector length must be n+p+2.');
            end
            if any(diff(U) < 0)
                error('NURBSCurve:Constructor', 'Knot vector must be nondecreasing.');
            end

            obj.P = P;
            obj.W = W;
            obj.U = U;
            obj.p = p;
            obj.validate();
        end

        function v = get.n(obj)
            v = size(obj.P,1) - 1;
        end

        function v = get.m(obj)
            v = numel(obj.U) - 1;
        end

        function v = get.domain(obj)
            v = [obj.U(obj.p+1), obj.U(end-obj.p)];
        end

        function v = get.Pw(obj)
            v = [obj.P .* obj.W, obj.W];
        end
    end

    methods
        function validate(obj)
            if size(obj.P,1) ~= numel(obj.W)
                error('NURBSCurve:Validate', 'P and W sizes are inconsistent.');
            end
            if numel(obj.U) ~= obj.n + obj.p + 2
                error('NURBSCurve:Validate', 'Knot-vector length is inconsistent with n and p.');
            end
            if any(diff(obj.U) < 0)
                error('NURBSCurve:Validate', 'Knot vector must be nondecreasing.');
            end
            if any(obj.W <= 0)
                error('NURBSCurve:Validate', 'Weights must be strictly positive.');
            end
        end

        function tf = isClamped(obj)
            tf = all(abs(obj.U(1:obj.p+1) - obj.U(1)) < 1e-12) && ...
                 all(abs(obj.U(end-obj.p:end) - obj.U(end)) < 1e-12);
        end

        function s = knotMultiplicity(obj, u, tol)
            if nargin < 3 || isempty(tol), tol = 1e-12; end
            s = sum(abs(obj.U - u) < tol);
        end

        function C = evaluate(obj, u)
            u = u(:)';
            C = zeros(numel(u), 3);
            for k = 1:numel(u)
                uk = obj.clamp(u(k));
                span = geom.BasisFunctions.FindSpan(obj.n, obj.p, uk, obj.U);
                N = geom.BasisFunctions.BasisFuns(span, uk, obj.p, obj.U);
                Cw = zeros(1,4);
                for j = 0:obj.p
                    idx = span - obj.p + j;
                    Cw = Cw + N(j+1) * obj.Pw(idx,:);
                end
                C(k,:) = Cw(1:3) / Cw(4);
            end
        end

        function CK = derivatives(obj, u, d)
        % Exact rational derivatives C^(0)...C^(d) at scalar u.
            if nargin < 3 || isempty(d), d = 1; end
            if ~isscalar(u)
                error('NURBSCurve:derivatives', 'u must be scalar.');
            end
            u = obj.clamp(u);
            d = min(d, obj.p);

            span = geom.BasisFunctions.FindSpan(obj.n, obj.p, u, obj.U);
            Nders = geom.BasisFunctions.DersBasisFuns(span, u, obj.p, d, obj.U);

            Aders = zeros(d+1, 3);
            wders = zeros(d+1, 1);
            for j = 0:obj.p
                idx = span - obj.p + j;
                wj = obj.W(idx);
                Pj = obj.P(idx,:);
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
            for i = 1:numel(u)
                CK = obj.derivatives(u(i), k);
                D(i,:) = CK(k+1,:);
            end
        end

        function T = tangent(obj, u)
            D = obj.derivative(u, 1);
            nrm = sqrt(sum(D.^2,2));
            nrm(nrm < eps) = 1;
            T = D ./ nrm;
        end

        function [T, Nvec, B] = frenetFrame(obj, u)
            D1 = obj.derivative(u, 1);
            D2 = obj.derivative(u, 2);
            T = D1 ./ max(vecnorm(D1,2,2), eps);
            Nvec = D2 - sum(D2 .* T, 2) .* T;
            Nvec = Nvec ./ max(vecnorm(Nvec,2,2), eps);
            B = cross(T, Nvec, 2);
        end

        function kappa = curvature(obj, u)
            D1 = obj.derivative(u,1);
            D2 = obj.derivative(u,2);
            num = vecnorm(cross(D1,D2,2),2,2);
            den = max(vecnorm(D1,2,2).^3, eps);
            kappa = num ./ den;
        end

        function L = arcLength(obj, u0, u1, nGauss)
            if nargin < 4 || isempty(nGauss), nGauss = 5; end
            if nargin < 3 || isempty(u1), u1 = obj.domain(2); end
            if nargin < 2 || isempty(u0), u0 = obj.domain(1); end
            u0 = obj.clamp(u0); u1 = obj.clamp(u1);
            if u1 < u0, tmp = u0; u0 = u1; u1 = tmp; end

            [xi, wi] = geom.NURBSCurve.gaussLegendre(nGauss);
            K = unique(obj.U);
            K = K(K >= u0 & K <= u1);
            if isempty(K) || K(1) > u0, K = [u0, K]; end
            if K(end) < u1, K = [K, u1]; end

            L = 0;
            for s = 1:numel(K)-1
                a = K(s); b = K(s+1);
                if b <= a, continue; end
                uq = 0.5*(b-a)*xi + 0.5*(a+b);
                dC = obj.derivative(uq, 1);
                L = L + 0.5*(b-a) * dot(wi, vecnorm(dC,2,2));
            end
        end

        function u = arcLengthParam(obj, nPts)
            us = linspace(obj.domain(1), obj.domain(2), max(200, 10*nPts));
            pts = obj.evaluate(us);
            s = [0; cumsum(vecnorm(diff(pts,1,1),2,2))];
            s = s / s(end);
            u = interp1(s, us(:), linspace(0,1,nPts)', 'pchip')';
        end
    end

    methods
        function C2 = insertKnot(obj, u, r)
        % Exact knot insertion by repeated refinement.
            if nargin < 3 || isempty(r), r = 1; end
            if r < 0 || r ~= floor(r)
                error('NURBSCurve:insertKnot', 'r must be a nonnegative integer.');
            end
            if r == 0
                C2 = geom.NURBSCurve(obj.P, obj.p, obj.U, obj.W);
                return;
            end
            X = repmat(obj.clamp(u), 1, r);
            C2 = obj.refine(X);
        end

        function C2 = refine(obj, X)
        % Exact curve knot-refinement (A5.4 style vector insertion).
            X = sort(X(:)');
            if isempty(X)
                C2 = geom.NURBSCurve(obj.P, obj.p, obj.U, obj.W);
                return;
            end

            r = numel(X) - 1;
            n = obj.n;
            p = obj.p;
            U = obj.U;
            Pw = obj.Pw;
            m = obj.m;

            a = geom.BasisFunctions.FindSpan(n, p, X(1), U);
            b = geom.BasisFunctions.FindSpan(n, p, X(end), U) + 1;

            nq = n + r + 1;
            mq = m + r + 1;
            Ub = zeros(1, mq+1);
            Qw = zeros(nq+1, 4);

            for j = 0:(a-p)
                Qw(j+1,:) = Pw(j+1,:);
            end
            for j = (b-1):n
                Qw(j+r+2,:) = Pw(j+1,:);
            end
            for j = 0:a
                Ub(j+1) = U(j+1);
            end
            for j = (b+p):m
                Ub(j+r+2) = U(j+1);
            end

            i = b + p - 1;
            k = b + p + r;
            for j = r:-1:0
                while X(j+1) <= U(i+1) && i > a
                    Qw(k-p+1,:) = Pw(i-p+1,:);
                    Ub(k+1) = U(i+1);
                    k = k - 1;
                    i = i - 1;
                end
                Qw(k-p+1,:) = Qw(k-p+2,:);
                for l = 1:p
                    ind = k - p + l;
                    alpha = Ub(k+l+1) - X(j+1);
                    if abs(alpha) < 1e-14
                        Qw(ind+1,:) = Qw(ind+2,:);
                    else
                        alpha = alpha / (Ub(k+l+1) - U(i-p+l+1));
                        Qw(ind+1,:) = alpha * Qw(ind+1,:) + (1-alpha) * Qw(ind+2,:);
                    end
                end
                Ub(k+1) = X(j+1);
                k = k - 1;
            end

            W2 = Qw(:,4);
            P2 = Qw(:,1:3) ./ W2;
            C2 = geom.NURBSCurve(P2, p, Ub, W2);
        end

        function [Cleft, Cright] = split(obj, u0)
        % Exact split by knot insertion to full multiplicity.
            u0 = obj.clamp(u0);
            tol = 1e-12;
            if u0 <= obj.domain(1)+tol || u0 >= obj.domain(2)-tol
                error('NURBSCurve:split', 'Split parameter must be strictly interior to the active domain.');
            end
            s = obj.knotMultiplicity(u0, tol);
            if s < obj.p
                Cref = obj.insertKnot(u0, obj.p - s);
            else
                Cref = obj;
            end
            U2 = Cref.U;
            first = find(abs(U2 - u0) < tol, 1, 'first');
            last  = find(abs(U2 - u0) < tol, 1, 'last');
            idxSplit = last - obj.p;

            P1 = Cref.P(1:idxSplit,:);
            W1 = Cref.W(1:idxSplit);
            U1 = [U2(1:last), u0];
            U1 = geom.NURBSCurve.normalizeKnotVector(U1);

            P2 = Cref.P(idxSplit:end,:);
            W2 = Cref.W(idxSplit:end);
            U3 = [u0, U2(first:end)];
            U3 = geom.NURBSCurve.normalizeKnotVector(U3);

            Cleft  = geom.NURBSCurve(P1, obj.p, U1, W1);
            Cright = geom.NURBSCurve(P2, obj.p, U3, W2);
        end

        function parts = decomposeBezier(obj)
        % Exact Bezier decomposition by inserting each internal knot to multiplicity p.
            Ui = unique(obj.U(obj.p+2:end-obj.p-1));
            Cw = obj;
            for i = 1:numel(Ui)
                s = Cw.knotMultiplicity(Ui(i));
                if s < Cw.p
                    Cw = Cw.insertKnot(Ui(i), Cw.p - s);
                end
            end

            U = Cw.U;
            breaks = unique(U);
            parts = cell(0,1);
            seg = 0;
            for k = 1:numel(breaks)-1
                if breaks(k+1) <= breaks(k), continue; end
                seg = seg + 1;
                i0 = (seg-1)*Cw.p + 1;
                i1 = i0 + Cw.p;
                Pseg = Cw.P(i0:i1,:);
                Wseg = Cw.W(i0:i1);
                Useg = [zeros(1,Cw.p+1), ones(1,Cw.p+1)];
                Cseg = geom.NURBSCurve(Pseg, Cw.p, Useg, Wseg);
                parts{seg,1} = struct('curve', Cseg, 'u0', breaks(k), 'u1', breaks(k+1));
            end
        end

        function C2 = elevate(obj, t)
        % Exact degree elevation via Bezier decomposition and reassembly.
            if nargin < 2 || isempty(t), t = 1; end
            if t < 0 || t ~= floor(t)
                error('NURBSCurve:elevate', 't must be a nonnegative integer.');
            end
            if t == 0
                C2 = geom.NURBSCurve(obj.P, obj.p, obj.U, obj.W);
                return;
            end

            parts = obj.decomposeBezier();
            ph = obj.p + t;
            nSeg = numel(parts);
            if nSeg == 0
                error('NURBSCurve:elevate', 'Bezier decomposition failed.');
            end

            % Elevate each Bezier segment in homogeneous space.
            PwAll = [];
            breaks = zeros(nSeg+1,1);
            for s = 1:nSeg
                Cseg = parts{s}.curve;
                PwSeg = Cseg.Pw;
                Pwh = geom.NURBSCurve.elevateBezierHomogeneous(PwSeg, obj.p, t);
                if s == 1
                    PwAll = Pwh;
                    breaks(1) = parts{s}.u0;
                else
                    PwAll = [PwAll; Pwh(2:end,:)]; %#ok<AGROW>
                end
                breaks(s+1) = parts{s}.u1;
            end

            Unew = [repmat(breaks(1), 1, ph+1)];
            for s = 2:nSeg
                Unew = [Unew, repmat(breaks(s), 1, ph)]; %#ok<AGROW>
            end
            Unew = [Unew, repmat(breaks(end), 1, ph+1)];

            W2 = PwAll(:,4);
            P2 = PwAll(:,1:3) ./ W2;
            C2 = geom.NURBSCurve(P2, ph, Unew, W2);
        end

        function C2 = reverse(obj)
            U2 = (obj.U(1) + obj.U(end)) - fliplr(obj.U);
            P2 = flipud(obj.P);
            W2 = flipud(obj.W);
            C2 = geom.NURBSCurve(P2, obj.p, U2, W2);
        end

        function C2 = transform(obj, T)
            Ph = [obj.P, ones(size(obj.P,1),1)];
            Qt = (T * Ph')';
            Q = Qt(:,1:3) ./ Qt(:,4);
            C2 = geom.NURBSCurve(Q, obj.p, obj.U, obj.W);
        end

        function C2 = translate(obj, v)
            Q = obj.P + repmat(v(:)', size(obj.P,1), 1);
            C2 = geom.NURBSCurve(Q, obj.p, obj.U, obj.W);
        end

        function bbox = boundingBox(obj, nPts)
            if nargin < 2 || isempty(nPts), nPts = 200; end
            pts = obj.evaluate(linspace(obj.domain(1), obj.domain(2), nPts));
            bbox = [min(pts); max(pts)]';
        end
    end

    methods
        function [u_c, pt_c, d_c] = closestPoint(obj, P, u0)
            P = P(:)';
            if nargin < 3 || isempty(u0)
                u0 = obj.coarseSearchPoint(P, 64);
            end
            u = obj.clamp(u0);
            for iter = 1:50
                CK = obj.derivatives(u, 2);
                C  = CK(1,:);
                D1 = CK(2,:);
                D2 = CK(3,:);
                r  = C - P;
                f1 = dot(r, D1);
                f2 = dot(D1, D1) + dot(r, D2);
                if abs(f2) < 1e-14, break; end
                un = obj.clamp(u - f1/f2);
                if abs(un - u) < 1e-12
                    u = un;
                    break;
                end
                u = un;
            end
            pt_c = obj.evaluate(u);
            d_c  = norm(pt_c - P);
            u_c  = u;
        end

        function [u_all, pt_all, d_all] = closestPointBatch(obj, pts, nCoarse)
            if nargin < 3 || isempty(nCoarse), nCoarse = 64; end
            N = size(pts,1);
            u_all = zeros(N,1);
            pt_all = zeros(N,3);
            d_all = zeros(N,1);
            seeds = linspace(obj.domain(1), obj.domain(2), nCoarse);
            samp = obj.evaluate(seeds);
            for i = 1:N
                d2 = sum((samp - pts(i,:)).^2, 2);
                [~, idx] = min(d2);
                [u_all(i), pt_all(i,:), d_all(i)] = obj.closestPoint(pts(i,:), seeds(idx));
            end
        end

        function [u_inv, pt_inv, err] = invertPoint(obj, P, u0)
            if nargin < 3, u0 = []; end
            [u_inv, pt_inv, err] = obj.closestPoint(P, u0);
        end
    end

    methods
        function plot(obj, nPts, varargin)
            if nargin < 2 || isempty(nPts), nPts = 200; end
            pa = inputParser;
            addParameter(pa, 'ShowCP', true);
            addParameter(pa, 'ShowKnots', false);
            addParameter(pa, 'Color', [0 0.4470 0.7410]);
            addParameter(pa, 'LineWidth', 1.5);
            parse(pa, varargin{:});
            opt = pa.Results;

            u = linspace(obj.domain(1), obj.domain(2), nPts);
            pts = obj.evaluate(u);
            hold on;
            plot3(pts(:,1), pts(:,2), pts(:,3), 'Color', opt.Color, 'LineWidth', opt.LineWidth);
            if opt.ShowCP
                plot3(obj.P(:,1), obj.P(:,2), obj.P(:,3), 'o--', 'Color', [0.6 0.6 0.6], ...
                    'MarkerFaceColor', 'w', 'MarkerSize', 5);
            end
            if opt.ShowKnots
                uk = unique(obj.U);
                uk = uk(uk > obj.domain(1) & uk < obj.domain(2));
                if ~isempty(uk)
                    kp = obj.evaluate(uk);
                    plot3(kp(:,1), kp(:,2), kp(:,3), 'k^', 'MarkerFaceColor', 'k');
                end
            end
            axis equal; grid on; xlabel('X'); ylabel('Y'); zlabel('Z');
        end
    end

    methods (Static)
        function U = averagingKnotVector(params, p)
            params = params(:)';
            m = numel(params) - 1;
            n = m;
            U = zeros(1, n+p+2);
            U(1:p+1) = params(1);
            U(end-p:end) = params(end);
            for j = 1:(n-p)
                U(j+p+1) = sum(params(j+1:j+p)) / p;
            end
        end

        function u = parameterizeData(Q, method)
            if nargin < 2 || isempty(method), method = 'centripetal'; end
            Q = geom.NURBSCurve.ensure3D(Q);
            m = size(Q,1) - 1;
            d = zeros(m,1);
            for k = 1:m
                dk = norm(Q(k+1,:) - Q(k,:));
                switch lower(method)
                    case 'uniform'
                        d(k) = 1;
                    case 'chord'
                        d(k) = dk;
                    otherwise
                        d(k) = sqrt(dk);
                end
            end
            total = sum(d);
            u = zeros(m+1,1);
            if total < eps
                u = linspace(0,1,m+1)';
                return;
            end
            for k = 2:m+1
                u(k) = u(k-1) + d(k-1)/total;
            end
            u(end) = 1;
        end

        function C = globalInterp(Q, p, method)
        % Global interpolation of data points with unit weights.
            if nargin < 3 || isempty(method), method = 'centripetal'; end
            Q = geom.NURBSCurve.ensure3D(Q);
            n = size(Q,1) - 1;
            if p > n
                error('NURBSCurve:globalInterp', 'Degree exceeds number of data points minus one.');
            end
            u = geom.NURBSCurve.parameterizeData(Q, method);
            U = geom.NURBSCurve.averagingKnotVector(u, p);
            A = zeros(n+1, n+1);
            for k = 1:n+1
                span = geom.BasisFunctions.FindSpan(n, p, u(k), U);
                N = geom.BasisFunctions.BasisFuns(span, u(k), p, U);
                A(k, span-p+1:span+1) = N;
            end
            P = A \ Q;
            C = geom.NURBSCurve(P, p, U, ones(n+1,1));
        end

        function C = globalLeastSquaresFit(Q, p, nCtrl, method)
        % Global least-squares fit of data Q with nCtrl control points.
            if nargin < 4 || isempty(method), method = 'centripetal'; end
            Q = geom.NURBSCurve.ensure3D(Q);
            m = size(Q,1) - 1;
            n = nCtrl - 1;
            if n < p
                error('NURBSCurve:globalLeastSquaresFit', 'Need nCtrl >= p+1.');
            end
            u = geom.NURBSCurve.parameterizeData(Q, method);
            U = geom.NURBSCurve.averagingKnotVector(linspace(0,1,n+1)', p);
            Nmat = zeros(m+1, n+1);
            for k = 1:m+1
                span = geom.BasisFunctions.FindSpan(n, p, u(k), U);
                N = geom.BasisFunctions.BasisFuns(span, u(k), p, U);
                Nmat(k, span-p+1:span+1) = N;
            end
            P = Nmat \ Q;
            C = geom.NURBSCurve(P, p, U, ones(n+1,1));
        end
    end

    methods (Static, Access = private)
        function PwElev = elevateBezierHomogeneous(Pw, p, t)
            ph = p + t;
            PwElev = zeros(ph+1, size(Pw,2));
            for i = 0:ph
                j0 = max(0, i-t);
                j1 = min(p, i);
                for j = j0:j1
                    coeff = nchoosek(p,j) * nchoosek(t, i-j) / nchoosek(ph,i);
                    PwElev(i+1,:) = PwElev(i+1,:) + coeff * Pw(j+1,:);
                end
            end
        end

        function U = normalizeKnotVector(U)
            a = U(1); b = U(end);
            if abs(b-a) < eps
                U = zeros(size(U));
            else
                U = (U - a) / (b - a);
            end
        end

        function P = ensure3D(P)
            if size(P,2) == 2
                P = [P, zeros(size(P,1),1)];
            end
        end

        function u = clamp_static(u, dom)
            u = max(dom(1), min(dom(2), u));
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
                    x = [-0.861136311594953; -0.339981043584856; 0.339981043584856; 0.861136311594953];
                    w = [0.347854845137454; 0.652145154862546; 0.652145154862546; 0.347854845137454];
                otherwise
                    x = [-0.906179845938664; -0.538469310105683; 0; 0.538469310105683; 0.906179845938664];
                    w = [0.236926885056189; 0.478628670499366; 0.568888888888889; 0.478628670499366; 0.236926885056189];
            end
        end
    end

    methods (Access = private)
        function u = clamp(obj, u)
            u = max(obj.domain(1), min(obj.domain(2), u));
        end

        function u0 = coarseSearchPoint(obj, P, n)
            us = linspace(obj.domain(1), obj.domain(2), n);
            pts = obj.evaluate(us);
            [~, idx] = min(sum((pts - P).^2,2));
            u0 = us(idx);
        end
    end


    methods
        function c = continuityAt(obj, u, tol)
        % CONTINUITYAT  Parametric continuity order at parameter u.
        %   For an interior knot of multiplicity s, continuity is C^(p-s).
        %   Outside interior knots, returns Inf.
            if nargin < 3 || isempty(tol), tol = 1e-12; end
            dom = obj.domain;
            if abs(u - dom(1)) < tol || abs(u - dom(2)) < tol
                c = Inf;
                return;
            end
            s = obj.knotMultiplicity(u, tol);
            if s == 0
                c = Inf;
            else
                c = obj.p - s;
            end
        end

        function tf = isClosed(obj, tol)
        % ISCLOSED  True if curve endpoints coincide geometrically.
            if nargin < 2 || isempty(tol), tol = 1e-9; end
            C0 = obj.evaluate(obj.domain(1));
            C1 = obj.evaluate(obj.domain(2));
            tf = norm(C1 - C0) <= tol;
        end

        function ug = grevilleAbscissae(obj)
        % GREVILLEABSCISSAE  Greville abscissae of the B-spline basis.
            n = obj.n; p = obj.p; U = obj.U;
            ug = zeros(n+1,1);
            if p == 0
                ug(:) = U(1:n+1);
                return;
            end
            for i = 0:n
                ug(i+1) = sum(U(i+2:i+p+1)) / p;
            end
        end

        function [Q, errMax, info] = removeKnot(obj, u, numRemove, tol, nSample)
        % REMOVEKNOT  Tolerance-controlled numerical knot removal.
        %
        %   [Q, errMax, info] = removeKnot(u)
        %   [Q, errMax, info] = removeKnot(u, numRemove, tol, nSample)
        %
        %   Attempts to remove up to numRemove copies of knot u while keeping
        %   the geometric deviation below tol. This is a global numerical
        %   removal on the homogeneous curve and is intended as a robust
        %   cleanup tool when exact local removal logic is not available.
            if nargin < 3 || isempty(numRemove), numRemove = 1; end
            if nargin < 4 || isempty(tol), tol = 1e-8; end
            if nargin < 5 || isempty(nSample), nSample = max(200, 20*(obj.n+1)); end
            if numRemove < 0 || numRemove ~= floor(numRemove)
                error('NURBSCurve:removeKnot', 'numRemove must be a nonnegative integer.');
            end

            Q = geom.NURBSCurve(obj.P, obj.p, obj.U, obj.W);
            errMax = 0;
            info = struct('requested', numRemove, 'removed', 0, 'finalMultiplicity', Q.knotMultiplicity(u), ...
                          'accepted', false, 'history', []);

            for k = 1:numRemove
                s = Q.knotMultiplicity(u);
                if s == 0
                    break;
                end
                Ucand = Q.U;
                idx = find(abs(Ucand - u) < 1e-12, 1, 'last');
                Ucand(idx) = [];
                % Build candidate by fitting the homogeneous curve on the new knot vector.
                Qcand = geom.NURBSCurve.fitHomogeneousToKnotVector(Q, Ucand, nSample);
                % Measure geometric deviation against current curve on the active domain.
                us = linspace(max(Q.domain(1), Qcand.domain(1)), min(Q.domain(2), Qcand.domain(2)), nSample);
                if isempty(us) || numel(us) < 2
                    break;
                end
                Pold = Q.evaluate(us);
                Pnew = Qcand.evaluate(us);
                e = sqrt(sum((Pnew - Pold).^2, 2));
                ek = max(e);
                info.history = [info.history; [k, ek]]; %#ok<AGROW>
                if ek <= tol
                    Q = Qcand;
                    errMax = ek;
                    info.removed = info.removed + 1;
                else
                    errMax = ek;
                    break;
                end
            end

            info.finalMultiplicity = Q.knotMultiplicity(u);
            info.accepted = info.removed > 0;
        end

        function [parts, breaks] = subdivideUniform(obj)
        % SUBDIVIDEUNIFORM  Split into knot spans as curve pieces.
            uk = unique(obj.U);
            uk = uk(uk >= obj.domain(1) & uk <= obj.domain(2));
            parts = {};
            breaks = uk(:);
            if numel(uk) < 2
                return;
            end
            Cwork = geom.NURBSCurve(obj.P, obj.p, obj.U, obj.W);
            parts = cell(numel(uk)-1, 1);
            for i = 1:numel(uk)-1
                if uk(i+1) <= uk(i), continue; end
                if i == 1
                    [left, right] = Cwork.split(uk(i+1));
                    parts{i} = left;
                    Cwork = right;
                elseif i < numel(uk)-1
                    local = (uk(i+1) - uk(i)) / (uk(end) - uk(i));
                    [left, right] = Cwork.split(local);
                    parts{i} = left;
                    Cwork = right;
                else
                    parts{i} = Cwork;
                end
            end
        end

        function S = sample(obj, u)
        % SAMPLE  Alias for evaluate; convenient for generic geometry code.
            S = obj.evaluate(u);
        end
    end

    methods (Static)
        function Q = fitHomogeneousToKnotVector(Csrc, Ucand, nSample)
        % FITHOMOGENEOUSTOKNOTVECTOR  Reconstruct homogeneous control points
        % for a target knot vector and the same degree by least squares.
            p = Csrc.p;
            n2 = numel(Ucand) - p - 2;
            if n2 < p
                error('NURBSCurve:fitHomogeneousToKnotVector', ...
                    'Target knot vector is invalid for degree p.');
            end
            % Sample on the overlapping active domain to avoid endpoint ambiguity.
            dom = [max(Csrc.U(Csrc.p+1), Ucand(p+1)), min(Csrc.U(end-Csrc.p), Ucand(end-p))];
            if dom(2) <= dom(1)
                error('NURBSCurve:fitHomogeneousToKnotVector', ...
                    'No overlapping active domain between source and candidate.');
            end
            us = linspace(dom(1), dom(2), nSample)';
            H = zeros(nSample, 4);
            Ps = Csrc.evaluate(us);
            Ws = zeros(nSample,1);
            for i = 1:nSample
                CK = Csrc.derivatives(us(i), 0);
                % derive weight by evaluating homogeneous denominator directly
                span = geom.BasisFunctions.FindSpan(Csrc.n, p, us(i), Csrc.U);
                N = geom.BasisFunctions.BasisFuns(span, us(i), p, Csrc.U);
                w = 0;
                for j = 0:p
                    idx = span - p + j;
                    w = w + N(j+1) * Csrc.W(idx);
                end
                Ws(i) = w;
                H(i,1:3) = Ps(i,:) * w;
                H(i,4) = w;
            end
            A = geom.NURBSCurve.buildBasisMatrixFromKnotVector(Ucand, p, us);
            Pw = A \ H;
            % Enforce positive weights.
            if any(Pw(:,4) <= 0)
                Pw(:,4) = max(Pw(:,4), 1e-12);
            end
            P = Pw(:,1:3) ./ Pw(:,4);
            Q = geom.NURBSCurve(P, p, Ucand, Pw(:,4));
        end

        function A = buildBasisMatrixFromKnotVector(U, p, us)
            n = numel(U) - p - 2;
            A = zeros(numel(us), n+1);
            for k = 1:numel(us)
                uk = min(max(us(k), U(p+1)), U(end-p));
                span = geom.BasisFunctions.FindSpan(n, p, uk, U);
                N = geom.BasisFunctions.BasisFuns(span, uk, p, U);
                A(k, span-p+1:span+1) = N;
            end
        end
    end

end