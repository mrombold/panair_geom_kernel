classdef NURBSCurve < handle
% NURBSCURVE  Rational B-spline / NURBS curve kernel.
%
% Core capabilities:
%   - exact rational evaluation
%   - exact rational derivatives
%   - tangent / curvature / Frenet frame
%   - arc length
%   - exact knot insertion / refinement / split / Bezier decomposition
%   - degree elevation
%   - approximate degree reduction in homogeneous space
%   - closest-point projection
%   - global interpolation / least-squares fitting
%
% Notes:
%   - Control points are stored in Cartesian form P plus weights W.
%   - Homogeneous control points are formed internally as Pw = [w*x w*y w*z w].
%   - Parameter domain is [U(p+1), U(end-p)] for clamped/open curves.

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
                U = U(:).';
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
        function tf = validate(obj)
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
            tf = true;
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
            u = u(:).';
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
            d = min(max(0, floor(d)), obj.p);

            span  = geom.BasisFunctions.FindSpan(obj.n, obj.p, u, obj.U);
            Nders = geom.BasisFunctions.DersBasisFuns(span, u, obj.p, d, obj.U);

            Aders = zeros(d+1, 3);
            wders = zeros(d+1, 1);

            for j = 0:obj.p
                idx = span - obj.p + j;
                wj = obj.W(idx);
                Pj = obj.P(idx,:);
                for kk = 0:d
                    Aders(kk+1,:) = Aders(kk+1,:) + Nders(kk+1,j+1) * wj * Pj;
                    wders(kk+1)   = wders(kk+1)   + Nders(kk+1,j+1) * wj;
                end
            end

            CK = zeros(d+1, 3);
            CK(1,:) = Aders(1,:) / wders(1);
            for kk = 1:d
                v = Aders(kk+1,:);
                for i = 1:kk
                    v = v - nchoosek(kk,i) * wders(i+1) * CK(kk-i+1,:);
                end
                CK(kk+1,:) = v / wders(1);
            end
        end

        function D = derivative(obj, u, k)
            if nargin < 3 || isempty(k), k = 1; end
            u = u(:).';
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

            u0 = obj.clamp(u0);
            u1 = obj.clamp(u1);
            if u1 < u0
                tmp = u0; u0 = u1; u1 = tmp;
            end

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
            if s(end) < eps
                u = linspace(obj.domain(1), obj.domain(2), nPts);
                return;
            end
            s = s / s(end);
            u = interp1(s, us(:), linspace(0,1,nPts)', 'pchip').';
        end
    end

    methods
        function C2 = insertKnot(obj, u, r)
        % Exact knot insertion by repeated single insertion.
            if nargin < 3 || isempty(r), r = 1; end
            if r < 0 || r ~= floor(r)
                error('NURBSCurve:insertKnot', 'r must be a nonnegative integer.');
            end

            C2 = geom.NURBSCurve(obj.P, obj.p, obj.U, obj.W);
            for ii = 1:r
                C2 = C2.insertKnotOnce(obj.clamp(u));
            end
        end

        function C2 = refine(obj, X)
        % Exact curve knot refinement by repeated single-knot insertion.
            X = sort(X(:).');
            C2 = geom.NURBSCurve(obj.P, obj.p, obj.U, obj.W);
            for ii = 1:numel(X)
                C2 = C2.insertKnotOnce(X(ii));
            end
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
            rows = find(abs(U2 - u0) < tol);
            first = rows(1);
            last  = rows(end);

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

        function segs = decomposeToBezier(obj)
        % Compatibility alias for older demos expecting a cell array of curves.
            parts = obj.decomposeBezier();
            segs = cell(size(parts));
            for i = 1:numel(parts)
                segs{i} = parts{i}.curve;
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
            Qt = (T * Ph.').';
            Q = Qt(:,1:3) ./ Qt(:,4);
            C2 = geom.NURBSCurve(Q, obj.p, obj.U, obj.W);
        end

        function C2 = translate(obj, v)
            Q = obj.P + repmat(v(:).', size(obj.P,1), 1);
            C2 = geom.NURBSCurve(Q, obj.p, obj.U, obj.W);
        end

        function bbox = boundingBox(obj, nPts)
            if nargin < 2 || isempty(nPts), nPts = 200; end
            pts = obj.evaluate(linspace(obj.domain(1), obj.domain(2), nPts));
            bbox = [min(pts); max(pts)].';
        end
    end

    methods
        function H = evaluateHomogeneous(obj, u)
        % Exact homogeneous 4D B-spline evaluation.
            u = u(:).';
            H = zeros(numel(u), 4);

            for k = 1:numel(u)
                uk = obj.clamp(u(k));
                span = geom.BasisFunctions.FindSpan(obj.n, obj.p, uk, obj.U);
                N = geom.BasisFunctions.BasisFuns(span, uk, obj.p, obj.U);
                Cw = zeros(1,4);
                for j = 0:obj.p
                    idx = span - obj.p + j;
                    Cw = Cw + N(j+1) * obj.Pw(idx,:);
                end
                H(k,:) = Cw;
            end
        end

        function g = grevilleAbscissae(obj)
            g = geom.NURBSCurve.grevilleFromKnotVector(obj.U, obj.p);
        end

        function tf = isClosed(obj, tol)
            if nargin < 2 || isempty(tol), tol = 1e-10; end
            a = obj.evaluate(obj.domain(1));
            b = obj.evaluate(obj.domain(2));
            tf = norm(a - b) <= tol;
        end

        function k = continuityAt(obj, u, tol)
            if nargin < 3 || isempty(tol), tol = 1e-12; end
            if abs(u - obj.domain(1)) < tol || abs(u - obj.domain(2)) < tol
                k = -inf;
                return;
            end
            s = obj.knotMultiplicity(u, tol);
            k = obj.p - s;
        end

        function [C2, removed, maxErr] = removeKnot(obj, u, numRemove, tol, nSample)
        % Remove one or more knot copies by reconstruction / validation.
            if nargin < 3 || isempty(numRemove), numRemove = 1; end
            if nargin < 4 || isempty(tol), tol = 1e-10; end
            if nargin < 5 || isempty(nSample), nSample = max(200, 20*(obj.n+1)); end

            Ccur = geom.NURBSCurve(obj.P, obj.p, obj.U, obj.W);
            removed = 0;
            maxErr = 0;

            for it = 1:numRemove
                [ok, Cnext, err] = Ccur.tryRemoveOneKnot(u, tol, nSample);
                if ~ok
                    break;
                end
                Ccur = Cnext;
                removed = removed + 1;
                maxErr = max(maxErr, err);
            end

            C2 = Ccur;
        end

        function [C2, maxErr] = reduceDegree(obj, numTimes, tol, nSample)
        % Degree reduction by Bezier decomposition + homogeneous LSQ reduction.
            if nargin < 2 || isempty(numTimes), numTimes = 1; end
            if nargin < 3 || isempty(tol), tol = inf; end
            if nargin < 4 || isempty(nSample), nSample = max(300, 30*(obj.n+1)); end

            if numTimes < 0 || numTimes ~= floor(numTimes)
                error('NURBSCurve:reduceDegree', 'numTimes must be a nonnegative integer.');
            end

            Ccur = geom.NURBSCurve(obj.P, obj.p, obj.U, obj.W);
            maxErr = 0;

            for step = 1:numTimes
                if Ccur.p <= 0
                    error('NURBSCurve:reduceDegree', 'Cannot reduce degree below zero.');
                end

                parts = Ccur.decomposeBezier();
                ph = Ccur.p - 1;
                nSeg = numel(parts);

                PwAll = [];
                breaks = zeros(nSeg+1,1);
                segErr = 0;

                for s = 1:nSeg
                    Cseg = parts{s}.curve;
                    PwSeg = Cseg.Pw;
                    [PwRed, eSeg] = geom.NURBSCurve.reduceBezierHomogeneousLSQ(PwSeg);
                    segErr = max(segErr, eSeg);

                    if s == 1
                        PwAll = PwRed;
                        breaks(1) = parts{s}.u0;
                    else
                        PwAll = [PwAll; PwRed(2:end,:)]; %#ok<AGROW>
                    end
                    breaks(s+1) = parts{s}.u1;
                end

                Unew = [repmat(breaks(1), 1, ph+1)];
                for s = 2:nSeg
                    Unew = [Unew, repmat(breaks(s), 1, ph)]; %#ok<AGROW>
                end
                Unew = [Unew, repmat(breaks(end), 1, ph+1)];

                W2 = PwAll(:,4);
                if any(W2 <= 0)
                    error('NURBSCurve:reduceDegree', 'Degree reduction produced nonpositive weights.');
                end

                P2 = PwAll(:,1:3) ./ W2;
                Cnext = geom.NURBSCurve(P2, ph, Unew, W2);

                us = geom.NURBSCurve.validationParams(Ccur.U, Cnext.U, nSample);
                err = max(vecnorm(Ccur.evaluate(us) - Cnext.evaluate(us), 2, 2));

                maxErr = max(maxErr, max(segErr, err));
                if maxErr > tol
                    error('NURBSCurve:reduceDegree', ...
                        'Degree reduction exceeded tolerance: %g > %g.', maxErr, tol);
                end

                Ccur = Cnext;
            end

            C2 = Ccur;
        end
    end

    methods
        function [u_c, pt_c, d_c] = closestPoint(obj, P, u0)
            P = P(:).';
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

                if abs(f2) < 1e-14
                    break;
                end

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
                plot3(obj.P(:,1), obj.P(:,2), obj.P(:,3), 'o--', ...
                    'Color', [0.6 0.6 0.6], 'MarkerFaceColor', 'w', 'MarkerSize', 5);
            end

            if opt.ShowKnots
                uk = unique(obj.U);
                uk = uk(uk > obj.domain(1) & uk < obj.domain(2));
                if ~isempty(uk)
                    kp = obj.evaluate(uk);
                    plot3(kp(:,1), kp(:,2), kp(:,3), 'k^', 'MarkerFaceColor', 'k');
                end
            end

            axis equal;
            grid on;
            xlabel('X'); ylabel('Y'); zlabel('Z');
        end
    end

    methods (Static)
        function U = averagingKnotVector(params, p)
            params = params(:).';
            n = numel(params) - 1;
    
            U = zeros(1, n + p + 2);
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
                u = linspace(0,1,m+1).';
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
                error('NURBSCurve:globalInterp', ...
                    'Degree exceeds number of data points minus one.');
            end
    
            u = geom.NURBSCurve.parameterizeData(Q, method);
            U = geom.NURBSCurve.averagingKnotVector(u, p);
    
            A = zeros(n+1, n+1);
            for r = 1:n+1
                span = geom.BasisFunctions.FindSpan(n, p, u(r), U);
                N = geom.BasisFunctions.BasisFuns(span, u(r), p, U);
                cols = (span-p):span;
                A(r, cols) = N;
            end
    
            P = A \ Q;
    
            if size(P,1) ~= n+1 || size(P,2) ~= 3
                error('NURBSCurve:globalInterp', ...
                    'Interpolation solve returned a control-point array of size %dx%d, expected %dx3.', ...
                    size(P,1), size(P,2), n+1);
            end
    
            C = geom.NURBSCurve(P, p, U, ones(n+1,1));
        end

            function C = globalLeastSquaresFit(Q, p, nCtrl, method)
        % Global least-squares fit of data Q with nCtrl control points.
    
            if nargin < 4 || isempty(method), method = 'centripetal'; end
    
            Q = geom.NURBSCurve.ensure3D(Q);
    
            m = size(Q,1) - 1;     % data index
            n = nCtrl - 1;         % control-point index
    
            if nCtrl < p + 1
                error('NURBSCurve:globalLeastSquaresFit', ...
                    'Need at least p+1 control points.');
            end
            if m < n
                error('NURBSCurve:globalLeastSquaresFit', ...
                    'Need at least as many data points as control points.');
            end
    
            u = geom.NURBSCurve.parameterizeData(Q, method);
    
            % Averaged knot vector for LSQ fitting
            U = zeros(1, n + p + 2);
            U(1:p+1) = 0;
            U(end-p:end) = 1;
    
            if (n-p) >= 1
                d = (m + 1) / (n - p + 1);
                for j = 1:(n-p)
                    jd = j * d;
                    i = floor(jd);
                    alpha = jd - i;
    
                    i = max(1, min(i, m));
                    U(j+p+1) = (1-alpha) * u(i) + alpha * u(i+1);
                end
            end
    
            Nmat = zeros(m+1, n+1);
            for r = 1:m+1
                span = geom.BasisFunctions.FindSpan(n, p, u(r), U);
                N = geom.BasisFunctions.BasisFuns(span, u(r), p, U);
                cols = (span-p):span;
                Nmat(r, cols) = N;
            end
    
            P = Nmat \ Q;
    
            if size(P,1) ~= (n+1) || size(P,2) ~= 3
                error('NURBSCurve:globalLeastSquaresFit', ...
                    'Least-squares solve returned size %dx%d, expected %dx3.', ...
                    size(P,1), size(P,2), n+1);
            end
    
            C = geom.NURBSCurve(P, p, U, ones(n+1,1));
    end

        function C = line(P0, P1)
            P0 = geom.NURBSCurve.ensure3D(P0);
            P1 = geom.NURBSCurve.ensure3D(P1);
            P = [P0(1,:); P1(1,:)];
            C = geom.NURBSCurve(P, 1, [0 0 1 1], [1; 1]);
        end

        function C = quadraticConic(P0, P1, P2, w1)
            if nargin < 4 || isempty(w1), w1 = 1; end
            if w1 <= 0
                error('NURBSCurve:quadraticConic', 'Middle weight must be positive.');
            end
            P0 = geom.NURBSCurve.ensure3D(P0);
            P1 = geom.NURBSCurve.ensure3D(P1);
            P2 = geom.NURBSCurve.ensure3D(P2);
            P = [P0(1,:); P1(1,:); P2(1,:)];
            W = [1; w1; 1];
            U = [0 0 0 1 1 1];
            C = geom.NURBSCurve(P, 2, U, W);
        end
    end

    methods (Static)
        function g = grevilleFromKnotVector(U, p)
            n = numel(U) - p - 2;
            g = zeros(n+1,1);
            if p == 0
                g(:) = U(1:n+1).';
                return;
            end
            for i = 1:n+1
                g(i) = sum(U(i+1:i+p)) / p;
            end
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

        function [PwRed, maxErr] = reduceBezierHomogeneousLSQ(Pw)
            p = size(Pw,1) - 1;
            if p < 1
                error('NURBSCurve:reduceBezierHomogeneousLSQ', 'Bezier degree must be at least 1.');
            end

            ph = p - 1;

            E = zeros(p+1, ph+1);
            E(1,1) = 1;
            E(end,end) = 1;
            for i = 1:p-1
                alpha = i / p;
                E(i+1, i)   = alpha;
                E(i+1, i+1) = 1 - alpha;
            end

            PwRed = zeros(ph+1, size(Pw,2));
            PwRed(1,:)   = Pw(1,:);
            PwRed(end,:) = Pw(end,:);

            if ph > 1
                A = E(:,2:end-1);
                rhs = Pw - E(:,1)*PwRed(1,:) - E(:,end)*PwRed(end,:);
                X = A \ rhs;
                PwRed(2:end-1,:) = X;
            end

            PwElev = E * PwRed;
            maxErr = max(vecnorm(Pw - PwElev, 2, 2));
        end

        function us = validationParams(U1, U2, nSample)
            if nargin < 3 || isempty(nSample), nSample = 400; end
            b = unique([U1(:); U2(:)]).';
            mids = 0.5 * (b(1:end-1) + b(2:end));
            a = max(U1(1), U2(1));
            z = min(U1(end), U2(end));
            us = unique([linspace(a, z, nSample), b, mids]);
            us = us(~isnan(us));
            us = us(:).';
        end

        function U = normalizeKnotVector(U)
            a = U(1);
            b = U(end);
            if abs(b-a) < eps
                U = zeros(size(U));
            else
                U = (U - a) / (b - a);
            end
        end

        function P = ensure3D(P)
            if size(P,2) == 2
                P = [P, zeros(size(P,1),1)];
            elseif size(P,2) == 3
                % ok
            elseif isvector(P) && numel(P) == 2
                P = [reshape(P,1,2), 0];
            elseif isvector(P) && numel(P) == 3
                P = reshape(P,1,3);
            else
                error('NURBSCurve:ensure3D', 'Input must be Nx2, Nx3, 1x2, or 1x3.');
            end
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
        function C2 = insertKnotOnce(obj, u)
        % Exact single knot insertion (Piegl & Tiller, 1-based MATLAB form).
    
            p  = obj.p;
            U  = obj.U;
            Pw = obj.Pw;
            n  = obj.n;
    
            u = obj.clamp(u);
    
            % span is used in the same 1-based convention as the rest of this class
            span = geom.BasisFunctions.FindSpan(n, p, u, U);
            s    = obj.knotMultiplicity(u, 1e-12);
    
            if s >= p
                C2 = geom.NURBSCurve(obj.P, obj.p, obj.U, obj.W);
                return;
            end
    
            % New knot vector: insert u after U(span)
            Up = zeros(1, numel(U) + 1);
            Up(1:span) = U(1:span);
            Up(span+1) = u;
            Up(span+2:end) = U(span+1:end);
    
            % New homogeneous control net
            Qw = zeros(size(Pw,1) + 1, 4);
    
            % Unaffected left block
            Qw(1:span-p, :) = Pw(1:span-p, :);
    
            % Unaffected right block
            Qw(span-s+1:end, :) = Pw(span-s:end, :);
    
            % Updated middle block
            for j = (span-p+1):(span-s)
                denom = U(j+p) - U(j);
                if abs(denom) < 1e-14
                    alpha = 0.0;
                else
                    alpha = (u - U(j)) / denom;
                end
                Qw(j,:) = alpha * Pw(j,:) + (1-alpha) * Pw(j-1,:);
            end
    
            if any(diff(Up) < -1e-14)
                error('NURBSCurve:insertKnotOnce', 'Internal error: generated knot vector is not nondecreasing.');
            end
    
            W2 = Qw(:,4);
            if any(W2 <= 0)
                error('NURBSCurve:insertKnotOnce', 'Knot insertion produced nonpositive weights.');
            end
    
            P2 = Qw(:,1:3) ./ W2;
            C2 = geom.NURBSCurve(P2, p, Up, W2);
        end

        function [ok, Ccand, maxErr] = tryRemoveOneKnot(obj, u, tol, nSample)
            ok = false;
            Ccand = obj;
            maxErr = inf;
    
            tolK = 1e-12;
            idx = find(abs(obj.U - u) < tolK, 1, 'first');
            if isempty(idx)
                return;
            end
    
            % Do not remove end-clamping copies here
            if idx <= obj.p+1 || idx >= numel(obj.U)-obj.p
                return;
            end
    
            Ucand = obj.U;
            Ucand(idx) = [];
    
            n2 = numel(Ucand) - obj.p - 2;
            if n2 < obj.p
                return;
            end
    
            % Greville points for the candidate knot vector
            g = geom.NURBSCurve.grevilleFromKnotVector(Ucand, obj.p);
    
            % Sample original curve in homogeneous space at those parameters
            H = obj.evaluateHomogeneous(g);
    
            % Build square interpolation matrix for the candidate curve
            A = zeros(n2+1, n2+1);
            for rr = 1:n2+1
                span = geom.BasisFunctions.FindSpan(n2, obj.p, g(rr), Ucand);
                N = geom.BasisFunctions.BasisFuns(span, g(rr), obj.p, Ucand);
    
                cols = (span-obj.p):span;   % 1-based convention used elsewhere
                if cols(1) < 1 || cols(end) > (n2+1)
                    return;
                end
    
                A(rr, cols) = N;
            end
    
            if size(A,1) ~= size(A,2)
                return;
            end
            if rcond(A) < 1e-13
                return;
            end
    
            Qw = A \ H;
            if any(Qw(:,4) <= 0)
                return;
            end
    
            Q = Qw(:,1:3) ./ Qw(:,4);
            Ctest = geom.NURBSCurve(Q, obj.p, Ucand, Qw(:,4));
    
            us = geom.NURBSCurve.validationParams(obj.U, Ucand, nSample);
            maxErr = max(vecnorm(obj.evaluate(us) - Ctest.evaluate(us), 2, 2));
    
            if maxErr <= tol
                ok = true;
                Ccand = Ctest;
            end
        end

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
end