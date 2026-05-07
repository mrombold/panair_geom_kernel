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
%   - weighted and constrained least-squares fitting
%   - multiple parameterization methods
%
% Notes:
%   - Control points are stored in Cartesian form P plus weights W.
%   - Homogeneous control points are formed internally as Pw = [w*x w*y w*z w].
%   - Parameter domain is [U(p+1), U(end-p)] for clamped/open curves.

    properties (Access = private)
        Pw_
        U_
        p_
    end

    properties (Dependent)
        P
        W
        U
        p
        n
        m
        domain
        Pw
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

            obj.Pw_ = [P .* W, W];
            obj.U_ = U;
            obj.p_ = p;
            obj.validate();
        end

        function v = get.P(obj)
            w = obj.Pw_(:,4);
            v = obj.Pw_(:,1:3) ./ w;
        end
        
        function set.P(obj, P)
            [npts, dim] = size(P);
            if dim == 2
                P = [P, zeros(npts,1)];
            elseif dim ~= 3
                error('NURBSCurve:set.P', 'P must be Nx2 or Nx3.');
            end
            if npts ~= size(obj.Pw_,1)
                error('NURBSCurve:set.P', 'P size must match existing control-point count.');
            end
            w = obj.Pw_(:,4);
            obj.Pw_ = [P .* w, w];
        end
        
        function v = get.W(obj)
            v = obj.Pw_(:,4);
        end
        
        function set.W(obj, W)
            W = W(:);
            if numel(W) ~= size(obj.Pw_,1)
                error('NURBSCurve:set.W', 'W must match the number of control points.');
            end
            if any(W <= 0)
                error('NURBSCurve:set.W', 'Weights must be strictly positive.');
            end
            P = obj.P;
            obj.Pw_ = [P .* W, W];
        end
        
        function v = get.U(obj)
            v = obj.U_;
        end
        
        function set.U(obj, U)
            U = U(:).';
            if any(diff(U) < 0)
                error('NURBSCurve:set.U', 'Knot vector must be nondecreasing.');
            end
            obj.U_ = U;
        end
        
        function v = get.p(obj)
            v = obj.p_;
        end
        
        function set.p(obj, p)
            if ~isscalar(p) || p < 0 || p ~= floor(p)
                error('NURBSCurve:set.p', 'Degree p must be a nonnegative integer.');
            end
            obj.p_ = p;
        end
        
        function v = get.n(obj)
            v = size(obj.Pw_,1) - 1;
        end
        
        function v = get.m(obj)
            v = numel(obj.U_) - 1;
        end
        
        function v = get.domain(obj)
            v = [obj.U_(obj.p_+1), obj.U_(end-obj.p_)];
        end
        
        function v = get.Pw(obj)
            v = obj.Pw_;
        end

    end

    methods
        function tf = validate(obj)
            if size(obj.Pw_,2) ~= 4
                error('NURBSCurve:Validate', 'Pw_ must be [n+1 x 4].');
            end
            if numel(obj.U_) ~= obj.n + obj.p_ + 2
                error('NURBSCurve:Validate', 'Knot-vector length is inconsistent with n and p.');
            end
            if any(diff(obj.U_) < 0)
                error('NURBSCurve:Validate', 'Knot vector must be nondecreasing.');
            end
            if any(obj.Pw_(:,4) <= 0)
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
        %INSERTKNOT Exact curve knot insertion.
        % Implements Piegl & Tiller Algorithm A5.1, CurveKnotIns,
        % in homogeneous coordinates.
        
            if nargin < 3 || isempty(r), r = 1; end
        
            r = floor(r);
            if r < 0
                error('NURBSCurve:insertKnot', ...
                    'r must be a nonnegative integer.');
            end
        
            if r == 0
                C2 = geom.NURBSCurve(obj.P, obj.p, obj.U, obj.W);
                return;
            end
        
            u = obj.clamp(u);
        
            s = obj.knotMultiplicity(u, 1e-12);
            if u > obj.domain(1) + 1e-12 && u < obj.domain(2) - 1e-12
                if s + r > obj.p
                    error('NURBSCurve:insertKnot', ...
                        'Cannot insert knot: multiplicity s+r exceeds degree p.');
                end
            end
        
            [Qw, Uq] = geom.NURBSCurve.curveKnotInsertHomogeneous( ...
                obj.p, obj.U, obj.Pw, u, r);
        
            W2 = Qw(:,4);
            if any(W2 <= 0)
                error('NURBSCurve:insertKnot', ...
                    'Knot insertion produced nonpositive weights.');
            end
        
            P2 = Qw(:,1:3) ./ W2;
            C2 = geom.NURBSCurve(P2, obj.p, Uq, W2);
        end

        function C2 = refine(obj, X)
        %REFINE Exact curve knot refinement by repeated A5.1 insertion.
        
            if nargin < 2 || isempty(X)
                C2 = geom.NURBSCurve(obj.P, obj.p, obj.U, obj.W);
                return;
            end
        
            X = sort(X(:).');
            C2 = geom.NURBSCurve(obj.P, obj.p, obj.U, obj.W);
        
            for ii = 1:numel(X)
                C2 = C2.insertKnot(X(ii), 1);
            end
        end

        function [Cleft, Cright] = split(obj, u0)
        %SPLIT Exact curve split by full-multiplicity knot insertion.
        %
        % Keeps absolute knot domains:
        %   Cleft.domain  = [old_start, u0]
        %   Cright.domain = [u0, old_end]
        
            u0 = obj.clamp(u0);
            tol = 1e-12;
        
            if u0 <= obj.domain(1)+tol || u0 >= obj.domain(2)-tol
                error('NURBSCurve:split', ...
                    'Split parameter must be strictly interior to the active domain.');
            end
        
            s = obj.knotMultiplicity(u0, tol);
            if s > obj.p
                error('NURBSCurve:split', ...
                    'Split knot multiplicity exceeds degree.');
            end
        
            if s < obj.p
                Cref = obj.insertKnot(u0, obj.p - s);
            else
                Cref = geom.NURBSCurve(obj.P, obj.p, obj.U, obj.W);
            end
        
            U2 = Cref.U;
            rows = find(abs(U2 - u0) < tol);
        
            if numel(rows) < obj.p
                error('NURBSCurve:split', ...
                    'Failed to create full-multiplicity split knot.');
            end
        
            first = rows(1);
            last  = rows(end);
        
            idxSplit = last - obj.p;
        
            P1 = Cref.P(1:idxSplit,:);
            W1 = Cref.W(1:idxSplit);
            U1 = [U2(1:last), u0];
        
            P2 = Cref.P(idxSplit:end,:);
            W2 = Cref.W(idxSplit:end);
            U3 = [u0, U2(first:end)];
        
            % Clamp child endpoint knots at the split while preserving absolute domains.
            U1(end-obj.p:end) = u0;
            U3(1:obj.p+1) = u0;
        
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
        %REMOVEKNOT Exact curve knot removal using Piegl & Tiller Algorithm A5.8.
        %
        % Attempts to remove up to numRemove copies of knot u.
        % Works in homogeneous coordinates.
        %
        % Outputs:
        %   C2      resulting curve
        %   removed number of copies actually removed
        %   maxErr  sampled Cartesian validation error after removals
        
            if nargin < 3 || isempty(numRemove), numRemove = 1; end
            if nargin < 4 || isempty(tol), tol = 1e-10; end
            if nargin < 5 || isempty(nSample), nSample = max(200, 20*(obj.n+1)); end
        
            numRemove = floor(numRemove);
            if numRemove < 0
                error('NURBSCurve:removeKnot', ...
                    'numRemove must be a nonnegative integer.');
            end
        
            if numRemove == 0
                C2 = geom.NURBSCurve(obj.P, obj.p, obj.U, obj.W);
                removed = 0;
                maxErr = 0;
                return;
            end
        
            Ccur = geom.NURBSCurve(obj.P, obj.p, obj.U, obj.W);
            removed = 0;
        
            for rr = 1:numRemove
                [ok, Unew, Pwnew] = geom.NURBSCurve.removeOneKnotA58( ...
                    Ccur.p, Ccur.U, Ccur.Pw, u, tol);
        
                if ~ok
                    break;
                end
        
                Wnew = Pwnew(:,4);
                if any(Wnew <= 0)
                    break;
                end
        
                Pnew = Pwnew(:,1:3) ./ Wnew;
                Ccur = geom.NURBSCurve(Pnew, Ccur.p, Unew, Wnew);
                removed = removed + 1;
            end
        
            C2 = Ccur;
        
            if removed == 0
                maxErr = 0;
            else
                us = geom.NURBSCurve.validationParams(obj.U, C2.U, nSample);
                maxErr = max(vecnorm(obj.evaluate(us) - C2.evaluate(us), 2, 2));
            end
        end


        function [C2, maxErr] = reduceDegree(obj, numTimes, tol, nSample)
        %REDUCEDEGREE Exact/tolerance-based degree reduction.
        %
        % Reduces degree by one per pass using exact Bezier degree-reduction
        % identities in homogeneous coordinates.
        %
        % If the curve is not reducible within tol, this errors.
        % This replaces the old LSQ approximate degree reduction.
        
            if nargin < 2 || isempty(numTimes), numTimes = 1; end
            if nargin < 3 || isempty(tol), tol = 1e-10; end
            if nargin < 4 || isempty(nSample), nSample = max(400, 40*(obj.n+1)); end
        
            numTimes = floor(numTimes);
            if numTimes < 0
                error('NURBSCurve:reduceDegree', ...
                    'numTimes must be a nonnegative integer.');
            end
        
            Ccur = geom.NURBSCurve(obj.P, obj.p, obj.U, obj.W);
            maxErr = 0;
        
            for step = 1:numTimes
                if Ccur.p <= 0
                    error('NURBSCurve:reduceDegree', ...
                        'Cannot reduce degree below zero.');
                end
        
                pOld = Ccur.p;
                pNew = pOld - 1;
        
                parts = Ccur.decomposeBezier();
                nSeg = numel(parts);
        
                if nSeg == 0
                    error('NURBSCurve:reduceDegree', ...
                        'Bezier decomposition failed.');
                end
        
                PwAll = [];
                breaks = zeros(nSeg+1, 1);
                localErr = 0;
        
                for s = 1:nSeg
                    Cseg = parts{s}.curve;
                    PwSeg = Cseg.Pw;
        
                    [PwRed, eSeg] = geom.NURBSCurve.reduceBezierHomogeneousExact(PwSeg, tol);
                    localErr = max(localErr, eSeg);
        
                    if s == 1
                        PwAll = PwRed;
                        breaks(1) = parts{s}.u0;
                    else
                        PwAll = [PwAll; PwRed(2:end,:)]; %#ok<AGROW>
                    end
        
                    breaks(s+1) = parts{s}.u1;
                end
        
                Unew = repmat(breaks(1), 1, pNew+1);
                for s = 2:nSeg
                    Unew = [Unew, repmat(breaks(s), 1, pNew)]; %#ok<AGROW>
                end
                Unew = [Unew, repmat(breaks(end), 1, pNew+1)];
        
                W2 = PwAll(:,4);
                if any(W2 <= 0)
                    error('NURBSCurve:reduceDegree', ...
                        'Degree reduction produced nonpositive weights.');
                end
        
                P2 = PwAll(:,1:3) ./ W2;
                Cnext = geom.NURBSCurve(P2, pNew, Unew, W2);
        
                us = geom.NURBSCurve.validationParams(Ccur.U, Cnext.U, nSample);
                err = max(vecnorm(Ccur.evaluate(us) - Cnext.evaluate(us), 2, 2));
        
                maxErr = max(maxErr, max(localErr, err));
        
                if maxErr > tol
                    error('NURBSCurve:reduceDegree', ...
                        'Curve is not degree-reducible within tolerance: %.3e > %.3e.', ...
                        maxErr, tol);
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
            view(45,30)
            axis equal;
            grid on;
            xlabel('X'); ylabel('Y'); zlabel('Z');
        end
    end

    methods (Static)
        function U = averagingKnotVector(params, p)
            params = params(:).';
            n = numel(params) - 1;

            U = zeros(1, n+p+2);
            U(1:p+1) = params(1);
            U(end-p:end) = params(end);

            for j = 1:(n-p)
                U(j+p+1) = sum(params(j+1:j+p)) / p;
            end
        end

        function u = parameterizeData(Q, method, varargin)
            if nargin < 2 || isempty(method), method = 'centripetal'; end
            Q = geom.NURBSCurve.ensure3D(Q);

            m = size(Q,1) - 1;
            method = lower(string(method));

            switch method
                case "uniform"
                    u = linspace(0, 1, m+1).';
                    return

                case "chord"
                    d = zeros(m,1);
                    for k = 1:m
                        d(k) = norm(Q(k+1,:) - Q(k,:));
                    end

                case "centripetal"
                    d = zeros(m,1);
                    for k = 1:m
                        d(k) = sqrt(norm(Q(k+1,:) - Q(k,:)));
                    end

                case "arc_length"
                    use_degree = [];
                    if ~isempty(varargin)
                        use_degree = varargin{1};
                    end
                    if isempty(use_degree)
                        use_degree = min(3, m);
                    end
                    use_degree = max(1, min(use_degree, m));

                    u0 = geom.NURBSCurve.parameterizeData(Q, 'chord');
                    C0 = geom.NURBSCurve.globalInterpWithParams(Q, use_degree, u0);

                    us = u0;
                    pts = C0.evaluate(us);
                    s = [0; cumsum(vecnorm(diff(pts,1,1),2,2))];
                    if s(end) < eps
                        u = linspace(0,1,m+1).';
                    else
                        u = s / s(end);
                    end
                    return

                otherwise
                    error('NURBSCurve:parameterizeData', ...
                        'Unknown parameterization method "%s".', method);
            end

            total = sum(d);
            u = zeros(m+1,1);

            if total < eps
                u = linspace(0,1,m+1).';
                return
            end

            for k = 2:m+1
                u(k) = u(k-1) + d(k-1)/total;
            end
            u(end) = 1;
        end

        function C = globalInterpWithParams(Q, p, u)
        % Global interpolation of Q using user-supplied parameters u.
            Q = geom.NURBSCurve.ensure3D(Q);
            u = u(:);

            n = size(Q,1) - 1;
            if numel(u) ~= n+1
                error('NURBSCurve:globalInterpWithParams', ...
                    'u must have one parameter per data point.');
            end
            if p > n
                error('NURBSCurve:globalInterpWithParams', ...
                    'Degree exceeds number of data points minus one.');
            end

            U = geom.NURBSCurve.averagingKnotVector(u, p);

            A = zeros(n+1, n+1);
            for r = 1:n+1
                span = geom.BasisFunctions.FindSpan(n, p, u(r), U);
                N = geom.BasisFunctions.BasisFuns(span, u(r), p, U);
                cols = (span-p):span;
                A(r, cols) = N;
            end

            P = A \ Q;
            C = geom.NURBSCurve(P, p, U, ones(n+1,1));
        end

        function C = globalInterp(Q, p, method)
            if nargin < 3 || isempty(method), method = 'centripetal'; end

            Q = geom.NURBSCurve.ensure3D(Q);
            n = size(Q,1) - 1;

            if p > n
                error('NURBSCurve:globalInterp', ...
                    'Degree exceeds number of data points minus one.');
            end

            u = geom.NURBSCurve.parameterizeData(Q, method, p);
            C = geom.NURBSCurve.globalInterpWithParams(Q, p, u);
        end

        function C = globalLeastSquaresFit(Q, p, nCtrl, method)
        % Global least-squares fit of data Q with nCtrl control points.
            if nargin < 4 || isempty(method), method = 'centripetal'; end

            Q = geom.NURBSCurve.ensure3D(Q);

            m = size(Q,1) - 1;
            n = nCtrl - 1;

            if nCtrl < p + 1
                error('NURBSCurve:globalLeastSquaresFit', ...
                    'Need at least p+1 control points.');
            end
            if m < n
                error('NURBSCurve:globalLeastSquaresFit', ...
                    'Need at least as many data points as control points.');
            end

            u = geom.NURBSCurve.parameterizeData(Q, method, p);
            U = geom.NURBSCurve.makeApproximationKnotVector(u, n, p);

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

        function C = globalWeightedLeastSquaresFit(Q, p, nCtrl, weights, method)
        % Weighted global least-squares fit.
            if nargin < 5 || isempty(method), method = 'centripetal'; end

            Q = geom.NURBSCurve.ensure3D(Q);
            weights = weights(:);

            m = size(Q,1) - 1;
            n = nCtrl - 1;

            if numel(weights) ~= m+1
                error('NURBSCurve:globalWeightedLeastSquaresFit', ...
                    'weights must match the number of data points.');
            end
            if any(weights <= 0)
                error('NURBSCurve:globalWeightedLeastSquaresFit', ...
                    'weights must be strictly positive.');
            end
            if nCtrl < p + 1
                error('NURBSCurve:globalWeightedLeastSquaresFit', ...
                    'Need at least p+1 control points.');
            end
            if m < n
                error('NURBSCurve:globalWeightedLeastSquaresFit', ...
                    'Need at least as many data points as control points.');
            end

            u = geom.NURBSCurve.parameterizeData(Q, method, p);
            U = geom.NURBSCurve.makeApproximationKnotVector(u, n, p);

            Nmat = zeros(m+1, n+1);
            for r = 1:m+1
                span = geom.BasisFunctions.FindSpan(n, p, u(r), U);
                N = geom.BasisFunctions.BasisFuns(span, u(r), p, U);
                cols = (span-p):span;
                Nmat(r, cols) = N;
            end

            Wsqrt = diag(sqrt(weights));
            A = Wsqrt * Nmat;
            B = Wsqrt * Q;

            P = A \ B;
            C = geom.NURBSCurve(P, p, U, ones(n+1,1));
        end

        function C = globalConstrainedLeastSquaresFit(Q, p, nCtrl, fixedIdx, fixedPts, method)
        % Constrained global LSQ fit.
            if nargin < 6 || isempty(method), method = 'centripetal'; end

            Q = geom.NURBSCurve.ensure3D(Q);
            fixedIdx = fixedIdx(:);
            fixedPts = geom.NURBSCurve.ensure3D(fixedPts);

            m = size(Q,1) - 1;
            n = nCtrl - 1;

            if nCtrl < p + 1
                error('NURBSCurve:globalConstrainedLeastSquaresFit', ...
                    'Need at least p+1 control points.');
            end
            if m < n
                error('NURBSCurve:globalConstrainedLeastSquaresFit', ...
                    'Need at least as many data points as control points.');
            end
            if size(fixedPts,1) ~= numel(fixedIdx)
                error('NURBSCurve:globalConstrainedLeastSquaresFit', ...
                    'fixedIdx and fixedPts must have matching lengths.');
            end
            if any(fixedIdx < 1 | fixedIdx > nCtrl)
                error('NURBSCurve:globalConstrainedLeastSquaresFit', ...
                    'fixedIdx out of bounds.');
            end
            if numel(unique(fixedIdx)) ~= numel(fixedIdx)
                error('NURBSCurve:globalConstrainedLeastSquaresFit', ...
                    'fixedIdx must be unique.');
            end

            u = geom.NURBSCurve.parameterizeData(Q, method, p);
            U = geom.NURBSCurve.makeApproximationKnotVector(u, n, p);

            Nmat = zeros(m+1, n+1);
            for r = 1:m+1
                span = geom.BasisFunctions.FindSpan(n, p, u(r), U);
                N = geom.BasisFunctions.BasisFuns(span, u(r), p, U);
                cols = (span-p):span;
                Nmat(r, cols) = N;
            end

            allIdx  = (1:nCtrl).';
            freeIdx = setdiff(allIdx, fixedIdx, 'stable');

            Afree = Nmat(:, freeIdx);
            Afixed = Nmat(:, fixedIdx);

            rhs = Q - Afixed * fixedPts;
            Pfree = Afree \ rhs;

            P = zeros(nCtrl, 3);
            P(fixedIdx, :) = fixedPts;
            P(freeIdx, :) = Pfree;

            C = geom.NURBSCurve(P, p, U, ones(nCtrl,1));
        end

        function C = globalLeastSquaresFitFixedEnds(Q, p, nCtrl, method)
        % LSQ fit with first and last control points fixed to the data endpoints.
            if nargin < 4 || isempty(method), method = 'centripetal'; end

            Q = geom.NURBSCurve.ensure3D(Q);
            fixedIdx = [1; nCtrl];
            fixedPts = [Q(1,:); Q(end,:)];

            C = geom.NURBSCurve.globalConstrainedLeastSquaresFit( ...
                Q, p, nCtrl, fixedIdx, fixedPts, method);
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

        function [ok, Uout, Qw] = removeOneKnotA58(p, U, Pw, u, tol)
        %REMOVEONEKNOTA58 Remove one copy of knot u using Algorithm A5.8 logic.
        %
        % This is a one-copy removal wrapper around the A5.8 local update.
        % It works in homogeneous coordinates.
        
            ok = false;
            U = U(:).';
        
            n = size(Pw,1) - 1;
            m = n + p + 1;
        
            % Find the last occurrence r of knot u, using book-style zero-based r.
            hits = find(abs(U - u) < 1e-12);
            if isempty(hits)
                Uout = U;
                Qw = Pw;
                return;
            end
        
            % Do not remove end knots here.
            if hits(1) <= p+1 || hits(end) >= numel(U)-p
                Uout = U;
                Qw = Pw;
                return;
            end
        
            s = numel(hits);
            r = hits(end) - 1;     % zero-based index of last occurrence
        
            % Algorithm A5.8 local index range for one attempted removal.
            ord = p + 1;
            first = r - p;         % zero-based
            last  = r - s;         % zero-based
        
            if first < 0 || last+1 > n
                Uout = U;
                Qw = Pw;
                return;
            end
        
            temp = zeros(last - first + 2, size(Pw,2));
        
            % MATLAB indices corresponding to book indices:
            % temp(1)            stores Pw(first-1)
            % temp(end)          stores Pw(last+1)
            temp(1,:)   = Pw(first,:);      % book Pw[first-1] -> MATLAB first
            temp(end,:) = Pw(last+2,:);     % book Pw[last+1]  -> MATLAB last+2
        
            i  = first;
            j  = last;
            ii = 1;
            jj = size(temp,1);
        
            remflag = false;
        
            while (j - i) > 0
                % Book:
                % alfi = (u - U[i]) / (U[i+ord] - U[i])
                % alfj = (u - U[j]) / (U[j+ord] - U[j])
                %
                % with zero-based U/P indices converted to MATLAB by +1.
        
                deni = U(i + ord + 1) - U(i + 1);
                denj = U(j + ord + 1) - U(j + 1);
        
                if abs(deni) < eps || abs(denj) < eps
                    Uout = U;
                    Qw = Pw;
                    return;
                end
        
                alfi = (u - U(i + 1)) / deni;
                alfj = (u - U(j + 1)) / denj;
        
                if abs(alfi) < eps || abs(1 - alfj) < eps
                    Uout = U;
                    Qw = Pw;
                    return;
                end
        
                % Book:
                % temp[ii] = (Pw[i] - (1-alfi)*temp[ii-1]) / alfi
                % temp[jj] = (Pw[j] - alfj*temp[jj+1]) / (1-alfj)
                %
                % MATLAB:
                temp(ii+1,:) = (Pw(i+1,:) - (1-alfi)*temp(ii,:)) / alfi;
                temp(jj-1,:) = (Pw(j+1,:) - alfj*temp(jj,:)) / (1-alfj);
        
                i  = i + 1;
                ii = ii + 1;
                j  = j - 1;
                jj = jj - 1;
            end
        
            if j < i
                % Even case: compare the two independently computed middle points.
                if norm(temp(ii,:) - temp(jj,:)) <= tol
                    remflag = true;
                end
            else
                % Odd case: compare reconstructed point to existing middle point.
                den = U(i + ord + 1) - U(i + 1);
                if abs(den) < eps
                    Uout = U;
                    Qw = Pw;
                    return;
                end
        
                alfa = (u - U(i + 1)) / den;
                test = alfa * temp(ii+1,:) + (1-alfa) * temp(ii,:);
        
                if norm(Pw(i+1,:) - test) <= tol
                    remflag = true;
                end
            end
        
            if ~remflag
                Uout = U;
                Qw = Pw;
                return;
            end
        
            % Accepted: build new control polygon.
            Qw = zeros(n, size(Pw,2));  % one fewer control point
        
            % Unaffected before local region.
            if first >= 1
                Qw(1:first,:) = Pw(1:first,:);
            end
        
            % Updated local region.
            for k = 1:(size(temp,1)-2)
                dst = first + k;
                if dst >= 1 && dst <= size(Qw,1)
                    Qw(dst,:) = temp(k+1,:);
                end
            end
        
            % Unaffected after local region.
            srcStart = last + 2;   % MATLAB index in old Pw
            dstStart = srcStart - 1;
            if srcStart <= size(Pw,1)
                Qw(dstStart:end,:) = Pw(srcStart:end,:);
            end
        
            % Remove one copy of u from the knot vector.
            removeIdx = hits(end);
            Uout = U;
            Uout(removeIdx) = [];
        
            if any(diff(Uout) < -1e-14)
                ok = false;
                return;
            end
        
            ok = true;
        end

        function [Qw, Uq] = curveKnotInsertHomogeneous(p, U, Pw, u, r)
        %CURVEKNOTINSERTHOMOGENEOUS Exact homogeneous curve knot insertion.
        % Implements Piegl & Tiller Algorithm A5.1, CurveKnotIns.
        %
        % Inputs:
        %   p  degree
        %   U  knot vector
        %   Pw homogeneous control points [n+1 x dim]
        %   u  knot to insert
        %   r  number of insertions
        %
        % Outputs:
        %   Qw new homogeneous control points
        %   Uq new knot vector
        
            if nargin < 5 || isempty(r), r = 1; end
            r = floor(r);
        
            if r < 0
                error('NURBSCurve:curveKnotInsertHomogeneous', ...
                    'r must be a nonnegative integer.');
            end
        
            U = U(:).';
        
            if r == 0
                Qw = Pw;
                Uq = U;
                return;
            end
        
            n  = size(Pw,1) - 1;
            mp = n + p + 1;
        
            % FindSpan in this codebase returns a 1-based MATLAB span.
            k1 = geom.BasisFunctions.FindSpan(n, p, u, U);
            k  = k1 - 1;   % book zero-based span
        
            s = sum(abs(U - u) < 1e-12);
        
            if s + r > p
                error('NURBSCurve:curveKnotInsertHomogeneous', ...
                    'Cannot insert knot: multiplicity s+r exceeds degree p.');
            end
        
            nq = n + r;
        
            Uq = zeros(1, mp + r + 1);
            Qw = zeros(nq + 1, size(Pw,2));
        
            % Load new knot vector.
            for i = 0:k
                Uq(i+1) = U(i+1);
            end
        
            for i = 1:r
                Uq(k+i+1) = u;
            end
        
            for i = k+1:mp
                Uq(i+r+1) = U(i+1);
            end
        
            % Save unaltered control points.
            for i = 0:k-p
                Qw(i+1,:) = Pw(i+1,:);
            end
        
            for i = k-s:n
                Qw(i+r+1,:) = Pw(i+1,:);
            end
        
            % Local affected control points.
            Rw = zeros(p-s+1, size(Pw,2));
            for i = 0:p-s
                Rw(i+1,:) = Pw(k-p+i+1,:);
            end
        
            % Insert knot r times.
            L = 0;
            for j = 1:r
                L = k - p + j;
        
                for i = 0:p-j-s
                    denom = U(i+k+2) - U(L+i+1);
                    if abs(denom) < 1e-14
                        alpha = 0.0;
                    else
                        alpha = (u - U(L+i+1)) / denom;
                    end
        
                    Rw(i+1,:) = alpha * Rw(i+2,:) + ...
                                (1.0 - alpha) * Rw(i+1,:);
                end
        
                Qw(L+1,:) = Rw(1,:);
                Qw(k+r-j-s+1,:) = Rw(p-j-s+1,:);
            end
        
            % Load remaining control points.
            for i = L+1:k-s-1
                Qw(i+1,:) = Rw(i-L+1,:);
            end
        end

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

        function [Qw, maxErr] = reduceBezierHomogeneousExact(Pw, tol)
        %REDUCEBEZIERHOMOGENEOUSEXACT Exact/tolerance Bezier degree reduction.
        %
        % Reduces a Bezier control polygon from degree p to p-1 using the
        % degree-elevation identities in reverse.
        %
        % If the control polygon is not reducible within tol, this errors.
        
            if nargin < 2 || isempty(tol), tol = 1e-10; end
        
            p = size(Pw,1) - 1;
            dim = size(Pw,2);
        
            if p < 1
                error('NURBSCurve:reduceBezierHomogeneousExact', ...
                    'Bezier degree must be at least 1.');
            end
        
            q = p - 1;
            Qw = zeros(q+1, dim);
        
            Qw(1,:)   = Pw(1,:);
            Qw(end,:) = Pw(end,:);
        
            if p == 1
                maxErr = 0;
                return;
            end
        
            if mod(p,2) == 1
                % Odd p: both recurrences compute the same middle reduced point.
                r = (p - 1) / 2;
        
                Qleft = zeros(q+1, dim);
                Qright = zeros(q+1, dim);
                Qleft(1,:) = Pw(1,:);
                Qright(end,:) = Pw(end,:);
        
                % Left recurrence.
                for i = 1:r
                    alpha = i / p;
                    Qleft(i+1,:) = (Pw(i+1,:) - alpha * Qleft(i,:)) / (1 - alpha);
                end
        
                % Right recurrence.
                for i = p-1:-1:r+1
                    alpha = i / p;
                    Qright(i,:) = (Pw(i+1,:) - (1 - alpha) * Qright(i+1,:)) / alpha;
                end
        
                midErr = norm(Qleft(r+1,:) - Qright(r+1,:));
                maxErr = midErr;
        
                if midErr > tol
                    error('NURBSCurve:reduceBezierHomogeneousExact', ...
                        'Bezier segment is not degree-reducible: %.3e > %.3e.', ...
                        midErr, tol);
                end
        
                Qw(1:r+1,:) = Qleft(1:r+1,:);
                Qw(r+2:end,:) = Qright(r+2:end,:);
        
            else
                % Even p: left and right recurrences meet through a middle
                % original Bezier point consistency check.
                r = p / 2;
        
                % Left recurrence computes Q_1 ... Q_{r-1}.
                for i = 1:r-1
                    alpha = i / p;
                    Qw(i+1,:) = (Pw(i+1,:) - alpha * Qw(i,:)) / (1 - alpha);
                end
        
                % Right recurrence computes Q_{p-2} ... Q_r.
                for i = p-1:-1:r+1
                    alpha = i / p;
                    Qw(i,:) = (Pw(i+1,:) - (1 - alpha) * Qw(i+1,:)) / alpha;
                end
        
                alpha = r / p;
                Pmid = alpha * Qw(r,:) + (1 - alpha) * Qw(r+1,:);
                midErr = norm(Pw(r+1,:) - Pmid);
                maxErr = midErr;
        
                if midErr > tol
                    error('NURBSCurve:reduceBezierHomogeneousExact', ...
                        'Bezier segment is not degree-reducible: %.3e > %.3e.', ...
                        midErr, tol);
                end
            end
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

        function U = makeApproximationKnotVector(u, n, p)
            u = u(:);
            m = numel(u) - 1;

            U = zeros(1, n + p + 2);
            U(1:p+1) = 0;
            U(end-p:end) = 1;

            if (n - p) >= 1
                d = (m + 1) / (n - p + 1);
                for j = 1:(n-p)
                    jd = j * d;
                    i = floor(jd);
                    alpha = jd - i;

                    i = max(1, min(i, m));
                    U(j+p+1) = (1-alpha) * u(i) + alpha * u(i+1);
                end
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
        %INSERTKNOTONCE Compatibility wrapper.
            C2 = obj.insertKnot(u, 1);
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