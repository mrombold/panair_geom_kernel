classdef NURBSSurface < handle
% NURBSSURFACE  Tensor-product rational B-spline surface kernel.
%
% Features
%   - exact rational evaluation
%   - exact rational derivative tensor
%   - analytic first/second partials
%   - normals and curvatures
%   - closest-point projection
%   - refinement, exact splitting, iso-curves, subpatch extraction
%   - Bezier patch decomposition
%   - degree elevation / reduction in u and v
%   - knot removal in u and v
%   - rectangular-net interpolation / least-squares fitting
%   - ruled and bilinear constructors
%   - Coons patch
%   - Gordon surface
%   - loft / skinning with compatible knot merging
%   - surface of revolution
%   - swept surface
%   - trimmed-surface support in UV parameter space
%   - network compatibility helpers

    properties
        P
        W
        U
        V
        p
        q

        % Trim support:
        % trimOuterLoops / trimInnerLoops are cell arrays of loops.
        % Each loop is one of:
        %   - Nx2 numeric UV polyline (closed or open, auto-closed)
        %   - geom.NURBSCurve living in UV plane (z ignored)
        %
        % Semantics:
        %   - If trimOuterLoops is empty, the full param domain is active.
        %   - If trimOuterLoops is nonempty, a point must lie inside at least
        %     one outer loop to be kept.
        %   - If trimInnerLoops is nonempty, a point inside any inner loop
        %     is excluded.
        trimOuterLoops = {}
        trimInnerLoops = {}
        trimTolerance = 1e-10
    end

    properties (Dependent)
        n
        m
        domainU
        domainV
    end

    methods
        function obj = NURBSSurface(P, p, q, U, V, W)
            if nargin < 3
                error('NURBSSurface:Constructor', 'Requires at least P, p, q.');
            end

            if isstruct(P)
                P = cat(3, P.x, P.y, P.z);
            end

            sz = size(P);
            if numel(sz) ~= 3 || sz(3) ~= 3
                error('NURBSSurface:Constructor', 'P must be an (n+1)x(m+1)x3 array.');
            end

            if ~isscalar(p) || p < 0 || p ~= floor(p)
                error('NURBSSurface:Constructor', 'p must be a nonnegative integer.');
            end
            if ~isscalar(q) || q < 0 || q ~= floor(q)
                error('NURBSSurface:Constructor', 'q must be a nonnegative integer.');
            end

            n = sz(1) - 1;
            m = sz(2) - 1;
            if p > n
                error('NURBSSurface:Constructor', 'Degree p=%d exceeds n=%d.', p, n);
            end
            if q > m
                error('NURBSSurface:Constructor', 'Degree q=%d exceeds m=%d.', q, m);
            end

            if nargin < 6 || isempty(W)
                W = ones(n+1, m+1);
            else
                if ~isequal(size(W), [n+1, m+1])
                    error('NURBSSurface:Constructor', 'W must be size (%d x %d).', n+1, m+1);
                end
                if any(W(:) <= 0)
                    error('NURBSSurface:Constructor', 'Weights must be strictly positive.');
                end
            end

            if nargin < 4 || isempty(U)
                U = geom.BasisFunctions.MakeUniformKnotVector(n, p);
            else
                U = U(:).';
            end
            if nargin < 5 || isempty(V)
                V = geom.BasisFunctions.MakeUniformKnotVector(m, q);
            else
                V = V(:).';
            end

            expectedU = n + p + 2;
            expectedV = m + q + 2;
            if numel(U) ~= expectedU
                error('NURBSSurface:Constructor', ...
                    'U length must be n+p+2 = %d, got %d.', expectedU, numel(U));
            end
            if numel(V) ~= expectedV
                error('NURBSSurface:Constructor', ...
                    'V length must be m+q+2 = %d, got %d.', expectedV, numel(V));
            end
            if any(diff(U) < 0)
                error('NURBSSurface:Constructor', 'U must be nondecreasing.');
            end
            if any(diff(V) < 0)
                error('NURBSSurface:Constructor', 'V must be nondecreasing.');
            end

            obj.P = P;
            obj.W = W;
            obj.U = U;
            obj.V = V;
            obj.p = p;
            obj.q = q;

            obj.validate();
        end

        function v = get.n(obj)
            v = size(obj.P, 1) - 1;
        end

        function v = get.m(obj)
            v = size(obj.P, 2) - 1;
        end

        function v = get.domainU(obj)
            v = [obj.U(obj.p+1), obj.U(end-obj.p)];
        end

        function v = get.domainV(obj)
            v = [obj.V(obj.q+1), obj.V(end-obj.q)];
        end
    end

    methods
        function tf = validate(obj)
            if size(obj.P, 1) ~= size(obj.W, 1) || size(obj.P, 2) ~= size(obj.W, 2)
                error('NURBSSurface:Validate', 'P and W sizes are inconsistent.');
            end
            if numel(obj.U) ~= obj.n + obj.p + 2
                error('NURBSSurface:Validate', 'U length is inconsistent with n and p.');
            end
            if numel(obj.V) ~= obj.m + obj.q + 2
                error('NURBSSurface:Validate', 'V length is inconsistent with m and q.');
            end
            if any(diff(obj.U) < 0)
                error('NURBSSurface:Validate', 'U must be nondecreasing.');
            end
            if any(diff(obj.V) < 0)
                error('NURBSSurface:Validate', 'V must be nondecreasing.');
            end
            if any(obj.W(:) <= 0)
                error('NURBSSurface:Validate', 'Weights must be strictly positive.');
            end
            tf = true;
        end

        function C = evaluate(obj, u, v)
            u = u(:);
            v = v(:);
            if isscalar(u), u = repmat(u, numel(v), 1); end
            if isscalar(v), v = repmat(v, numel(u), 1); end
            assert(numel(u) == numel(v), ...
                'NURBSSurface:evaluate', 'u and v must match in length.');

            np = numel(u);
            C = zeros(np, 3);

            for k = 1:np
                uk = obj.clampU(u(k));
                vk = obj.clampV(v(k));

                uspan = geom.BasisFunctions.FindSpan(obj.n, obj.p, uk, obj.U);
                vspan = geom.BasisFunctions.FindSpan(obj.m, obj.q, vk, obj.V);

                Nu = geom.BasisFunctions.BasisFuns(uspan, uk, obj.p, obj.U);
                Nv = geom.BasisFunctions.BasisFuns(vspan, vk, obj.q, obj.V);

                Sw = zeros(1, 4);
                for l = 0:obj.q
                    temp = zeros(1, 4);
                    jj = vspan - obj.q + l;
                    for i = 0:obj.p
                        ii = uspan - obj.p + i;
                        wij = obj.W(ii, jj);
                        Pij = reshape(obj.P(ii, jj, :), 1, 3);
                        temp = temp + Nu(i+1) * [wij * Pij, wij];
                    end
                    Sw = Sw + Nv(l+1) * temp;
                end

                C(k,:) = Sw(1:3) / Sw(4);
            end
        end

        function C = evaluateTrimmed(obj, u, v)
        % Evaluate only if (u,v) lies inside the active trimmed region.
            if ~obj.isInsideTrim(u, v)
                error('NURBSSurface:evaluateTrimmed', ...
                    'Requested parameter point lies outside trimmed region.');
            end
            C = obj.evaluate(u, v);
        end

        function tf = isTrimmed(obj)
            tf = ~isempty(obj.trimOuterLoops) || ~isempty(obj.trimInnerLoops);
        end

        function obj2 = clearTrims(obj)
            obj2 = obj.copySurfaceOnly();
            obj2.trimOuterLoops = {};
            obj2.trimInnerLoops = {};
            obj2.trimTolerance = obj.trimTolerance;
        end

        function obj2 = setTrims(obj, outerLoops, innerLoops)
            if nargin < 2 || isempty(outerLoops), outerLoops = {}; end
            if nargin < 3 || isempty(innerLoops), innerLoops = {}; end

            obj2 = obj.copySurfaceOnly();
            obj2.trimOuterLoops = geom.NURBSSurface.normalizeTrimLoopSet(outerLoops);
            obj2.trimInnerLoops = geom.NURBSSurface.normalizeTrimLoopSet(innerLoops);
            obj2.trimTolerance = obj.trimTolerance;
        end

        function obj2 = addOuterTrimLoop(obj, loop)
            obj2 = obj.copySurfaceOnly();
            obj2.trimOuterLoops = [obj.trimOuterLoops, ...
                geom.NURBSSurface.normalizeTrimLoopSet({loop})];
            obj2.trimInnerLoops = obj.trimInnerLoops;
            obj2.trimTolerance = obj.trimTolerance;
        end

        function obj2 = addInnerTrimLoop(obj, loop)
            obj2 = obj.copySurfaceOnly();
            obj2.trimOuterLoops = obj.trimOuterLoops;
            obj2.trimInnerLoops = [obj.trimInnerLoops, ...
                geom.NURBSSurface.normalizeTrimLoopSet({loop})];
            obj2.trimTolerance = obj.trimTolerance;
        end

        function tf = isInsideTrim(obj, u, v)
        % True if (u,v) lies in the active kept region of the trimmed surface.
            u = obj.clampU(u);
            v = obj.clampV(v);

            tol = obj.trimTolerance;
            if ~obj.isTrimmed()
                tf = true;
                return;
            end

            % Outer logic
            if isempty(obj.trimOuterLoops)
                insideOuter = true;
            else
                insideOuter = false;
                for k = 1:numel(obj.trimOuterLoops)
                    loop = obj.trimOuterLoops{k};
                    if geom.NURBSSurface.pointInPolygon2D(loop, [u v], tol)
                        insideOuter = true;
                        break;
                    end
                end
            end

            if ~insideOuter
                tf = false;
                return;
            end

            % Inner holes
            insideInner = false;
            for k = 1:numel(obj.trimInnerLoops)
                loop = obj.trimInnerLoops{k};
                if geom.NURBSSurface.pointInPolygon2D(loop, [u v], tol)
                    insideInner = true;
                    break;
                end
            end

            tf = ~insideInner;
        end

        function mask = trimMask(obj, uGrid, vGrid)
            if ~isequal(size(uGrid), size(vGrid))
                error('NURBSSurface:trimMask', 'uGrid and vGrid must have same size.');
            end
            mask = false(size(uGrid));
            for i = 1:size(uGrid,1)
                for j = 1:size(uGrid,2)
                    mask(i,j) = obj.isInsideTrim(uGrid(i,j), vGrid(i,j));
                end
            end
        end

        function UV = sampleTrimLoops(obj, nPerLoop)
            if nargin < 2 || isempty(nPerLoop), nPerLoop = 120; end
            UV.outer = cell(size(obj.trimOuterLoops));
            UV.inner = cell(size(obj.trimInnerLoops));

            for k = 1:numel(obj.trimOuterLoops)
                UV.outer{k} = geom.NURBSSurface.resampleTrimLoop(obj.trimOuterLoops{k}, nPerLoop);
            end
            for k = 1:numel(obj.trimInnerLoops)
                UV.inner{k} = geom.NURBSSurface.resampleTrimLoop(obj.trimInnerLoops{k}, nPerLoop);
            end
        end

        function XYZ = sampleTrimLoopsXYZ(obj, nPerLoop)
            if nargin < 2 || isempty(nPerLoop), nPerLoop = 120; end
            UV = obj.sampleTrimLoops(nPerLoop);
            XYZ.outer = cell(size(UV.outer));
            XYZ.inner = cell(size(UV.inner));

            for k = 1:numel(UV.outer)
                uv = UV.outer{k};
                XYZ.outer{k} = obj.evaluate(uv(:,1), uv(:,2));
            end
            for k = 1:numel(UV.inner)
                uv = UV.inner{k};
                XYZ.inner{k} = obj.evaluate(uv(:,1), uv(:,2));
            end
        end

        function SKL = derivatives(obj, u, v, d)
            if nargin < 4 || isempty(d)
                d = 1;
            end
            if ~isscalar(u) || ~isscalar(v)
                error('NURBSSurface:derivatives', 'Use scalar u,v for derivatives().');
            end

            u = obj.clampU(u);
            v = obj.clampV(v);
            d = max(0, floor(d));

            du = min(d, obj.p);
            dv = min(d, obj.q);

            uspan = geom.BasisFunctions.FindSpan(obj.n, obj.p, u, obj.U);
            vspan = geom.BasisFunctions.FindSpan(obj.m, obj.q, v, obj.V);

            Nu = geom.BasisFunctions.DersBasisFuns(uspan, u, obj.p, du, obj.U);
            Nv = geom.BasisFunctions.DersBasisFuns(vspan, v, obj.q, dv, obj.V);

            Aders = cell(d+1, d+1);
            wders = zeros(d+1, d+1);
            for k = 0:d
                for l = 0:d
                    Aders{k+1,l+1} = zeros(1,3);
                end
            end

            for k = 0:du
                for l = 0:dv
                    Akl = zeros(1,3);
                    wkl = 0.0;
                    for i = 0:obj.p
                        ii = uspan - obj.p + i;
                        for j = 0:obj.q
                            jj = vspan - obj.q + j;
                            B = Nu(k+1, i+1) * Nv(l+1, j+1);
                            wij = obj.W(ii, jj);
                            Pij = reshape(obj.P(ii, jj, :), 1, 3);
                            Akl = Akl + B * wij * Pij;
                            wkl = wkl + B * wij;
                        end
                    end
                    Aders{k+1,l+1} = Akl;
                    wders(k+1,l+1) = wkl;
                end
            end

            SKL = cell(d+1, d+1);
            for k = 0:d
                for l = 0:d
                    SKL{k+1,l+1} = zeros(1,3);
                end
            end

            w00 = wders(1,1);
            if abs(w00) < eps
                error('NURBSSurface:ZeroWeight', ...
                    'Surface weight function vanished at (u,v).');
            end

            for k = 0:d
                for l = 0:(d-k)
                    vkl = Aders{k+1,l+1};

                    for j = 1:l
                        vkl = vkl - nchoosek(l,j) * wders(1,j+1) * SKL{k+1,l-j+1};
                    end

                    for i = 1:k
                        vkl = vkl - nchoosek(k,i) * wders(i+1,1) * SKL{k-i+1,l+1};
                        for j = 1:l
                            vkl = vkl - nchoosek(k,i) * nchoosek(l,j) * ...
                                wders(i+1,j+1) * SKL{k-i+1,l-j+1};
                        end
                    end

                    SKL{k+1,l+1} = vkl / w00;
                end
            end
        end

        function [Su, Sv] = partialDerivatives(obj, u, v)
            SKL = obj.derivatives(u, v, 1);
            Su = SKL{2,1};
            Sv = SKL{1,2};
        end

        function [S, Su, Sv] = firstPartials(obj, u, v)
            SKL = obj.derivatives(u, v, 1);
            S  = SKL{1,1};
            Su = SKL{2,1};
            Sv = SKL{1,2};
        end

        function [Suu, Suv, Svv] = secondPartials(obj, u, v)
            SKL = obj.derivatives(u, v, 2);
            Suu = SKL{3,1};
            Suv = SKL{2,2};
            Svv = SKL{1,3};
        end

        function [Su, Sv, Suu, Suv, Svv] = firstSecond(obj, u, v)
            SKL = obj.derivatives(u, v, 2);
            Su  = SKL{2,1};
            Sv  = SKL{1,2};
            Suu = SKL{3,1};
            Suv = SKL{2,2};
            Svv = SKL{1,3};
        end

        function N = normal(obj, u, v)
            [Su, Sv] = obj.partialDerivatives(u, v);
            N = cross(Su, Sv);
            nrm = norm(N);
            if nrm > eps
                N = N / nrm;
            else
                N = [0 0 0];
            end
        end

        function [K, H, Nhat] = curvatures(obj, u, v)
            [Su, Sv, Suu, Suv, Svv] = obj.firstSecond(u, v);

            N = cross(Su, Sv);
            nrm = norm(N);
            if nrm < eps
                Nhat = [0 0 0];
                K = NaN;
                H = NaN;
                return;
            end
            Nhat = N / nrm;

            E = dot(Su, Su);
            F = dot(Su, Sv);
            G = dot(Sv, Sv);
            L = dot(Suu, Nhat);
            M = dot(Suv, Nhat);
            N2 = dot(Svv, Nhat);

            den = E * G - F * F;
            if abs(den) < eps
                K = NaN;
                H = NaN;
                return;
            end

            K = (L * N2 - M * M) / den;
            H = (E * N2 - 2 * F * M + G * L) / (2 * den);
        end

        function [u_c, v_c, pt_c, d_c] = closestPoint(obj, P, u0, v0)
            P = P(:).';

            if nargin < 3 || isempty(u0) || nargin < 4 || isempty(v0)
                [u0, v0] = obj.coarseSearchPoint(P, 12, 12);
            end

            u = obj.clampU(u0);
            v = obj.clampV(v0);

            for iter = 1:50
                S = obj.evaluate(u, v);
                [Su, Sv, Suu, Suv, Svv] = obj.firstSecond(u, v);

                r = S - P;
                F = [dot(r, Su); dot(r, Sv)];
                J = [dot(Su, Su) + dot(r, Suu), dot(Su, Sv) + dot(r, Suv); ...
                     dot(Su, Sv) + dot(r, Suv), dot(Sv, Sv) + dot(r, Svv)];

                if rcond(J) < eps
                    break;
                end

                delta = -J \ F;
                u_new = obj.clampU(u + delta(1));
                v_new = obj.clampV(v + delta(2));

                if norm([u_new-u, v_new-v]) < 1e-12
                    u = u_new;
                    v = v_new;
                    break;
                end

                u = u_new;
                v = v_new;
            end

            pt_c = obj.evaluate(u, v);
            d_c = norm(pt_c - P);
            u_c = u;
            v_c = v;
        end

        function [u_all, v_all, pt_all, d_all] = closestPointBatch(obj, pts, nu_c, nv_c)
            if nargin < 3 || isempty(nu_c), nu_c = 16; end
            if nargin < 4 || isempty(nv_c), nv_c = 16; end

            N = size(pts, 1);
            u_all = zeros(N, 1);
            v_all = zeros(N, 1);
            pt_all = zeros(N, 3);
            d_all = zeros(N, 1);

            [u_g, v_g, pts_g] = obj.isoGrid(nu_c, nv_c);
            pts_flat = reshape(pts_g, [], 3);

            [UG, VG] = meshgrid(u_g, v_g);
            u_flat = UG.';
            v_flat = VG.';
            u_flat = u_flat(:);
            v_flat = v_flat(:);

            for k = 1:N
                P = pts(k, :);
                d2 = sum((pts_flat - P).^2, 2);
                [~, idx] = min(d2);
                [u_all(k), v_all(k), pt_all(k,:), d_all(k)] = ...
                    obj.closestPoint(P, u_flat(idx), v_flat(idx));
            end
        end

        function [u_inv, v_inv, pt_inv, err] = invertPoint(obj, P, u0, v0)
            if nargin < 3, u0 = []; end
            if nargin < 4, v0 = []; end
            [u_inv, v_inv, pt_inv, err] = obj.closestPoint(P, u0, v0);
        end
    end

    methods
        function [u_iso, v_iso, pts] = isoGrid(obj, nu, nv)
            u_iso = linspace(obj.domainU(1), obj.domainU(2), nu);
            v_iso = linspace(obj.domainV(1), obj.domainV(2), nv);
            pts = zeros(nu, nv, 3);

            for i = 1:nu
                for j = 1:nv
                    pts(i, j, :) = obj.evaluate(u_iso(i), v_iso(j));
                end
            end
        end

        function mesh = isoMesh(obj, nu, nv, varargin)
            pa = inputParser;
            addParameter(pa, 'SpacingU', 'uniform');
            addParameter(pa, 'SpacingV', 'uniform');
            addParameter(pa, 'RespectTrim', true);
            parse(pa, varargin{:});
            opts = pa.Results;

            u_iso = geom.NURBSSurface.makeSpacing(obj.domainU(1), obj.domainU(2), nu, opts.SpacingU);
            v_iso = geom.NURBSSurface.makeSpacing(obj.domainV(1), obj.domainV(2), nv, opts.SpacingV);

            pts = zeros(nu, nv, 3);
            for i = 1:nu
                for j = 1:nv
                    pts(i, j, :) = obj.evaluate(u_iso(i), v_iso(j));
                end
            end

            mesh.X = squeeze(pts(:,:,1));
            mesh.Y = squeeze(pts(:,:,2));
            mesh.Z = squeeze(pts(:,:,3));
            mesh.nu = nu;
            mesh.nv = nv;

            [mesh.u, mesh.v] = meshgrid(u_iso, v_iso);
            mesh.u = mesh.u.';
            mesh.v = mesh.v.';

            if opts.RespectTrim && obj.isTrimmed()
                mesh.trimMask = obj.trimMask(mesh.u, mesh.v);
                X = mesh.X; Y = mesh.Y; Z = mesh.Z;
                X(~mesh.trimMask) = NaN;
                Y(~mesh.trimMask) = NaN;
                Z(~mesh.trimMask) = NaN;
                mesh.X = X;
                mesh.Y = Y;
                mesh.Z = Z;
            else
                mesh.trimMask = true(size(mesh.u));
            end

            Normals = zeros(nu, nv, 3);
            for i = 1:nu
                for j = 1:nv
                    if mesh.trimMask(i,j)
                        Normals(i, j, :) = obj.normal(u_iso(i), v_iso(j));
                    else
                        Normals(i, j, :) = [NaN NaN NaN];
                    end
                end
            end
            mesh.normals = Normals;

            nquads = (nu - 1) * (nv - 1);
            conn = zeros(nquads, 4);
            qid = 1;
            for i = 1:nu-1
                for j = 1:nv-1
                    n1 = (i-1) * nv + j;
                    n2 = i * nv + j;
                    n3 = i * nv + j + 1;
                    n4 = (i-1) * nv + j + 1;
                    conn(qid, :) = [n1, n2, n3, n4];
                    qid = qid + 1;
                end
            end
            mesh.connectivity = conn;
        end

        function S2 = refine(obj, Xu, Xv)
            if nargin < 2, Xu = []; end
            if nargin < 3, Xv = []; end

            P2 = obj.P;
            W2 = obj.W;
            U2 = obj.U;
            V2 = obj.V;

            if ~isempty(Xu)
                Xu = sort(Xu(:).');
                p = obj.p;
                U_orig = U2;

                for j = 1:size(P2, 2)
                    Prow = squeeze(P2(:, j, :));
                    Wrow = W2(:, j);

                    Ctmp = geom.NURBSCurve(Prow, p, U_orig, Wrow);
                    Cref = Ctmp.refine(Xu);

                    if j == 1
                        n_new = size(Cref.P, 1);
                        P2_new = zeros(n_new, size(P2, 2), 3);
                        W2_new = zeros(n_new, size(P2, 2));
                        U2 = Cref.U;
                    end

                    P2_new(:, j, :) = Cref.P;
                    W2_new(:, j) = Cref.W;
                end

                P2 = P2_new;
                W2 = W2_new;
            end

            if ~isempty(Xv)
                Xv = sort(Xv(:).');
                q = obj.q;
                V_orig = V2;

                for i = 1:size(P2, 1)
                    Pcol = squeeze(P2(i, :, :));
                    Wcol = W2(i, :).';

                    Ctmp = geom.NURBSCurve(Pcol, q, V_orig, Wcol);
                    Cref = Ctmp.refine(Xv);

                    if i == 1
                        m_new = size(Cref.P, 1);
                        P2_new = zeros(size(P2, 1), m_new, 3);
                        W2_new = zeros(size(P2, 1), m_new);
                        V2 = Cref.U;
                    end

                    P2_new(i, :, :) = Cref.P;
                    W2_new(i, :) = Cref.W.';
                end

                P2 = P2_new;
                W2 = W2_new;
            end

            S2 = geom.NURBSSurface(P2, obj.p, obj.q, U2, V2, W2);
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance = obj.trimTolerance;
        end

        function S2 = flipNormals(obj)
            P2 = flipud(obj.P);
            W2 = flipud(obj.W);
            U2 = (obj.U(1) + obj.U(end)) - fliplr(obj.U);
            S2 = geom.NURBSSurface(P2, obj.p, obj.q, U2, obj.V, W2);
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance = obj.trimTolerance;
        end

        function [S_lo, S_hi] = splitU(obj, u0)
            tol = 1e-12;
            u0 = obj.clampU(u0);

            if u0 <= obj.domainU(1) + tol || u0 >= obj.domainU(2) - tol
                error('NURBSSurface:splitU', 'u0 must lie strictly inside the active U domain.');
            end

            s = sum(abs(obj.U - u0) < tol);
            if s > obj.p
                error('NURBSSurface:splitU', 'Split knot multiplicity exceeds degree.');
            end

            reps = obj.p - s;
            if reps > 0
                S_ref = obj.refine(repmat(u0, 1, reps), []);
            else
                S_ref = obj;
            end

            U2 = S_ref.U;
            first = find(abs(U2 - u0) < tol, 1, 'first');
            last  = find(abs(U2 - u0) < tol, 1, 'last');

            if isempty(first) || isempty(last)
                error('NURBSSurface:splitU', 'Failed to create split knot.');
            end

            cp_split = last - obj.p;

            P_lo = S_ref.P(1:cp_split, :, :);
            W_lo = S_ref.W(1:cp_split, :);
            U_lo = [U2(1:last), u0];
            U_lo = geom.NURBSSurface.normalizeKnotVector(U_lo);

            P_hi = S_ref.P(cp_split:end, :, :);
            W_hi = S_ref.W(cp_split:end, :);
            U_hi = [u0, U2(first:end)];
            U_hi = geom.NURBSSurface.normalizeKnotVector(U_hi);

            S_lo = geom.NURBSSurface(P_lo, obj.p, obj.q, U_lo, obj.V, W_lo);
            S_hi = geom.NURBSSurface(P_hi, obj.p, obj.q, U_hi, obj.V, W_hi);

            if obj.isTrimmed()
                [outerLo, outerHi] = geom.NURBSSurface.splitTrimLoopsU(obj.trimOuterLoops, u0);
                [innerLo, innerHi] = geom.NURBSSurface.splitTrimLoopsU(obj.trimInnerLoops, u0);

                S_lo.trimOuterLoops = outerLo;
                S_lo.trimInnerLoops = innerLo;
                S_hi.trimOuterLoops = outerHi;
                S_hi.trimInnerLoops = innerHi;
                S_lo.trimTolerance = obj.trimTolerance;
                S_hi.trimTolerance = obj.trimTolerance;
            end
        end

        function [S_lo, S_hi] = splitV(obj, v0)
            St = obj.swapUV();
            [A, B] = St.splitU(v0);
            S_lo = A.swapUV();
            S_hi = B.swapUV();

            S_lo.trimOuterLoops = obj.trimOuterLoops;
            S_lo.trimInnerLoops = obj.trimInnerLoops;
            S_hi.trimOuterLoops = obj.trimOuterLoops;
            S_hi.trimInnerLoops = obj.trimInnerLoops;
            S_lo.trimTolerance = obj.trimTolerance;
            S_hi.trimTolerance = obj.trimTolerance;
        end

        function C = isoCurveU(obj, v0)
            v0 = obj.clampV(v0);
            vspan = geom.BasisFunctions.FindSpan(obj.m, obj.q, v0, obj.V);
            Nv = geom.BasisFunctions.BasisFuns(vspan, v0, obj.q, obj.V);

            n = obj.n;
            Pw = zeros(n + 1, 4);

            for i = 0:n
                for l = 0:obj.q
                    jj = vspan - obj.q + l;
                    wij = obj.W(i+1, jj);
                    Pij = reshape(obj.P(i+1, jj, :), 1, 3);
                    Pw(i+1, :) = Pw(i+1, :) + Nv(l+1) * [wij * Pij, wij];
                end
            end

            W2 = Pw(:,4);
            P2 = Pw(:,1:3) ./ W2;
            C = geom.NURBSCurve(P2, obj.p, obj.U, W2);
        end

        function C = isoCurveV(obj, u0)
            u0 = obj.clampU(u0);
            uspan = geom.BasisFunctions.FindSpan(obj.n, obj.p, u0, obj.U);
            Nu = geom.BasisFunctions.BasisFuns(uspan, u0, obj.p, obj.U);

            m = obj.m;
            Pw = zeros(m + 1, 4);

            for j = 0:m
                for l = 0:obj.p
                    ii = uspan - obj.p + l;
                    wij = obj.W(ii, j+1);
                    Pij = reshape(obj.P(ii, j+1, :), 1, 3);
                    Pw(j+1, :) = Pw(j+1, :) + Nu(l+1) * [wij * Pij, wij];
                end
            end

            W2 = Pw(:,4);
            P2 = Pw(:,1:3) ./ W2;
            C = geom.NURBSCurve(P2, obj.q, obj.V, W2);
        end

        function Ssub = extractUV(obj, urange, vrange)
            ua = urange(1);
            ub = urange(2);
            va = vrange(1);
            vb = vrange(2);

            if ~(ua < ub && va < vb)
                error('NURBSSurface:extractUV', 'Require ua < ub and va < vb.');
            end

            ua = obj.clampU(ua);
            ub = obj.clampU(ub);
            va = obj.clampV(va);
            vb = obj.clampV(vb);

            Swork = obj;
            tol = 1e-12;

            if ua > Swork.domainU(1) + tol
                [~, Swork] = Swork.splitU(ua);
            end

            if ub < obj.domainU(2) - tol
                ub_local = (ub - ua) / (obj.domainU(2) - ua);
                ub_local = max(0, min(1, ub_local));
                [Swork, ~] = Swork.splitU(ub_local);
            end

            if va > Swork.domainV(1) + tol
                [~, Swork] = Swork.splitV(va);
            end

            if vb < obj.domainV(2) - tol
                vb_local = (vb - va) / (obj.domainV(2) - va);
                vb_local = max(0, min(1, vb_local));
                [Swork, ~] = Swork.splitV(vb_local);
            end

            Ssub = Swork;
        end

        function B = boundaryCurves(obj)
            B.u0 = obj.isoCurveV(obj.domainU(1));
            B.u1 = obj.isoCurveV(obj.domainU(2));
            B.v0 = obj.isoCurveU(obj.domainV(1));
            B.v1 = obj.isoCurveU(obj.domainV(2));
        end

        function C = extractInnerSpanCurve(obj)
            u_pick = obj.domainU(1) + 0.05 * (obj.domainU(2) - obj.domainU(1));
            C = obj.isoCurveV(u_pick);
        end

        function C = extractInnerCircCurve(obj)
            v_pick = obj.domainV(1) + 0.05 * (obj.domainV(2) - obj.domainV(1));
            C = obj.isoCurveU(v_pick);
        end

        function St = swapUV(obj)
            P2 = permute(obj.P, [2 1 3]);
            W2 = obj.W.';
            St = geom.NURBSSurface(P2, obj.q, obj.p, obj.V, obj.U, W2);
            St.trimOuterLoops = obj.trimOuterLoops;
            St.trimInnerLoops = obj.trimInnerLoops;
            St.trimTolerance = obj.trimTolerance;
        end

        function patches = decomposeBezier(obj)
            Su = obj;

            Ui = unique(Su.U(Su.p+2:end-Su.p-1));
            for k = 1:numel(Ui)
                s = sum(abs(Su.U - Ui(k)) < 1e-12);
                if s < Su.p
                    Su = Su.refine(repmat(Ui(k), 1, Su.p - s), []);
                end
            end

            Vi = unique(Su.V(Su.q+2:end-Su.q-1));
            for k = 1:numel(Vi)
                s = sum(abs(Su.V - Vi(k)) < 1e-12);
                if s < Su.q
                    Su = Su.refine([], repmat(Vi(k), 1, Su.q - s));
                end
            end

            uBreaks = unique(Su.U);
            vBreaks = unique(Su.V);

            nuSeg = numel(uBreaks) - 1;
            nvSeg = numel(vBreaks) - 1;

            patches = cell(nuSeg, nvSeg);

            for iu = 1:nuSeg
                for jv = 1:nvSeg
                    i0 = (iu-1)*Su.p + 1;
                    i1 = i0 + Su.p;
                    j0 = (jv-1)*Su.q + 1;
                    j1 = j0 + Su.q;

                    Pij = Su.P(i0:i1, j0:j1, :);
                    Wij = Su.W(i0:i1, j0:j1);
                    Uij = [zeros(1, Su.p+1), ones(1, Su.p+1)];
                    Vij = [zeros(1, Su.q+1), ones(1, Su.q+1)];

                    patches{iu,jv} = struct( ...
                        'surface', geom.NURBSSurface(Pij, Su.p, Su.q, Uij, Vij, Wij), ...
                        'u0', uBreaks(iu), 'u1', uBreaks(iu+1), ...
                        'v0', vBreaks(jv), 'v1', vBreaks(jv+1));
                end
            end
        end

        function S2 = elevateU(obj, t)
            if nargin < 2 || isempty(t), t = 1; end
            if t < 0 || t ~= floor(t)
                error('NURBSSurface:elevateU', 't must be a nonnegative integer.');
            end
            if t == 0
                S2 = obj.copySurfaceOnly();
                return;
            end

            for j = 1:size(obj.P,2)
                Cj = geom.NURBSCurve(squeeze(obj.P(:,j,:)), obj.p, obj.U, obj.W(:,j));
                Ce = Cj.elevate(t);
                if j == 1
                    nnew = size(Ce.P,1);
                    P2 = zeros(nnew, size(obj.P,2), 3);
                    W2 = zeros(nnew, size(obj.P,2));
                    U2 = Ce.U;
                end
                P2(:,j,:) = Ce.P;
                W2(:,j) = Ce.W;
            end

            S2 = geom.NURBSSurface(P2, obj.p + t, obj.q, U2, obj.V, W2);
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance = obj.trimTolerance;
        end

        function S2 = elevateV(obj, t)
            if nargin < 2 || isempty(t), t = 1; end
            if t < 0 || t ~= floor(t)
                error('NURBSSurface:elevateV', 't must be a nonnegative integer.');
            end
            if t == 0
                S2 = obj.copySurfaceOnly();
                return;
            end

            for i = 1:size(obj.P,1)
                Ci = geom.NURBSCurve(squeeze(obj.P(i,:,:)), obj.q, obj.V, obj.W(i,:).');
                Ce = Ci.elevate(t);
                if i == 1
                    mnew = size(Ce.P,1);
                    P2 = zeros(size(obj.P,1), mnew, 3);
                    W2 = zeros(size(obj.P,1), mnew);
                    V2 = Ce.U;
                end
                P2(i,:,:) = Ce.P;
                W2(i,:) = Ce.W.';
            end

            S2 = geom.NURBSSurface(P2, obj.p, obj.q + t, obj.U, V2, W2);
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance = obj.trimTolerance;
        end

        function [S2, maxErr] = reduceDegreeU(obj, numTimes, tol, nSample)
            if nargin < 2 || isempty(numTimes), numTimes = 1; end
            if nargin < 3 || isempty(tol), tol = inf; end
            if nargin < 4 || isempty(nSample), nSample = 200; end

            maxErr = 0;
            Scur = obj;

            for step = 1:numTimes
                if Scur.p <= 0
                    error('NURBSSurface:reduceDegreeU', 'Cannot reduce p below zero.');
                end

                for j = 1:size(Scur.P,2)
                    Cj = geom.NURBSCurve(squeeze(Scur.P(:,j,:)), Scur.p, Scur.U, Scur.W(:,j));
                    [Cr, errj] = Cj.reduceDegree(1, tol, nSample);
                    if j == 1
                        nnew = size(Cr.P,1);
                        P2 = zeros(nnew, size(Scur.P,2), 3);
                        W2 = zeros(nnew, size(Scur.P,2));
                        U2 = Cr.U;
                    end
                    P2(:,j,:) = Cr.P;
                    W2(:,j) = Cr.W;
                    maxErr = max(maxErr, errj);
                end

                Scur = geom.NURBSSurface(P2, Scur.p - 1, Scur.q, U2, Scur.V, W2);
            end

            S2 = Scur;
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance = obj.trimTolerance;
        end

        function [S2, maxErr] = reduceDegreeV(obj, numTimes, tol, nSample)
            if nargin < 2 || isempty(numTimes), numTimes = 1; end
            if nargin < 3 || isempty(tol), tol = inf; end
            if nargin < 4 || isempty(nSample), nSample = 200; end

            maxErr = 0;
            Scur = obj;

            for step = 1:numTimes
                if Scur.q <= 0
                    error('NURBSSurface:reduceDegreeV', 'Cannot reduce q below zero.');
                end

                for i = 1:size(Scur.P,1)
                    Ci = geom.NURBSCurve(squeeze(Scur.P(i,:,:)), Scur.q, Scur.V, Scur.W(i,:).');
                    [Cr, erri] = Ci.reduceDegree(1, tol, nSample);
                    if i == 1
                        mnew = size(Cr.P,1);
                        P2 = zeros(size(Scur.P,1), mnew, 3);
                        W2 = zeros(size(Scur.P,1), mnew);
                        V2 = Cr.U;
                    end
                    P2(i,:,:) = Cr.P;
                    W2(i,:) = Cr.W.';
                    maxErr = max(maxErr, erri);
                end

                Scur = geom.NURBSSurface(P2, Scur.p, Scur.q - 1, Scur.U, V2, W2);
            end

            S2 = Scur;
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance = obj.trimTolerance;
        end

        function [S2, removed, maxErr] = removeKnotU(obj, u, numRemove, tol, nSample)
            if nargin < 3 || isempty(numRemove), numRemove = 1; end
            if nargin < 4 || isempty(tol), tol = 1e-10; end
            if nargin < 5 || isempty(nSample), nSample = 200; end

            Scur = obj;
            removed = inf;
            maxErr = 0;

            for j = 1:size(Scur.P,2)
                Cj = geom.NURBSCurve(squeeze(Scur.P(:,j,:)), Scur.p, Scur.U, Scur.W(:,j));
                [Cr, remj, errj] = Cj.removeKnot(u, numRemove, tol, nSample);

                if j == 1
                    removed = remj;
                    nnew = size(Cr.P,1);
                    P2 = zeros(nnew, size(Scur.P,2), 3);
                    W2 = zeros(nnew, size(Scur.P,2));
                    U2 = Cr.U;
                else
                    removed = min(removed, remj);
                end

                P2(:,j,:) = Cr.P;
                W2(:,j) = Cr.W;
                maxErr = max(maxErr, errj);
            end

            S2 = geom.NURBSSurface(P2, Scur.p, Scur.q, U2, Scur.V, W2);
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance = obj.trimTolerance;
        end

        function [S2, removed, maxErr] = removeKnotV(obj, v, numRemove, tol, nSample)
            if nargin < 3 || isempty(numRemove), numRemove = 1; end
            if nargin < 4 || isempty(tol), tol = 1e-10; end
            if nargin < 5 || isempty(nSample), nSample = 200; end

            Scur = obj;
            removed = inf;
            maxErr = 0;

            for i = 1:size(Scur.P,1)
                Ci = geom.NURBSCurve(squeeze(Scur.P(i,:,:)), Scur.q, Scur.V, Scur.W(i,:).');
                [Cr, remi, erri] = Ci.removeKnot(v, numRemove, tol, nSample);

                if i == 1
                    removed = remi;
                    mnew = size(Cr.P,1);
                    P2 = zeros(size(Scur.P,1), mnew, 3);
                    W2 = zeros(size(Scur.P,1), mnew);
                    V2 = Cr.U;
                else
                    removed = min(removed, remi);
                end

                P2(i,:,:) = Cr.P;
                W2(i,:) = Cr.W.';
                maxErr = max(maxErr, erri);
            end

            S2 = geom.NURBSSurface(P2, Scur.p, Scur.q, Scur.U, V2, W2);
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance = obj.trimTolerance;
        end
    end

    methods
        function plot(obj, nu, nv, varargin)
            if nargin < 2, nu = 30; end
            if nargin < 3, nv = 30; end

            pa = inputParser;
            addParameter(pa, 'ShowCP', false);
            addParameter(pa, 'ShowIso', false);
            addParameter(pa, 'ShowTrims', true);
            addParameter(pa, 'Alpha', 0.85);
            addParameter(pa, 'EdgeAlpha', 0.2);
            addParameter(pa, 'FaceColor', [0.3 0.6 0.9]);
            parse(pa, varargin{:});
            opts = pa.Results;

            mesh = obj.isoMesh(nu, nv, 'RespectTrim', true);
            surf(mesh.X, mesh.Y, mesh.Z, ...
                'FaceColor', opts.FaceColor, ...
                'FaceAlpha', opts.Alpha, ...
                'EdgeAlpha', opts.EdgeAlpha, ...
                'EdgeColor', 'k');
            hold on;

            if opts.ShowIso
                n_iso = 8;
                u_lines = linspace(obj.domainU(1), obj.domainU(2), n_iso + 2);
                u_lines = u_lines(2:end-1);
                v_lines = linspace(obj.domainV(1), obj.domainV(2), n_iso + 2);
                v_lines = v_lines(2:end-1);

                for i = 1:numel(u_lines)
                    C = obj.isoCurveV(u_lines(i));
                    C.plot(50, 'ShowCP', false, 'Color', [0.1 0.1 0.5], 'LineWidth', 0.7);
                end
                for j = 1:numel(v_lines)
                    C = obj.isoCurveU(v_lines(j));
                    C.plot(50, 'ShowCP', false, 'Color', [0.5 0.1 0.1], 'LineWidth', 0.7);
                end
            end

            if opts.ShowCP
                for i = 1:size(obj.P, 1)
                    plot3(squeeze(obj.P(i,:,1)), squeeze(obj.P(i,:,2)), squeeze(obj.P(i,:,3)), ...
                        '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
                end
                for j = 1:size(obj.P, 2)
                    plot3(squeeze(obj.P(:,j,1)), squeeze(obj.P(:,j,2)), squeeze(obj.P(:,j,3)), ...
                        '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
                end
                Pf = reshape(obj.P, [], 3);
                plot3(Pf(:,1), Pf(:,2), Pf(:,3), 'o', ...
                    'MarkerSize', 4, 'MarkerFaceColor', 'w', 'Color', [0.3 0.3 0.3]);
            end

            if opts.ShowTrims && obj.isTrimmed()
                T = obj.sampleTrimLoopsXYZ(160);
                for k = 1:numel(T.outer)
                    xyz = T.outer{k};
                    plot3(xyz(:,1), xyz(:,2), xyz(:,3), '-', ...
                        'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
                end
                for k = 1:numel(T.inner)
                    xyz = T.inner{k};
                    plot3(xyz(:,1), xyz(:,2), xyz(:,3), '-', ...
                        'Color', [0.1 0.1 0.9], 'LineWidth', 2.0);
                end
            end

            axis equal;
            grid on;
            lighting phong;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
        end

        function plotNormals(obj, nu, nv, scale)
            if nargin < 2, nu = 8; end
            if nargin < 3, nv = 8; end
            if nargin < 4, scale = 0.05; end

            mesh = obj.isoMesh(nu, nv, 'RespectTrim', true);
            quiver3(mesh.X, mesh.Y, mesh.Z, ...
                scale * mesh.normals(:,:,1), ...
                scale * mesh.normals(:,:,2), ...
                scale * mesh.normals(:,:,3), ...
                0, 'b');
        end

        function plotTrimUV(obj, nPerLoop)
            if nargin < 2 || isempty(nPerLoop), nPerLoop = 160; end
            UV = obj.sampleTrimLoops(nPerLoop);
            hold on; grid on; axis equal;
            xlabel('u'); ylabel('v');
            xlim(obj.domainU);
            ylim(obj.domainV);

            rectangle('Position', [obj.domainU(1), obj.domainV(1), ...
                                   obj.domainU(2)-obj.domainU(1), ...
                                   obj.domainV(2)-obj.domainV(1)], ...
                      'EdgeColor', [0.5 0.5 0.5], 'LineStyle', '--');

            for k = 1:numel(UV.outer)
                uv = UV.outer{k};
                plot(uv(:,1), uv(:,2), 'r-', 'LineWidth', 2.0);
            end
            for k = 1:numel(UV.inner)
                uv = UV.inner{k};
                plot(uv(:,1), uv(:,2), 'b-', 'LineWidth', 2.0);
            end
            title('Trim loops in UV space');
        end
    end

    methods
        function obj2 = copySurfaceOnly(obj)
            obj2 = geom.NURBSSurface(obj.P, obj.p, obj.q, obj.U, obj.V, obj.W);
            obj2.trimOuterLoops = obj.trimOuterLoops;
            obj2.trimInnerLoops = obj.trimInnerLoops;
            obj2.trimTolerance = obj.trimTolerance;
        end
    end

    methods (Access = private)
        function u = clampU(obj, u)
            u = max(obj.domainU(1), min(obj.domainU(2), u));
        end

        function v = clampV(obj, v)
            v = max(obj.domainV(1), min(obj.domainV(2), v));
        end

        function [u0, v0] = coarseSearchPoint(obj, P, nu, nv)
            u_c = linspace(obj.domainU(1), obj.domainU(2), nu);
            v_c = linspace(obj.domainV(1), obj.domainV(2), nv);

            best_d2 = inf;
            u0 = u_c(1);
            v0 = v_c(1);

            for i = 1:nu
                for j = 1:nv
                    pt = obj.evaluate(u_c(i), v_c(j));
                    d2 = sum((pt - P).^2);
                    if d2 < best_d2
                        best_d2 = d2;
                        u0 = u_c(i);
                        v0 = v_c(j);
                    end
                end
            end
        end
    end

    methods (Static)
        function t = makeSpacing(t0, t1, n, type)
            switch lower(type)
                case 'cosine'
                    beta = linspace(0, pi, n);
                    s = 0.5 * (1 - cos(beta));
                case 'cosine_half'
                    beta = linspace(0, pi/2, n);
                    s = sin(beta);
                case 'cosine_wrap'
                    n1 = ceil(n/2);
                    n2 = n - n1 + 1;
                    b1 = linspace(0, pi, n1);
                    s1 = 0.5 * (1 - cos(b1)) * 0.5;
                    b2 = linspace(0, pi, n2);
                    s2 = 0.5 + 0.5 * (1 - cos(b2)) * 0.5;
                    s = [s1, s2(2:end)];
                otherwise
                    s = linspace(0, 1, n);
            end
            t = t0 + s * (t1 - t0);
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

        function loops = normalizeTrimLoopSet(loopsIn)
            if isempty(loopsIn)
                loops = {};
                return;
            end
            if ~iscell(loopsIn)
                loopsIn = {loopsIn};
            end
            loops = cell(size(loopsIn));
            for k = 1:numel(loopsIn)
                loops{k} = geom.NURBSSurface.loopToPolylineUV(loopsIn{k}, 200);
            end
        end

        function uv = loopToPolylineUV(loop, nSamp)
            if nargin < 2 || isempty(nSamp), nSamp = 200; end

            if isa(loop, 'geom.NURBSCurve')
                u = linspace(loop.domain(1), loop.domain(2), nSamp);
                xyz = loop.evaluate(u);
                uv = xyz(:,1:2);
            elseif isnumeric(loop) && size(loop,2) == 2
                uv = loop;
            else
                error('NURBSSurface:loopToPolylineUV', ...
                    'Trim loop must be Nx2 numeric or geom.NURBSCurve in UV plane.');
            end

            if size(uv,1) < 3
                error('NURBSSurface:loopToPolylineUV', ...
                    'Trim loop needs at least 3 points.');
            end

            if norm(uv(1,:) - uv(end,:)) > 1e-12
                uv = [uv; uv(1,:)];
            end
        end

        function tf = pointInPolygon2D(poly, pt, tol)
            if nargin < 3 || isempty(tol), tol = 1e-10; end
            x = pt(1);
            y = pt(2);
            xv = poly(:,1);
            yv = poly(:,2);

            % Boundary check
            for i = 1:size(poly,1)-1
                a = poly(i,:);
                b = poly(i+1,:);
                if geom.NURBSSurface.pointOnSegment2D(pt, a, b, tol)
                    tf = true;
                    return;
                end
            end

            tf = inpolygon(x, y, xv, yv);
        end

        function tf = pointOnSegment2D(p, a, b, tol)
            ap = p - a;
            ab = b - a;
            lab2 = dot(ab, ab);
            if lab2 < eps
                tf = norm(ap) <= tol;
                return;
            end
            t = dot(ap, ab) / lab2;
            if t < -tol || t > 1+tol
                tf = false;
                return;
            end
            proj = a + min(max(t,0),1) * ab;
            tf = norm(p - proj) <= tol;
        end

        function uv = resampleTrimLoop(loop, nPerLoop)
            if isa(loop, 'geom.NURBSCurve')
                us = linspace(loop.domain(1), loop.domain(2), nPerLoop);
                xyz = loop.evaluate(us);
                uv = xyz(:,1:2);
                if norm(uv(1,:) - uv(end,:)) > 1e-12
                    uv = [uv; uv(1,:)];
                end
            else
                uv0 = geom.NURBSSurface.loopToPolylineUV(loop, nPerLoop);
                d = [0; cumsum(vecnorm(diff(uv0,1,1),2,2))];
                if d(end) < eps
                    uv = uv0;
                    return;
                end
                s = linspace(0, d(end), nPerLoop).';
                uv = [interp1(d, uv0(:,1), s, 'linear'), ...
                      interp1(d, uv0(:,2), s, 'linear')];
                if norm(uv(1,:) - uv(end,:)) > 1e-12
                    uv = [uv; uv(1,:)];
                end
            end
        end

        function [leftLoops, rightLoops] = splitTrimLoopsU(loops, u0)
            % Practical placeholder: preserve loops on both children by
            % renormalizing coordinates. This is conservative and useful for
            % demos / masking, though it is not a full topological trim split.
            leftLoops = {};
            rightLoops = {};
            for k = 1:numel(loops)
                uv = loops{k};
                uvL = uv;
                uvR = uv;
                uvL(:,1) = min(uvL(:,1), u0);
                uvR(:,1) = max(uvR(:,1), u0);

                % map to [0,1]
                if u0 > 0
                    uvL(:,1) = uvL(:,1) / u0;
                else
                    uvL(:,1) = 0;
                end
                if (1-u0) > 0
                    uvR(:,1) = (uvR(:,1) - u0) / (1-u0);
                else
                    uvR(:,1) = 0;
                end
                leftLoops{end+1} = uvL; %#ok<AGROW>
                rightLoops{end+1} = uvR; %#ok<AGROW>
            end
        end

        function S = globalInterpNet(Q, p, q, methodU, methodV)
            if nargin < 4 || isempty(methodU), methodU = 'centripetal'; end
            if nargin < 5 || isempty(methodV), methodV = 'centripetal'; end

            if size(Q,3) ~= 3
                error('NURBSSurface:globalInterpNet', 'Q must be nu x nv x 3.');
            end

            nv = size(Q,2);

            for j = 1:nv
                Cj = geom.NURBSCurve.globalInterp(squeeze(Q(:,j,:)), p, methodU);
                if j == 1
                    Pu = zeros(size(Cj.P,1), nv, 3);
                    Wu = zeros(size(Cj.P,1), nv);
                    U = Cj.U;
                end
                Pu(:,j,:) = Cj.P;
                Wu(:,j) = Cj.W;
            end

            for i = 1:size(Pu,1)
                Ci = geom.NURBSCurve.globalInterp(squeeze(Pu(i,:,:)), q, methodV);
                if i == 1
                    P = zeros(size(Pu,1), size(Ci.P,1), 3);
                    W = zeros(size(Pu,1), size(Ci.P,1));
                    V = Ci.U;
                end
                P(i,:,:) = Ci.P;
                W(i,:) = Ci.W.';
            end

            S = geom.NURBSSurface(P, p, q, U, V, W);
        end

        function S = globalLeastSquaresFitNet(Q, p, q, nCtrlU, nCtrlV, methodU, methodV)
            if nargin < 6 || isempty(methodU), methodU = 'centripetal'; end
            if nargin < 7 || isempty(methodV), methodV = 'centripetal'; end

            if size(Q,3) ~= 3
                error('NURBSSurface:globalLeastSquaresFitNet', 'Q must be nu x nv x 3.');
            end

            nv = size(Q,2);

            for j = 1:nv
                Cj = geom.NURBSCurve.globalLeastSquaresFit(squeeze(Q(:,j,:)), p, nCtrlU, methodU);
                if j == 1
                    Pu = zeros(size(Cj.P,1), nv, 3);
                    Wu = zeros(size(Cj.P,1), nv);
                    U = Cj.U;
                end
                Pu(:,j,:) = Cj.P;
                Wu(:,j) = Cj.W;
            end

            for i = 1:size(Pu,1)
                Ci = geom.NURBSCurve.globalLeastSquaresFit(squeeze(Pu(i,:,:)), q, nCtrlV, methodV);
                if i == 1
                    P = zeros(size(Pu,1), size(Ci.P,1), 3);
                    W = zeros(size(Pu,1), size(Ci.P,1));
                    V = Ci.U;
                end
                P(i,:,:) = Ci.P;
                W(i,:) = Ci.W.';
            end

            S = geom.NURBSSurface(P, p, q, U, V, W);
        end

        function S = ruled(C1, C2)
            if C1.p ~= C2.p || numel(C1.U) ~= numel(C2.U) || any(abs(C1.U - C2.U) > 1e-12)
                error('NURBSSurface:ruled', 'Curves must have matching degree and knot vector.');
            end
            if size(C1.P,1) ~= size(C2.P,1)
                error('NURBSSurface:ruled', 'Curves must have same number of control points.');
            end

            n = size(C1.P,1);
            P = zeros(n, 2, 3);
            W = zeros(n, 2);

            P(:,1,:) = C1.P;
            P(:,2,:) = C2.P;
            W(:,1) = C1.W;
            W(:,2) = C2.W;

            U = C1.U;
            V = [0 0 1 1];

            S = geom.NURBSSurface(P, C1.p, 1, U, V, W);
        end

        function S = bilinearCoons(P00, P10, P01, P11)
            P = zeros(2,2,3);
            P(1,1,:) = reshape(P00,1,1,3);
            P(2,1,:) = reshape(P10,1,1,3);
            P(1,2,:) = reshape(P01,1,1,3);
            P(2,2,:) = reshape(P11,1,1,3);

            W = ones(2,2);
            U = [0 0 1 1];
            V = [0 0 1 1];

            S = geom.NURBSSurface(P, 1, 1, U, V, W);
        end

        function [curvesOut, pmax, Uunion] = makeCompatibleCurves(curves)
            if ~iscell(curves), curves = {curves}; end
            if isempty(curves)
                error('NURBSSurface:makeCompatibleCurves', 'Need at least one curve.');
            end

            pmax = max(cellfun(@(c) c.p, curves));
            curvesOut = curves;

            for k = 1:numel(curvesOut)
                if curvesOut{k}.p < pmax
                    curvesOut{k} = curvesOut{k}.elevate(pmax - curvesOut{k}.p);
                end
            end

            Uunion = curvesOut{1}.U;
            for k = 2:numel(curvesOut)
                Uunion = union(Uunion, curvesOut{k}.U);
            end
            Uunion = sort(Uunion(:).');

            for val = Uunion
                multTarget = max(cellfun(@(c) sum(abs(c.U - val) < 1e-12), curvesOut));
                for k = 1:numel(curvesOut)
                    multK = sum(abs(curvesOut{k}.U - val) < 1e-12);
                    if multK < multTarget
                        curvesOut{k} = curvesOut{k}.insertKnot(val, multTarget - multK);
                    end
                end
            end

            Uunion = curvesOut{1}.U;
        end

        function [profilesC, guidesC, uPar, vPar] = makeCompatibleNetwork(profileCurves, guideCurves)
        % Practical compatibility helper for rectangular curve networks.
        %
        % Returns:
        %   profilesC : compatible profile curves
        %   guidesC   : compatible guide curves
        %   uPar      : normalized guide stations
        %   vPar      : normalized profile stations
            if ~iscell(profileCurves), profileCurves = {profileCurves}; end
            if ~iscell(guideCurves), guideCurves = {guideCurves}; end
            if isempty(profileCurves) || isempty(guideCurves)
                error('NURBSSurface:makeCompatibleNetwork', ...
                    'Need profile and guide curves.');
            end

            [profilesC, ~, ~] = geom.NURBSSurface.makeCompatibleCurves(profileCurves);
            [guidesC, ~, ~] = geom.NURBSSurface.makeCompatibleCurves(guideCurves);

            uPar = linspace(0,1,numel(guideCurves)).';
            vPar = linspace(0,1,numel(profileCurves)).';
        end

        function v = sectionParameters(curves, method)
            if nargin < 2 || isempty(method), method = 'centripetal'; end
            if ~iscell(curves), curves = {curves}; end

            nsec = numel(curves);
            Q = zeros(nsec, 3);
            for k = 1:nsec
                uk = 0.5 * (curves{k}.domain(1) + curves{k}.domain(2));
                Q(k,:) = curves{k}.evaluate(uk);
            end
            v = geom.NURBSCurve.parameterizeData(Q, method, min(3, nsec-1));
        end

        function S = loft(curves, q, method, sectionParams)
            if nargin < 2 || isempty(q), q = 3; end
            if nargin < 3 || isempty(method), method = 'centripetal'; end
            if ~iscell(curves), curves = {curves}; end
            if numel(curves) < 2
                error('NURBSSurface:loft', 'Need at least two section curves.');
            end

            [curvesC, pU, U] = geom.NURBSSurface.makeCompatibleCurves(curves);
            nsec = numel(curvesC);
            nCtrlU = size(curvesC{1}.P,1);

            if q > nsec-1
                error('NURBSSurface:loft', ...
                    'Loft degree q=%d exceeds number of sections-1 = %d.', q, nsec-1);
            end

            if nargin < 4 || isempty(sectionParams)
                vpar = geom.NURBSSurface.sectionParameters(curvesC, method);
            else
                vpar = sectionParams(:);
                if numel(vpar) ~= nsec
                    error('NURBSSurface:loft', 'sectionParams must match number of curves.');
                end
            end

            V = geom.NURBSCurve.averagingKnotVector(vpar, q);

            A = zeros(nsec, nsec);
            nV = nsec - 1;
            for r = 1:nsec
                span = geom.BasisFunctions.FindSpan(nV, q, vpar(r), V);
                N = geom.BasisFunctions.BasisFuns(span, vpar(r), q, V);
                cols = (span-q):span;
                A(r, cols) = N;
            end

            P = zeros(nCtrlU, nsec, 3);
            W = zeros(nCtrlU, nsec);

            for i = 1:nCtrlU
                H = zeros(nsec, 4);
                for s = 1:nsec
                    w = curvesC{s}.W(i);
                    H(s,:) = [w * curvesC{s}.P(i,:), w];
                end

                Qh = A \ H;
                if any(Qh(:,4) <= 0)
                    error('NURBSSurface:loft', ...
                        'Lofting produced nonpositive weights.');
                end

                W(i,:) = Qh(:,4).';
                P(i,:,:) = Qh(:,1:3) ./ Qh(:,4);
            end

            S = geom.NURBSSurface(P, pU, q, U, V, W);
        end

        function S = coons(Cu0, Cu1, Cv0, Cv1, p, q, nu, nv)
            if nargin < 5 || isempty(p), p = 3; end
            if nargin < 6 || isempty(q), q = 3; end
            if nargin < 7 || isempty(nu), nu = 21; end
            if nargin < 8 || isempty(nv), nv = 21; end

            P00 = Cu0.evaluate(Cu0.domain(1));
            P10 = Cu0.evaluate(Cu0.domain(2));
            P01 = Cu1.evaluate(Cu1.domain(1));
            P11 = Cu1.evaluate(Cu1.domain(2));

            Q00 = Cv0.evaluate(Cv0.domain(1));
            Q01 = Cv0.evaluate(Cv0.domain(2));
            Q10 = Cv1.evaluate(Cv1.domain(1));
            Q11 = Cv1.evaluate(Cv1.domain(2));

            if max([norm(P00-Q00), norm(P01-Q01), norm(P10-Q10), norm(P11-Q11)]) > 1e-6
                error('NURBSSurface:coons', ...
                    'Boundary curve corners are inconsistent.');
            end

            ug = linspace(0,1,nu);
            vg = linspace(0,1,nv);
            Q = zeros(nu, nv, 3);

            for i = 1:nu
                u = ug(i);

                Pu0 = Cu0.evaluate(Cu0.domain(1) + u*(Cu0.domain(2)-Cu0.domain(1)));
                Pu1 = Cu1.evaluate(Cu1.domain(1) + u*(Cu1.domain(2)-Cu1.domain(1)));

                for j = 1:nv
                    v = vg(j);

                    Pv0 = Cv0.evaluate(Cv0.domain(1) + v*(Cv0.domain(2)-Cv0.domain(1)));
                    Pv1 = Cv1.evaluate(Cv1.domain(1) + v*(Cv1.domain(2)-Cv1.domain(1)));

                    Suv = (1-v) * Pu0 + v * Pu1;
                    Svv = (1-u) * Pv0 + u * Pv1;
                    Sbil = (1-u)*(1-v)*P00 + u*(1-v)*P10 + (1-u)*v*P01 + u*v*P11;

                    Q(i,j,:) = Suv + Svv - Sbil;
                end
            end

            S = geom.NURBSSurface.globalInterpNet(Q, p, q, 'chord', 'chord');
        end

        function S = gordon(profileCurves, guideCurves, p, q, nu, nv)
            if nargin < 3 || isempty(p), p = 3; end
            if nargin < 4 || isempty(q), q = 3; end
            if nargin < 5 || isempty(nu), nu = 25; end
            if nargin < 6 || isempty(nv), nv = 25; end

            if ~iscell(profileCurves), profileCurves = {profileCurves}; end
            if ~iscell(guideCurves), guideCurves = {guideCurves}; end
            if numel(profileCurves) < 2 || numel(guideCurves) < 2
                error('NURBSSurface:gordon', ...
                    'Need at least two profile curves and two guide curves.');
            end

            Sp = geom.NURBSSurface.loft(profileCurves, q, 'centripetal');
            Sg = geom.NURBSSurface.loft(guideCurves, p, 'centripetal').swapUV();

            vp = linspace(0,1,numel(profileCurves));
            ug = linspace(0,1,numel(guideCurves));
            X = zeros(numel(ug), numel(vp), 3);

            for i = 1:numel(ug)
                for j = 1:numel(vp)
                    X(i,j,:) = Sp.evaluate(ug(i), vp(j));
                end
            end

            Si = geom.NURBSSurface.globalInterpNet(X, p, q, 'uniform', 'uniform');

            ugrid = linspace(0,1,nu);
            vgrid = linspace(0,1,nv);
            Q = zeros(nu, nv, 3);

            for i = 1:nu
                for j = 1:nv
                    Pp = Sp.evaluate(ugrid(i), vgrid(j));
                    Pg = Sg.evaluate(ugrid(i), vgrid(j));
                    Pi = Si.evaluate(ugrid(i), vgrid(j));
                    Q(i,j,:) = Pp + Pg - Pi;
                end
            end

            S = geom.NURBSSurface.globalInterpNet(Q, p, q, 'chord', 'chord');
        end

        function S = multiGordon(profileFamilies, guideFamilies, weights, p, q, nu, nv)
        % Higher-order multi-network blending variant.
        %
        % profileFamilies, guideFamilies are cell arrays of curve-network
        % pairs. Each pair is blended Gordon-style, then combined with weights.
            if nargin < 3 || isempty(weights)
                weights = ones(numel(profileFamilies),1);
            end
            if nargin < 4 || isempty(p), p = 3; end
            if nargin < 5 || isempty(q), q = 3; end
            if nargin < 6 || isempty(nu), nu = 25; end
            if nargin < 7 || isempty(nv), nv = 25; end

            weights = weights(:);
            nfam = numel(profileFamilies);
            if numel(guideFamilies) ~= nfam || numel(weights) ~= nfam
                error('NURBSSurface:multiGordon', ...
                    'Families and weights must have matching lengths.');
            end

            Sfam = cell(nfam,1);
            for k = 1:nfam
                Sfam{k} = geom.NURBSSurface.gordon( ...
                    profileFamilies{k}, guideFamilies{k}, p, q, nu, nv);
            end

            ug = linspace(0,1,nu);
            vg = linspace(0,1,nv);
            Q = zeros(nu, nv, 3);

            ws = sum(weights);
            if abs(ws) < eps
                error('NURBSSurface:multiGordon', 'Weights sum to zero.');
            end

            for i = 1:nu
                for j = 1:nv
                    acc = zeros(1,3);
                    for k = 1:nfam
                        acc = acc + weights(k) * Sfam{k}.evaluate(ug(i), vg(j));
                    end
                    Q(i,j,:) = acc / ws;
                end
            end

            S = geom.NURBSSurface.globalInterpNet(Q, p, q, 'chord', 'chord');
        end

        function S = revolve(C, axisPoint, axisDir, thetaTotal, q, nSections)
            if nargin < 4 || isempty(thetaTotal), thetaTotal = 2*pi; end
            if nargin < 5 || isempty(q), q = 3; end
            if nargin < 6 || isempty(nSections), nSections = 9; end

            axisPoint = axisPoint(:).';
            axisDir = axisDir(:).';
            axisDir = axisDir / norm(axisDir);

            th = linspace(0, thetaTotal, nSections);
            curves = cell(1, nSections);

            for k = 1:nSections
                T = geom.NURBSSurface.axisRotationTransform(axisPoint, axisDir, th(k));
                curves{k} = C.transform(T);
            end

            S = geom.NURBSSurface.loft(curves, q, 'chord');
        end

        function S = sweep(profileCurve, spineCurve, q, nStations, upVec)
            if nargin < 3 || isempty(q), q = 3; end
            if nargin < 4 || isempty(nStations), nStations = 9; end
            if nargin < 5 || isempty(upVec), upVec = [0 0 1]; end

            upVec = upVec(:).';
            upVec = upVec / norm(upVec);

            ts = spineCurve.arcLengthParam(nStations);
            sections = cell(1, nStations);

            for k = 1:nStations
                tk = ts(k);
                P0 = spineCurve.evaluate(tk);
                T0 = spineCurve.tangent(tk);
                T0 = T0 / norm(T0);

                B0 = cross(T0, upVec);
                if norm(B0) < 1e-8
                    alt = [1 0 0];
                    if abs(dot(T0, alt)) > 0.9
                        alt = [0 1 0];
                    end
                    B0 = cross(T0, alt);
                end
                B0 = B0 / norm(B0);
                N0 = cross(B0, T0);
                N0 = N0 / norm(N0);

                R = [T0(:), N0(:), B0(:)];
                H = eye(4);
                H(1:3,1:3) = R;
                H(1:3,4) = P0(:);

                sections{k} = profileCurve.transform(H);
            end

            S = geom.NURBSSurface.loft(sections, q, 'chord');
        end
    end

    methods (Static, Access = private)
        function T = axisRotationTransform(point, axisDir, theta)
            axisDir = axisDir(:) / norm(axisDir);
            x = axisDir(1); y = axisDir(2); z = axisDir(3);
            c = cos(theta);
            s = sin(theta);
            C = 1 - c;

            R = [x*x*C + c,   x*y*C - z*s, x*z*C + y*s; ...
                 y*x*C + z*s, y*y*C + c,   y*z*C - x*s; ...
                 z*x*C - y*s, z*y*C + x*s, z*z*C + c];

            p = point(:);
            t = p - R*p;

            T = eye(4);
            T(1:3,1:3) = R;
            T(1:3,4) = t;
        end
    end
end