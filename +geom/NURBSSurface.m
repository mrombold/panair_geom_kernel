classdef NURBSSurface < handle
% NURBSSURFACE  Non-Uniform Rational B-Spline surface.
%
%   Cleaned tensor-product NURBS surface class with exact rational
%   evaluation, exact rational derivative tensor, analytic first/second
%   partials, curvature, projection, exact split, iso-curves, and subpatch
%   extraction by repeated exact splitting.

    properties
        P       % [(n+1) x (m+1) x 3] Cartesian control net
        W       % [(n+1) x (m+1)] weights
        U       % u knot vector
        V       % v knot vector
        p       % degree in u
        q       % degree in v
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
                error('NURBSSurface: requires at least P, p, q.');
            end

            if isstruct(P)
                P = cat(3, P.x, P.y, P.z);
            end
            sz = size(P);
            if numel(sz) ~= 3 || sz(3) ~= 3
                error('NURBSSurface: P must be an (n+1)x(m+1)x3 array.');
            end
            if ~isscalar(p) || p < 0 || p ~= floor(p)
                error('NURBSSurface: p must be a nonnegative integer.');
            end
            if ~isscalar(q) || q < 0 || q ~= floor(q)
                error('NURBSSurface: q must be a nonnegative integer.');
            end

            n = sz(1) - 1;
            m = sz(2) - 1;
            if p > n
                error('NURBSSurface: degree p=%d exceeds n=%d.', p, n);
            end
            if q > m
                error('NURBSSurface: degree q=%d exceeds m=%d.', q, m);
            end

            if nargin < 6 || isempty(W)
                W = ones(n+1, m+1);
            else
                if ~isequal(size(W), [n+1, m+1])
                    error('NURBSSurface: W must be size (%d x %d).', n+1, m+1);
                end
                if any(W(:) <= 0)
                    error('NURBSSurface: weights must be strictly positive.');
                end
            end

            if nargin < 4 || isempty(U)
                U = geom.BasisFunctions.MakeUniformKnotVector(n, p);
            else
                U = U(:)';
            end
            if nargin < 5 || isempty(V)
                V = geom.BasisFunctions.MakeUniformKnotVector(m, q);
            else
                V = V(:)';
            end

            expectedU = n + p + 2;
            expectedV = m + q + 2;
            if numel(U) ~= expectedU
                error('NURBSSurface: U length must be n+p+2 = %d, got %d.', expectedU, numel(U));
            end
            if numel(V) ~= expectedV
                error('NURBSSurface: V length must be m+q+2 = %d, got %d.', expectedV, numel(V));
            end
            if any(diff(U) < 0)
                error('NURBSSurface: U must be nondecreasing.');
            end
            if any(diff(V) < 0)
                error('NURBSSurface: V must be nondecreasing.');
            end

            obj.P = P;
            obj.W = W;
            obj.U = U;
            obj.V = V;
            obj.p = p;
            obj.q = q;
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
        function C = evaluate(obj, u, v)
            u = u(:);
            v = v(:);
            if isscalar(u), u = repmat(u, numel(v), 1); end
            if isscalar(v), v = repmat(v, numel(u), 1); end
            assert(numel(u) == numel(v), ...
                'NURBSSurface.evaluate: u and v must match in length.');

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
                    for i = 0:obj.p
                        ii = uspan - obj.p + i;
                        jj = vspan - obj.q + l;
                        wij = obj.W(ii, jj);
                        Pij = squeeze(obj.P(ii, jj, :))';
                        temp = temp + Nu(i+1) * [Pij * wij, wij];
                    end
                    Sw = Sw + Nv(l+1) * temp;
                end
                C(k,:) = Sw(1:3) / Sw(4);
            end
        end

        function SKL = derivatives(obj, u, v, d)
        % DERIVATIVES  Exact rational surface derivatives up to order d.
        % Returns SKL{k+1,l+1} = d^(k+l)S / du^k dv^l.
            if nargin < 4 || isempty(d)
                d = 1;
            end
            if numel(u) ~= 1 || numel(v) ~= 1
                error('NURBSSurface:derivatives', 'Use scalar u,v for derivatives().');
            end

            u = obj.clampU(u);
            v = obj.clampV(v);

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
                            Pij = squeeze(obj.P(ii, jj, :))';
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
                error('geom.NURBSSurface:ZeroWeight', ...
                    'Surface weight function vanished at (u,v).');
            end

            for k = 0:d
                for l = 0:d-k
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
            Su = SKL{2,1};
            Sv = SKL{1,2};
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
                K = NaN; H = NaN;
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
                K = NaN; H = NaN;
                return;
            end

            K = (L * N2 - M * M) / den;
            H = (E * N2 - 2 * F * M + G * L) / (2 * den);
        end

        function [u_c, v_c, pt_c, d_c] = closestPoint(obj, P, u0, v0)
            P = P(:)';

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
            u_flat = UG';
            v_flat = VG';
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
            mesh.u = mesh.u';
            mesh.v = mesh.v';

            Normals = zeros(nu, nv, 3);
            for i = 1:nu
                for j = 1:nv
                    Normals(i, j, :) = obj.normal(u_iso(i), v_iso(j));
                end
            end
            mesh.normals = Normals;

            nquads = (nu - 1) * (nv - 1);
            conn = zeros(nquads, 4);
            q = 1;
            for i = 1:nu-1
                for j = 1:nv-1
                    n1 = (i-1) * nv + j;
                    n2 = i * nv + j;
                    n3 = i * nv + j + 1;
                    n4 = (i-1) * nv + j + 1;
                    conn(q, :) = [n1, n2, n3, n4];
                    q = q + 1;
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
                Xu = sort(Xu(:)');
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
                Xv = sort(Xv(:)');
                q = obj.q;
                V_orig = V2;
                for i = 1:size(P2, 1)
                    Pcol = squeeze(P2(i, :, :));
                    Wcol = W2(i, :)';
                    Ctmp = geom.NURBSCurve(Pcol, q, V_orig, Wcol);
                    Cref = Ctmp.refine(Xv);
                    if i == 1
                        m_new = size(Cref.P, 1);
                        P2_new = zeros(size(P2, 1), m_new, 3);
                        W2_new = zeros(size(P2, 1), m_new);
                        V2 = Cref.U;
                    end
                    P2_new(i, :, :) = Cref.P;
                    W2_new(i, :) = Cref.W';
                end
                P2 = P2_new;
                W2 = W2_new;
            end

            S2 = geom.NURBSSurface(P2, obj.p, obj.q, U2, V2, W2);
        end

        function S2 = flipNormals(obj)
            P2 = flipud(obj.P);
            W2 = flipud(obj.W);
            U2 = (obj.domainU(1) + obj.domainU(2)) - fliplr(obj.U);
            S2 = geom.NURBSSurface(P2, obj.p, obj.q, U2, obj.V, W2);
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
            last = find(abs(U2 - u0) < tol, 1, 'last');
            if isempty(first) || isempty(last)
                error('NURBSSurface:splitU', 'Failed to create split knot.');
            end

            cp_split = last - obj.p;

            P_lo = S_ref.P(1:cp_split, :, :);
            W_lo = S_ref.W(1:cp_split, :);
            U_lo = [U2(1:last), u0];
            U_lo = (U_lo - U_lo(1)) / (U_lo(end) - U_lo(1));

            P_hi = S_ref.P(cp_split:end, :, :);
            W_hi = S_ref.W(cp_split:end, :);
            U_hi = [u0, U2(first:end)];
            U_hi = (U_hi - U_hi(1)) / (U_hi(end) - U_hi(1));

            S_lo = geom.NURBSSurface(P_lo, obj.p, obj.q, U_lo, obj.V, W_lo);
            S_hi = geom.NURBSSurface(P_hi, obj.p, obj.q, U_hi, obj.V, W_hi);
        end

        function [S_lo, S_hi] = splitV(obj, v0)
            St = obj.swapUV();
            [A, B] = St.splitU(v0);
            S_lo = A.swapUV();
            S_hi = B.swapUV();
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
                    Pij = squeeze(obj.P(i+1, jj, :))';
                    Pw(i+1, :) = Pw(i+1, :) + Nv(l+1) * [Pij * wij, wij];
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
                    Pij = squeeze(obj.P(ii, j+1, :))';
                    Pw(j+1, :) = Pw(j+1, :) + Nu(l+1) * [Pij * wij, wij];
                end
            end
            W2 = Pw(:,4);
            P2 = Pw(:,1:3) ./ W2;
            C = geom.NURBSCurve(P2, obj.q, obj.V, W2);
        end

        function Ssub = extractUV(obj, urange, vrange)
            ua = urange(1); ub = urange(2);
            va = vrange(1); vb = vrange(2);
            if ~(ua < ub && va < vb)
                error('geom.NURBSSurface:extractUV', 'Require ua < ub and va < vb.');
            end

            Swork = obj;
            tol = 1e-12;
            if ua > Swork.domainU(1) + tol
                [~, Swork] = Swork.splitU(ua);
            end
            if ub < Swork.domainU(2) - tol
                ub_local = (ub - ua) / (obj.domainU(2) - ua);
                [Swork, ~] = Swork.splitU(ub_local);
            end
            if va > Swork.domainV(1) + tol
                [~, Swork] = Swork.splitV(va);
            end
            if vb < Swork.domainV(2) - tol
                vb_local = (vb - va) / (obj.domainV(2) - va);
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
        end
    end

    methods
        function plot(obj, nu, nv, varargin)
            if nargin < 2, nu = 30; end
            if nargin < 3, nv = 30; end

            pa = inputParser;
            addParameter(pa, 'ShowCP', false);
            addParameter(pa, 'ShowIso', false);
            addParameter(pa, 'Alpha', 0.85);
            addParameter(pa, 'EdgeAlpha', 0.2);
            addParameter(pa, 'FaceColor', [0.3 0.6 0.9]);
            parse(pa, varargin{:});
            opts = pa.Results;

            mesh = obj.isoMesh(nu, nv);
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

            axis equal;
            grid on;
            lighting phong;
            xlabel('X'); ylabel('Y'); zlabel('Z');
        end

        function plotNormals(obj, nu, nv, scale)
            if nargin < 2, nu = 8; end
            if nargin < 3, nv = 8; end
            if nargin < 4, scale = 0.05; end
            mesh = obj.isoMesh(nu, nv);
            quiver3(mesh.X, mesh.Y, mesh.Z, ...
                scale * mesh.normals(:,:,1), ...
                scale * mesh.normals(:,:,2), ...
                scale * mesh.normals(:,:,3), ...
                0, 'b');
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
    end
end