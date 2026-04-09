classdef NURBSSurface < handle
% NURBSSURFACE  Non-Uniform Rational B-Spline surface.
%
%   S(u,v) = sum_i sum_j [ N_{i,p}(u) * N_{j,q}(v) * w_{ij} * P_{ij} ]
%            / sum_i sum_j [ N_{i,p}(u) * N_{j,q}(v) * w_{ij} ]
%
%   Construction:
%     S = geom.NURBSSurface(P, p, q)              % uniform knots, unit weights
%     S = geom.NURBSSurface(P, p, q, U, V)        % user-supplied knot vectors
%     S = geom.NURBSSurface(P, p, q, U, V, W)     % full NURBS with weights
%
%   Inputs:
%     P  - [(n+1) x (m+1) x 3] control net
%          OR struct with fields .x .y .z each [n+1 x m+1]
%     p  - degree in u-direction
%     q  - degree in v-direction
%     U  - knot vector u-direction  (length n+p+2)
%     V  - knot vector v-direction  (length m+q+2)
%     W  - [(n+1) x (m+1)] weight array (optional)
%
%   Key methods:
%     pt            = S.evaluate(u, v)
%     [Su,Sv]       = S.partialDerivatives(u,v)
%     [Suu,Suv,Svv] = S.secondPartials(u,v)
%     [K,H]         = S.curvatures(u,v)
%     [A,B]         = S.splitU(u0)
%     [A,B]         = S.splitV(v0)
%
%   Reference: Piegl & Tiller, The NURBS Book, 2nd ed., Springer 1997.

    properties
        P       % [(n+1) x (m+1) x 3] control net
        W       % [(n+1) x (m+1)] weights
        U       % [1 x n+p+2] knot vector u
        V       % [1 x m+q+2] knot vector v
        p       % degree u
        q       % degree v
    end

    properties (Dependent)
        n       % last index in u-direction
        m       % last index in v-direction
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
    end

    methods
        function v = get.n(obj)
            v = size(obj.P, 1) - 1;
        end
        function v = get.m(obj)
            v = size(obj.P, 2) - 1;
        end
        function v = get.domainU(obj)
            v = [obj.U(1), obj.U(end)];
        end
        function v = get.domainV(obj)
            v = [obj.V(1), obj.V(end)];
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

            nu = numel(u);
            C = zeros(nu, 3);

            for k = 1:nu
                uk = obj.clampU(u(k));
                vk = obj.clampV(v(k));

                uspan = geom.BasisFunctions.FindSpan(obj.n, obj.p, uk, obj.U);
                vspan = geom.BasisFunctions.FindSpan(obj.m, obj.q, vk, obj.V);
                Nu    = geom.BasisFunctions.BasisFuns(uspan, uk, obj.p, obj.U);
                Nv    = geom.BasisFunctions.BasisFuns(vspan, vk, obj.q, obj.V);

                Sw = zeros(1, 4);
                for l = 0:obj.q
                    temp = zeros(1, 4);
                    for i_l = 0:obj.p
                        ii = uspan - obj.p + i_l;
                        jj = vspan - obj.q + l;
                        wij = obj.W(ii, jj);
                        Pij = squeeze(obj.P(ii, jj, :))';
                        temp = temp + Nu(i_l+1) * [Pij * wij, wij];
                    end
                    Sw = Sw + Nv(l+1) * temp;
                end
                C(k, :) = Sw(1:3) / Sw(4);
            end
        end

        function [Su, Sv] = partialDerivatives(obj, u, v)
            u = obj.clampU(u);
            v = obj.clampV(v);

            uspan = geom.BasisFunctions.FindSpan(obj.n, obj.p, u, obj.U);
            vspan = geom.BasisFunctions.FindSpan(obj.m, obj.q, v, obj.V);
            Nu_d  = geom.BasisFunctions.DersBasisFuns(uspan, u, obj.p, 1, obj.U);
            Nv_d  = geom.BasisFunctions.DersBasisFuns(vspan, v, obj.q, 1, obj.V);

            Sw  = zeros(1, 4);
            Swu = zeros(1, 4);
            Swv = zeros(1, 4);

            for l = 0:obj.q
                temp  = zeros(1, 4);
                tempu = zeros(1, 4);
                for i_l = 0:obj.p
                    ii  = uspan - obj.p + i_l;
                    jj  = vspan - obj.q + l;
                    wij = obj.W(ii, jj);
                    Pij = squeeze(obj.P(ii, jj, :))';
                    Pw  = [Pij * wij, wij];
                    temp  = temp  + Nu_d(1, i_l+1) * Pw;
                    tempu = tempu + Nu_d(2, i_l+1) * Pw;
                end
                Sw  = Sw  + Nv_d(1, l+1) * temp;
                Swu = Swu + Nv_d(1, l+1) * tempu;
                Swv = Swv + Nv_d(2, l+1) * temp;
            end

            w  = Sw(4);
            A  = Sw(1:3);
            Au = Swu(1:3);
            Av = Swv(1:3);
            wu = Swu(4);
            wv = Swv(4);

            Su = (Au - wu * A / w) / w;
            Sv = (Av - wv * A / w) / w;
        end

        function [Suu, Suv, Svv] = secondPartials(obj, u, v)
            if exist('geom.SurfaceEval', 'class') == 8
                [~, ~, Suu, Suv, Svv] = geom.SurfaceEval.firstSecond(obj, u, v);
                return;
            end

            du = 1e-6 * max(obj.U(end) - obj.U(1), 1);
            dv = 1e-6 * max(obj.V(end) - obj.V(1), 1);

            [Su_p, ~] = obj.partialDerivatives(obj.clampU(u + du), v);
            [Su_m, ~] = obj.partialDerivatives(obj.clampU(u - du), v);
            Suu = (Su_p - Su_m) / (2 * du);

            [Su_vp, ~] = obj.partialDerivatives(u, obj.clampV(v + dv));
            [Su_vm, ~] = obj.partialDerivatives(u, obj.clampV(v - dv));
            Suv = (Su_vp - Su_vm) / (2 * dv);

            [~, Sv_p] = obj.partialDerivatives(u, obj.clampV(v + dv));
            [~, Sv_m] = obj.partialDerivatives(u, obj.clampV(v - dv));
            Svv = (Sv_p - Sv_m) / (2 * dv);
        end

        function N = normal(obj, u, v)
            [Su, Sv] = obj.partialDerivatives(u, v);
            N = cross(Su, Sv);
            nrm = norm(N);
            if nrm > eps
                N = N / nrm;
            end
        end

        function [K, H] = curvatures(obj, u, v)
            if exist('geom.SurfaceEval', 'class') == 8
                [K, H] = geom.SurfaceEval.curvatures(obj, u, v);
                return;
            end

            [Su, Sv] = obj.partialDerivatives(u, v);
            N = cross(Su, Sv);
            N = N / norm(N);
            [Suu, Suv, Svv] = obj.secondPartials(u, v);

            E = dot(Su, Su);
            F = dot(Su, Sv);
            G = dot(Sv, Sv);

            L = dot(Suu, N);
            M = dot(Suv, N);
            Nf = dot(Svv, N);

            den = E * G - F^2;
            K = (L * Nf - M^2) / den;
            H = (E * Nf - 2 * F * M + G * L) / (2 * den);
        end
    end

    methods
        function [u_iso, v_iso, pts] = isoGrid(obj, nu, nv)
            u_iso = linspace(obj.U(1), obj.U(end), nu);
            v_iso = linspace(obj.V(1), obj.V(end), nv);
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

            u_iso = geom.NURBSSurface.makeSpacing(obj.U(1), obj.U(end), nu, opts.SpacingU);
            v_iso = geom.NURBSSurface.makeSpacing(obj.V(1), obj.V(end), nv, opts.SpacingV);

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
            U2 = (obj.U(1) + obj.U(end)) - fliplr(obj.U);
            S2 = geom.NURBSSurface(P2, obj.p, obj.q, U2, obj.V, W2);
        end

        function [S_lo, S_hi] = splitU(obj, u0)
            tol = 1e-12;
            u0 = obj.clampU(u0);
            if u0 <= obj.U(1) + tol || u0 >= obj.U(end) - tol
                error('NURBSSurface:splitU', ...
                      'u0 must lie strictly inside the U domain.');
            end

            s = sum(abs(obj.U - u0) < tol);
            if s > obj.p
                error('NURBSSurface:splitU', ...
                      'Split knot multiplicity exceeds degree.');
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
            U_lo = (U_lo - U_lo(1)) / (U_lo(end) - U_lo(1));

            P_hi = S_ref.P(cp_split:end, :, :);
            W_hi = S_ref.W(cp_split:end, :);
            U_hi = [u0, U2(first:end)];
            U_hi = (U_hi - U_hi(1)) / (U_hi(end) - U_hi(1));

            S_lo = geom.NURBSSurface(P_lo, obj.p, obj.q, U_lo, obj.V, W_lo);
            S_hi = geom.NURBSSurface(P_hi, obj.p, obj.q, U_hi, obj.V, W_hi);
        end

        function [S_lo, S_hi] = splitV(obj, v0)
            St = geom.NURBSSurface(permute(obj.P, [2 1 3]), obj.q, obj.p, obj.V, obj.U, obj.W.');
            [A, B] = St.splitU(v0);
            S_lo = geom.NURBSSurface(permute(A.P, [2 1 3]), obj.p, obj.q, A.V, A.U, A.W.');
            S_hi = geom.NURBSSurface(permute(B.P, [2 1 3]), obj.p, obj.q, B.V, B.U, B.W.');
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
            W2 = Pw(:, 4);
            P2 = bsxfun(@rdivide, Pw(:, 1:3), W2);
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
            W2 = Pw(:, 4);
            P2 = bsxfun(@rdivide, Pw(:, 1:3), W2);
            C = geom.NURBSCurve(P2, obj.q, obj.V, W2);
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
                u_lines = linspace(obj.U(1), obj.U(end), n_iso + 2);
                u_lines = u_lines(2:end-1);
                v_lines = linspace(obj.V(1), obj.V(end), n_iso + 2);
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
            u = max(obj.U(1), min(obj.U(end), u));
        end

        function v = clampV(obj, v)
            v = max(obj.V(1), min(obj.V(end), v));
        end
    end

    methods (Static, Access = public)
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