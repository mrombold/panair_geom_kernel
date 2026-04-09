classdef Projection
% PROJECTION  Closest-point projection onto NURBS curves and surfaces.
%
%   Newton-Raphson iteration for finding the parameter value(s) that
%   minimize distance from a query point to a curve or surface.
%   Used for wing-body intersection, point inversion, and mesh projection.
%
%   Static methods:
%     [u, pt, d]  = Projection.toCurve(C, P, u0)
%     [u, v, pt, d] = Projection.toSurface(S, P, u0, v0)
%     [u, pts, ds]  = Projection.toCurveBatch(C, pts)
%     [u, v, pts, ds] = Projection.toSurfaceBatch(S, pts)

    properties (Constant)
        MAX_ITER = 50
        TOL_U    = 1e-10   % parameter convergence tolerance
        TOL_D    = 1e-10   % distance convergence tolerance
    end

    methods (Static)

        function [u_out, pt_out, d_out] = toCurve(C, P, u0)
        % TOCURVE  Closest point on NURBSCurve to query point P.
        %
        %   [u, pt, d] = toCurve(C, P)
        %   [u, pt, d] = toCurve(C, P, u0)   % with initial guess
        %
        %   C   - geom.NURBSCurve
        %   P   - [1x3] query point
        %   u0  - initial parameter guess (default: coarse search)
        %
        %   Returns:
        %   u   - parameter of closest point
        %   pt  - [1x3] closest point coordinates
        %   d   - distance from P to pt

            P  = P(:)';
            dom = C.domain;

            % Initial guess via coarse search if not provided
            if nargin < 3 || isempty(u0)
                u0 = geom.Projection.coarseSearchCurve(C, P, 32);
            end
            u0 = max(dom(1), min(dom(2), u0));

            % Newton-Raphson:  minimize f(u) = |C(u) - P|^2
            % f'(u) = 2*(C(u)-P) . C'(u) = 0
            % Newton step: du = -f'(u)/f''(u)
            % f'(u)  = (C-P).C'
            % f''(u) = C'.C' + (C-P).C''

            u = u0;
            for iter = 1:geom.Projection.MAX_ITER
                Ck  = C.evaluate(u);
                D1  = C.derivative(u, 1);
                D2  = C.derivative(u, 2);

                r    = Ck - P;
                f1   = dot(r, D1);
                f2   = dot(D1,D1) + dot(r,D2);

                if abs(f2) < eps, break; end
                du = -f1/f2;

                % Clamp to domain
                u_new = max(dom(1), min(dom(2), u + du));
                du    = u_new - u;
                u     = u_new;

                if abs(du) < geom.Projection.TOL_U && norm(r) < geom.Projection.TOL_D
                    break;
                end
            end

            pt_out = C.evaluate(u);
            d_out  = norm(pt_out - P);
            u_out  = u;
        end

        function [u_out, v_out, pt_out, d_out] = toSurface(S, P, u0, v0)
        % TOSURFACE  Closest point on NURBSSurface to query point P.
        %
        %   [u, v, pt, d] = toSurface(S, P)
        %   [u, v, pt, d] = toSurface(S, P, u0, v0)
        %
        %   S   - geom.NURBSSurface
        %   P   - [1x3] query point
        %   u0, v0  - initial parameter guesses (default: coarse search)
        %
        %   Uses 2D Newton-Raphson on the system:
        %     F1(u,v) = (S(u,v) - P) . Su = 0
        %     F2(u,v) = (S(u,v) - P) . Sv = 0

            P = P(:)';

            if nargin < 3 || isempty(u0) || isempty(v0)
                [u0, v0] = geom.Projection.coarseSearchSurface(S, P, 12, 12);
            end
            u0 = max(S.domainU(1), min(S.domainU(2), u0));
            v0 = max(S.domainV(1), min(S.domainV(2), v0));

            u = u0; v = v0;
            du_step = 1e-6*(S.domainU(2)-S.domainU(1));
            dv_step = 1e-6*(S.domainV(2)-S.domainV(1));

            for iter = 1:geom.Projection.MAX_ITER
                Sk       = S.evaluate(u, v);
                [Su, Sv] = S.partialDerivatives(u, v);

                r  = Sk - P;
                F1 = dot(r, Su);
                F2 = dot(r, Sv);

                % Jacobian (numerical second derivatives)
                [Su_u,~] = S.partialDerivatives(u+du_step, v);
                [Su_m,~] = S.partialDerivatives(u-du_step, v);
                [~,Sv_v] = S.partialDerivatives(u, v+dv_step);
                [~,Sv_m] = S.partialDerivatives(u, v-dv_step);
                [Su_v,~] = S.partialDerivatives(u, v+dv_step);
                [Su_mv,~]= S.partialDerivatives(u, v-dv_step);

                Suu = (Su_u - Su_m)/(2*du_step);
                Svv = (Sv_v - Sv_m)/(2*dv_step);
                Suv = (Su_v - Su_mv)/(2*dv_step);

                J11 = dot(Su,Su) + dot(r,Suu);
                J12 = dot(Su,Sv) + dot(r,Suv);
                J21 = J12;
                J22 = dot(Sv,Sv) + dot(r,Svv);

                J   = [J11 J12; J21 J22];
                F   = [F1; F2];

                if rcond(J) < eps, break; end
                delta = -J \ F;

                u_new = max(S.domainU(1), min(S.domainU(2), u + delta(1)));
                v_new = max(S.domainV(1), min(S.domainV(2), v + delta(2)));

                conv = norm([u_new-u, v_new-v]);
                u = u_new; v = v_new;

                if conv < geom.Projection.TOL_U
                    break;
                end
            end

            pt_out = S.evaluate(u, v);
            d_out  = norm(pt_out - P);
            u_out  = u;
            v_out  = v;
        end

        function [u_all, pt_all, d_all] = toCurveBatch(C, pts)
        % TOCURVEBATCH  Closest-point projection for many query points.
        %
        %   pts  - [N x 3] query points
        %   Returns u_all [Nx1], pt_all [Nx3], d_all [Nx1]

            N     = size(pts, 1);
            u_all = zeros(N, 1);
            pt_all= zeros(N, 3);
            d_all = zeros(N, 1);

            % Coarse search: one pass for all points to get good initial guesses
            n_coarse = 64;
            u_coarse = linspace(C.U(1), C.U(end), n_coarse);
            pts_c    = C.evaluate(u_coarse);

            for k = 1:N
                P = pts(k,:);
                d2 = sum(bsxfun(@minus, pts_c, P).^2, 2);
                [~, idx] = min(d2);
                u0 = u_coarse(idx);

                [u_all(k), pt_all(k,:), d_all(k)] = ...
                    geom.Projection.toCurve(C, P, u0);
            end
        end

        function [u_all, v_all, pt_all, d_all] = toSurfaceBatch(S, pts, nu_c, nv_c)
        % TOSURFACEBATCH  Closest-point projection for many points onto surface.
        %
        %   pts  - [N x 3] query points
        %   nu_c, nv_c - coarse grid resolution for initial search (default 16x16)

            if nargin < 3, nu_c = 16; end
            if nargin < 4, nv_c = 16; end

            N     = size(pts, 1);
            u_all = zeros(N, 1);
            v_all = zeros(N, 1);
            pt_all= zeros(N, 3);
            d_all = zeros(N, 1);

            % Build coarse surface grid once
            [u_c, v_c, pts_c] = S.isoGrid(nu_c, nv_c);
            pts_flat = reshape(pts_c, [], 3);
            [ug, vg] = meshgrid(u_c, v_c);
            ug = ug'; vg = vg';
            u_flat = ug(:); v_flat = vg(:);

            for k = 1:N
                P  = pts(k,:);
                d2 = sum(bsxfun(@minus, pts_flat, P).^2, 2);
                [~, idx] = min(d2);
                u0 = u_flat(idx); v0 = v_flat(idx);

                [u_all(k), v_all(k), pt_all(k,:), d_all(k)] = ...
                    geom.Projection.toSurface(S, P, u0, v0);
            end
        end

        function [u_int, v_int, pt_int] = curveSurfaceIntersect(C, S, n_seed)
        % CURVESURFACEINTERSECT  Find intersection(s) of curve with surface.
        %
        %   Marches along the curve and detects sign changes in the
        %   signed distance  d(t) = N(u,v) . (C(t) - S(u,v))
        %   then refines with Newton iteration.
        %
        %   C  - geom.NURBSCurve (the piercing curve)
        %   S  - geom.NURBSSurface
        %   n_seed - number of seeds for initial search (default 200)
        %
        %   Returns arrays of intersection parameters u_int, v_int and points.

            if nargin < 3, n_seed = 200; end

            t_seed = linspace(C.U(1), C.U(end), n_seed);
            pts_c  = C.evaluate(t_seed);

            % For each curve point, find closest surface point
            [~, ~, surf_pts, ~] = geom.Projection.toSurfaceBatch(S, pts_c, 10, 10);

            % Signed distance: (C(t) - S_closest) . N_surface
            sgn = zeros(n_seed, 1);
            for k = 1:n_seed
                u_k = t_seed(k);
                % Re-project (use cached result here for efficiency)
                r = pts_c(k,:) - surf_pts(k,:);
                sgn(k) = norm(r);  % unsigned - sign change detection below
            end

            % Detect crossings by projecting and checking side
            % This is a simplified version - full implementation would
            % use bisection or Newton on the combined system.
            warning('Projection.curveSurfaceIntersect: basic seed-refine only. See README.');

            u_int = []; v_int = []; pt_int = zeros(0,3);

            for k = 1:n_seed-1
                P1 = pts_c(k,:);   P2 = pts_c(k+1,:);
                [~, ~, S1, d1] = geom.Projection.toSurface(S, P1);
                [~, ~, S2, d2] = geom.Projection.toSurface(S, P2);

                n1 = S.normal( ...
                    geom.Projection.toSurface(S, P1), ...
                    0.5);  % simplified

                side1 = dot(P1 - S1, S.normal(0.5, 0.5));
                side2 = dot(P2 - S2, S.normal(0.5, 0.5));

                if sign(side1) ~= sign(side2)
                    % Bisect to find crossing
                    ta = t_seed(k); tb = t_seed(k+1);
                    for bi = 1:20
                        tm = (ta+tb)/2;
                        Pm = C.evaluate(tm);
                        [u_m, v_m, Sm, ~] = geom.Projection.toSurface(S, Pm);
                        Nm = S.normal(u_m, v_m);
                        side_m = dot(Pm - Sm, Nm);
                        if sign(side_m) == sign(side1), ta = tm;
                        else, tb = tm; end
                    end
                    t_int = (ta+tb)/2;
                    Pint  = C.evaluate(t_int);
                    [u_s, v_s, ~, ~] = geom.Projection.toSurface(S, Pint);
                    u_int(end+1) = u_s; %#ok<AGROW>
                    v_int(end+1) = v_s; %#ok<AGROW>
                    pt_int(end+1,:) = S.evaluate(u_s, v_s); %#ok<AGROW>
                end
            end
        end

    end  % static methods

    % ------------------------------------------------------------------ %
    %  Private coarse search helpers
    % ------------------------------------------------------------------ %
    methods (Static, Access = private)

        function u0 = coarseSearchCurve(C, P, n)
            u_c  = linspace(C.U(1), C.U(end), n);
            pts  = C.evaluate(u_c);
            d2   = sum(bsxfun(@minus, pts, P).^2, 2);
            [~,i]= min(d2);
            u0   = u_c(i);
        end

        function [u0, v0] = coarseSearchSurface(S, P, nu, nv)
            u_c  = linspace(S.domainU(1), S.domainU(2), nu);
            v_c  = linspace(S.domainV(1), S.domainV(2), nv);
            best_d2 = inf;
            u0 = u_c(1); v0 = v_c(1);
            for i = 1:nu
                for j = 1:nv
                    pt = S.evaluate(u_c(i), v_c(j));
                    d2 = sum((pt-P).^2);
                    if d2 < best_d2
                        best_d2 = d2;
                        u0 = u_c(i); v0 = v_c(j);
                    end
                end
            end
        end

    end

end  % classdef Projection
