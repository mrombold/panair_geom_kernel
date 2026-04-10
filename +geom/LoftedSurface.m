
classdef LoftedSurface < handle
    % LOFTEDSURFACE Surface created by skinning through section curves.
    %
    % This version performs an actual spanwise interpolation solve for each
    % chordwise control-point column rather than simply stacking section
    % control points into the surface net.
    %
    % Reference:
    %   Piegl & Tiller, The NURBS Book, skinning / global interpolation.

    properties
        surface   % geom.NURBSSurface
        sections  % harmonized section curves
        vParams   % section parameters in v
    end

    methods
        function obj = LoftedSurface(curves, varargin)
            if ~iscell(curves)
                error('LoftedSurface: curves must be a cell array of NURBSCurve.');
            end
            nsec = numel(curves);
            if nsec < 2
                error('LoftedSurface: need at least 2 section curves.');
            end

            pa = inputParser;
            addParameter(pa, 'q', min(3, nsec - 1));
            addParameter(pa, 'method', 'chord');
            addParameter(pa, 'vParams', []);
            parse(pa, varargin{:});
            opts = pa.Results;

            q = opts.q;
            if q >= nsec
                q = nsec - 1;
                warning('LoftedSurface: degree q reduced to %d (nsec-1).', q);
            end

            curves = geom.LoftedSurface.harmonizeCurves(curves);
            obj.sections = curves;

            p = curves{1}.p;
            U = curves{1}.U;
            ncp_u = size(curves{1}.P, 1);

            for k = 2:nsec
                if curves{k}.p ~= p
                    error('LoftedSurface: harmonization failed to match section degree.');
                end
                if size(curves{k}.P,1) ~= ncp_u
                    error('LoftedSurface: harmonization failed to match section control-point count.');
                end
                if numel(curves{k}.U) ~= numel(U) || any(abs(curves{k}.U - U) > 1e-12)
                    error('LoftedSurface: harmonization failed to match section knot vectors.');
                end
            end

            if ~isempty(opts.vParams)
                t = opts.vParams(:);
                if numel(t) ~= nsec
                    error('LoftedSurface: vParams must have one value per section.');
                end
                if abs(t(end) - t(1)) < eps
                    t = linspace(0,1,nsec).';
                else
                    t = (t - t(1)) / (t(end) - t(1));
                end
            else
                t = geom.LoftedSurface.sectionParams(curves, opts.method);
            end
            obj.vParams = t;

            V = geom.LoftedSurface.globalKnotVector(t, q);

            % True skinning: for each chordwise control-point column, solve a
            % spanwise interpolation problem for Cartesian control points and
            % weights using the common section parameter values t.
            ncp_v = nsec;
            P_net = zeros(ncp_u, ncp_v, 3);
            W_net = zeros(ncp_u, ncp_v);

            for i = 1:ncp_u
                pts_col = zeros(nsec, 3);
                w_col   = zeros(nsec, 1);
                for k = 1:nsec
                    pts_col(k,:) = curves{k}.P(i,:);
                    w_col(k)     = curves{k}.W(i);
                end

                Pcol = geom.LoftedSurface.globalInterpValues(t, pts_col, q, V);
                Wcol = geom.LoftedSurface.globalInterpValues(t, w_col, q, V);

                P_net(i,:,:) = Pcol;
                W_net(i,:)   = max(Wcol(:).', eps);
            end

            obj.surface = geom.NURBSSurface(P_net, p, q, U, V, W_net);
        end
    end

    methods
        function C = evaluate(obj, u, v)
            C = obj.surface.evaluate(u, v);
        end

        function N = normal(obj, u, v)
            N = obj.surface.normal(u, v);
        end

        function mesh = isoMesh(obj, nu, nv, varargin)
            mesh = obj.surface.isoMesh(nu, nv, varargin{:});
        end

        function plot(obj, nu, nv, varargin)
            if nargin < 2, nu = 40; end
            if nargin < 3, nv = 40; end
            obj.surface.plot(nu, nv, varargin{:});
            title('Lofted Surface');
            hold on;
            for k = 1:numel(obj.sections)
                obj.sections{k}.plot(100, 'ShowCP', false, ...
                    'Color', [0.9 0.3 0.1], 'LineWidth', 1.5);
            end
        end
    end

    methods (Static)
        function curves = harmonizeCurves(curves)
            % Ensure common degree, knot vector, and control-point count.
            %
            % Strategy:
            %   1) elevate all section degrees to the maximum degree
            %   2) refine all sections to the union knot vector (with multiplicity)
            %   3) if control-point counts still differ, resample/refit all to
            %      the maximum control-point count and a common uniform knot vector

            nsec = numel(curves);
            p_max = max(cellfun(@(c) c.p, curves));

            % 1) Degree elevation to common degree.
            for k = 1:nsec
                if curves{k}.p < p_max
                    curves{k} = curves{k}.elevate(p_max - curves{k}.p);
                end
            end

            % 2) Refine to common union knot vector (respect multiplicities).
            U_ref = curves{1}.U;
            for k = 2:nsec
                U_ref = geom.LoftedSurface.knotUnionWithMultiplicity(U_ref, curves{k}.U);
            end

            for k = 1:nsec
                missing = geom.LoftedSurface.knotsToInsert(curves{k}.U, U_ref, p_max);
                if ~isempty(missing)
                    curves{k} = curves{k}.refine(missing);
                end
            end

            % 3) Fallback resample if any control-count mismatch remains.
            ncp_all = cellfun(@(c) size(c.P,1), curves);
            if any(ncp_all ~= ncp_all(1))
                ncp_max = max(ncp_all);
                U_common = geom.BasisFunctions.MakeUniformKnotVector(ncp_max - 1, p_max);
                for k = 1:nsec
                    curves{k} = geom.LoftedSurface.resampleCurve(curves{k}, ncp_max, p_max);
                    if numel(curves{k}.U) ~= numel(U_common) || any(abs(curves{k}.U - U_common) > 1e-12)
                        curves{k} = geom.NURBSCurve(curves{k}.P, p_max, U_common, curves{k}.W);
                    end
                end
            end
        end

        function Uout = knotUnionWithMultiplicity(Ua, Ub)
            % Union two clamped knot vectors while preserving the maximum
            % multiplicity of each distinct knot value.
            tol = 1e-12;
            vals = unique([Ua(:); Ub(:)]).';
            Uout = [];
            for a = vals
                ma = sum(abs(Ua - a) < tol);
                mb = sum(abs(Ub - a) < tol);
                Uout = [Uout, repmat(a, 1, max(ma, mb))]; %#ok<AGROW>
            end
            Uout = sort(Uout);
        end

        function Xin = knotsToInsert(Ucur, Uref, p)
            %#ok<INUSD>
            % Return the multiset of knots that must be inserted into Ucur so
            % that its multiplicities match Uref.
            tol = 1e-12;
            vals = unique(Uref);
            Xin = [];
            for a = vals
                mcur = sum(abs(Ucur - a) < tol);
                mref = sum(abs(Uref - a) < tol);
                if mref > mcur
                    Xin = [Xin, repmat(a, 1, mref - mcur)]; %#ok<AGROW>
                end
            end
            Xin = sort(Xin);
        end

        function C2 = resampleCurve(C, n_new, p_new)
            % Re-fit a curve using approximately arc-length-spaced samples.
            n_samp = max(3*n_new, 80);
            u_samp = C.arcLengthParam(n_samp);
            pts = C.evaluate(u_samp);
            C2 = geom.LoftedSurface.globalCurveInterp(pts, p_new, n_new);
        end

        function C = globalCurveInterp(pts, p, n_cp)
            % Global interpolation / least-squares fit through point data.
            m = size(pts,1);
            if nargin < 3 || isempty(n_cp)
                n_cp = m;
            end

            t = geom.LoftedSurface.dataParams(pts);
            U = geom.LoftedSurface.knotVectorFromParams(t, p, n_cp);

            P = geom.LoftedSurface.globalInterpValues(t, pts, p, U);
            C = geom.NURBSCurve(P, p, U);
        end

        function X = globalInterpValues(t, Y, p, U)
            % Solve for control values X given data Y(t).
            %
            % t : [m x 1] data parameters in [0,1]
            % Y : [m x d] data values
            % p : degree
            % U : knot vector of target curve
            %
            % Returns X : [n_cp x d] control values

            t = t(:);
            if isvector(Y)
                Y = Y(:);
            end
            m = size(Y,1);
            n_cp = numel(U) - p - 1;

            Nmat = zeros(m, n_cp);
            for row = 1:m
                span = geom.BasisFunctions.FindSpan(n_cp - 1, p, t(row), U);
                Nb = geom.BasisFunctions.BasisFuns(span, t(row), p, U);
                for j = 0:p
                    idx = span - p + j;
                    Nmat(row, idx) = Nb(j+1);
                end
            end

            if m == n_cp
                X = Nmat \ Y;
            else
                X = (Nmat' * Nmat) \ (Nmat' * Y);
            end
        end

        function t = dataParams(pts)
            % Chord-length parameterization of point data.
            m = size(pts,1);
            if m < 2
                t = 0;
                return;
            end
            dk = sqrt(sum(diff(pts,1,1).^2, 2));
            tot = sum(dk);
            if tot < eps
                t = linspace(0,1,m).';
            else
                t = [0; cumsum(dk)/tot];
            end
        end

        function U = knotVectorFromParams(t, p, n_cp)
            % Clamped knot vector from data parameters using averaging.
            m = numel(t);
            n = n_cp - 1;

            U = zeros(1, n_cp + p + 1);
            U(end-p:end) = 1;

            if n_cp <= p
                error('LoftedSurface: n_cp must be at least p+1.');
            end

            if n_cp == m
                for j = 1:(n - p)
                    U(j + p + 1) = mean(t(j+1:j+p));
                end
            else
                d = m / (n_cp - p);
                for j = 1:(n - p)
                    jf = floor(j*d);
                    alpha = j*d - jf;
                    i0 = max(1, min(m, jf));
                    i1 = max(1, min(m, jf + 1));
                    U(j + p + 1) = (1 - alpha) * t(i0) + alpha * t(i1);
                end
            end
        end

        function t = sectionParams(curves, method)
            nsec = numel(curves);
            switch lower(method)
                case 'uniform'
                    t = linspace(0, 1, nsec).';

                case {'chord','centripetal'}
                    reps = zeros(nsec, 3);
                    for k = 1:nsec
                        u_s = curves{k}.arcLengthParam(81);
                        pts = curves{k}.evaluate(u_s);

                        % Representative section position: average of sampled
                        % section points. Stable for closed or open sections.
                        reps(k,:) = mean(pts, 1);
                    end

                    dk = sqrt(sum(diff(reps,1,1).^2, 2));
                    if strcmpi(method, 'centripetal')
                        dk = sqrt(max(dk, 0));
                    end
                    tot = sum(dk);

                    if tot < eps
                        t = linspace(0,1,nsec).';
                    else
                        t = [0; cumsum(dk)/tot];
                    end

                otherwise
                    error('LoftedSurface: unknown method "%s".', method);
            end
        end

        function V = globalKnotVector(t, q)
            nsec = numel(t);
            n = nsec - 1;
            m = n + q + 1;
            V = zeros(1, m + 1);
            V(end-q:end) = 1;

            for j = 1:(n - q)
                V(j + q + 1) = mean(t(j+1:j+q));
            end
        end
    end
end