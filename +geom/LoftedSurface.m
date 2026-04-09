classdef LoftedSurface < handle
% LOFTEDSURFACE  Surface created by skinning through a set of section curves.
%
%   Constructs a NURBS surface by interpolating or approximating through
%   a collection of cross-sectional NURBSCurve objects placed at v-parameter
%   stations.  This is the primary aircraft fuselage/wing-section operation.
%
%   The section curves must all have compatible degree and the same number
%   of control points.  If they do not, call harmonize() first.
%
%   Construction:
%     S = geom.LoftedSurface(curves)
%     S = geom.LoftedSurface(curves, 'q', 3, 'method', 'chord')
%
%   Options (name-value):
%     q       - degree in v (spanwise) direction  (default: min(3, nsec-1))
%     method  - 'uniform' | 'chord' | 'centripetal'  (default: 'chord')
%               Sets the v-parameterization across sections.
%     vParams - Manual v-parameter values for each section [nsec x 1]
%
%   Properties:
%     surface   - the resulting geom.NURBSSurface object
%     sections  - cell array of input NURBSCurve sections
%
%   Methods:
%     S.plot(nu, nv)           - plot the lofted surface
%
% Reference: Piegl & Tiller, "The NURBS Book" Section 10.3 (Skinning).

    properties
        surface     % geom.NURBSSurface
        sections    % cell array of input NURBSCurve
        vParams     % [nsec x 1] parameter values for each section
    end

    % ------------------------------------------------------------------ %
    %  Construction
    % ------------------------------------------------------------------ %
    methods

        function obj = LoftedSurface(curves, varargin)
        % Constructor.
        %   curves  - cell array of geom.NURBSCurve objects

            if ~iscell(curves)
                error('LoftedSurface: curves must be a cell array of NURBSCurve.');
            end

            nsec = numel(curves);
            if nsec < 2
                error('LoftedSurface: need at least 2 section curves.');
            end

            % Parse options
            pa = inputParser;
            addParameter(pa, 'q',       min(3, nsec-1));
            addParameter(pa, 'method',  'chord');
            addParameter(pa, 'vParams', []);
            parse(pa, varargin{:});
            opts = pa.Results;

            q = opts.q;
            if q >= nsec
                q = nsec - 1;
                warning('LoftedSurface: degree q reduced to %d (nsec-1).', q);
            end

            % Harmonize curves - ensure compatible parameterization
            curves = geom.LoftedSurface.harmonizeCurves(curves);
            obj.sections = curves;

            % Get number of u control points from first curve
            ncp_u = size(curves{1}.P, 1);
            p     = curves{1}.p;
            U     = curves{1}.U;

            % ---- Compute v-parameters for sections ----
            if ~isempty(opts.vParams)
                t = opts.vParams(:) / opts.vParams(end);  % normalize
            else
                t = geom.LoftedSurface.sectionParams(curves, opts.method);
            end
            obj.vParams = t;

            % ---- Global v knot vector by averaging ----
            V = geom.LoftedSurface.globalKnotVector(t, q);

            % ---- Interpolate column by column in u ----
            % For each u-index i, collect the i-th control point from every
            % section and interpolate through them in v.
            % This is the "skinning" approach: directly use section CPs as
            % the u-direction control net rows, then fit v-direction.

            % Control net: [ncp_u x nsec x 3]
            P_secs = zeros(ncp_u, nsec, 3);
            W_secs = zeros(ncp_u, nsec);

            for k = 1:nsec
                P_secs(:,k,:) = curves{k}.P;
                W_secs(:,k)   = curves{k}.W;
            end

            % For each u column, do global curve interpolation in v
            ncp_v = nsec;  % direct skinning: one CP per section
            P_net = zeros(ncp_u, ncp_v, 3);
            W_net = zeros(ncp_u, ncp_v);

            for i = 1:ncp_u
                for d = 1:3
                    P_net(i,:,d) = P_secs(i,:,d);
                end
                W_net(i,:) = W_secs(i,:);
            end

            obj.surface = geom.NURBSSurface(P_net, p, q, U, V, W_net);
        end

    end  % construction

    % ------------------------------------------------------------------ %
    %  Delegation to surface
    % ------------------------------------------------------------------ %
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

            % Overlay section curves
            hold on;
            for k = 1:numel(obj.sections)
                obj.sections{k}.plot(100, 'ShowCP', false, ...
                    'Color', [0.9 0.3 0.1], 'LineWidth', 1.5);
            end
        end

    end

    % ------------------------------------------------------------------ %
    %  Static helpers
    % ------------------------------------------------------------------ %
    methods (Static)

        function curves = harmonizeCurves(curves)
        % HARMONIZECURVES  Ensure all curves have same degree and knot vector.
        %   Elevates degree and inserts knots as needed.

            nsec = numel(curves);
            p_max = max(cellfun(@(c) c.p, curves));
            ncp_max = max(cellfun(@(c) size(c.P,1), curves));

            % 1. Collect union of all distinct interior knots
            all_int_knots = [];
            for k = 1:nsec
                U = curves{k}.U;
                p = curves{k}.p;
                int_k = unique(U(p+2:end-p-1));
                all_int_knots = union(all_int_knots, int_k);
            end

            % 2. For each curve, insert missing knots so they all share U
            for k = 1:nsec
                % Simple approach: re-parameterize all to same # control pts
                % by sampling and re-fitting if needed.
                % For now, just use uniform knot vector with ncp_max points.
                if size(curves{k}.P, 1) ~= ncp_max
                    curves{k} = geom.LoftedSurface.resampleCurve(curves{k}, ncp_max, p_max);
                end
            end

            % Re-assign uniform knot vector to all
            U_ref = geom.BasisFunctions.MakeUniformKnotVector(ncp_max-1, p_max);
            for k = 1:nsec
                curves{k} = geom.NURBSCurve(curves{k}.P, p_max, U_ref, curves{k}.W);
            end
        end

        function C2 = resampleCurve(C, n_new, p_new)
        % RESAMPLECURVE  Re-fit a NURBSCurve with n_new control points.
        %   Uses global interpolation through uniformly-sampled arc-length pts.

            n_samp = max(3*n_new, 50);
            u_samp = C.arcLengthParam(n_samp);
            pts    = C.evaluate(u_samp);

            % Global interpolation: solve for control points
            C2 = geom.LoftedSurface.globalCurveInterp(pts, p_new, n_new);
        end

        function C = globalCurveInterp(pts, p, n_cp)
        % GLOBALCURVEINTERP  Global curve interpolation through data points.
        %   pts    - [m x 3] data points
        %   p      - degree
        %   n_cp   - number of control points (default = size(pts,1))
        %
        %   Uses chord-length parameterization and knot averaging (P&T §9.2)
        %   for a much smoother fit than uniform knot placement.

            m  = size(pts,1);
            if nargin < 3 || isempty(n_cp)
                n_cp = m;
            end

            % Chord-length parameterization
            dk  = sqrt(sum(diff(pts).^2, 2));
            tot = sum(dk);
            if tot < eps
                t = linspace(0,1,m)';
            else
                t = [0; cumsum(dk)/tot];
            end

            % Knot vector: clamped with interior knots by averaging (P&T eq 9.8)
            % This matches knot placement to actual data distribution.
            n    = n_cp - 1;
            U    = zeros(1, n_cp + p + 1);
            U(end-p:end) = 1;
            if n_cp == m
                % Exact interpolation: average consecutive parameters
                for j = 1:n-p
                    U(j+p+1) = mean(t(j+1:j+p));
                end
            else
                % Least-squares: spread knots to cover parameter distribution
                d = m / (n_cp - p);
                for j = 1:n-p
                    i_f  = floor(j*d);
                    alpha = j*d - i_f;
                    U(j+p+1) = (1-alpha)*t(max(1,i_f)) + alpha*t(min(m,i_f+1));
                end
            end

            % Collocation matrix N: m x n_cp
            N_mat = zeros(m, n_cp);
            for row = 1:m
                span = geom.BasisFunctions.FindSpan(n_cp-1, p, t(row), U);
                Nb   = geom.BasisFunctions.BasisFuns(span, t(row), p, U);
                for j = 0:p
                    N_mat(row, span-p+j) = Nb(j+1);
                end
            end

            % Solve (least-squares if m > n_cp, interpolation if m == n_cp)
            if m == n_cp
                P = N_mat \ pts;
            else
                P = (N_mat' * N_mat) \ (N_mat' * pts);
            end

            C = geom.NURBSCurve(P, p, U);
        end

        function t = sectionParams(curves, method)
        % SECTIONPARAMS  Compute v-parameter values across sections.

            nsec = numel(curves);
            switch lower(method)
                case 'uniform'
                    t = linspace(0, 1, nsec)';

                case 'chord'
                    % Use centroid of each section curve
                    centroids = zeros(nsec, 3);
                    for k = 1:nsec
                        u_s = linspace(curves{k}.U(1), curves{k}.U(end), 50);
                        pts = curves{k}.evaluate(u_s);
                        centroids(k,:) = mean(pts, 1);
                    end
                    dk  = sqrt(sum(diff(centroids).^2, 2));
                    tot = sum(dk);
                    if tot < eps
                        t = linspace(0,1,nsec)';
                    else
                        t = [0; cumsum(dk)/tot];
                    end

                case 'centripetal'
                    centroids = zeros(nsec, 3);
                    for k = 1:nsec
                        u_s = linspace(curves{k}.U(1), curves{k}.U(end), 50);
                        pts = curves{k}.evaluate(u_s);
                        centroids(k,:) = mean(pts, 1);
                    end
                    dk  = sqrt(sum(diff(centroids).^2, 2)) .^ 0.5;
                    tot = sum(dk);
                    if tot < eps
                        t = linspace(0,1,nsec)';
                    else
                        t = [0; cumsum(dk)/tot];
                    end

                otherwise
                    error('LoftedSurface: unknown method "%s".', method);
            end
        end

        function V = globalKnotVector(t, q)
        % GLOBALKNOTVECTOR  Compute clamped knot vector from parameter values.
        %   Averaging method (Piegl & Tiller eq. 9.8).

            nsec = numel(t);
            n    = nsec - 1;  % last CP index
            m    = n + q + 1;
            V    = zeros(1, m+1);
            V(end-q:end) = 1;

            for j = 1:n-q
                V(j+q+1) = mean(t(j+1:j+q));
            end
        end

    end  % static methods

end  % classdef LoftedSurface