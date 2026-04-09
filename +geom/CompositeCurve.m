classdef CompositeCurve < handle
% COMPOSITECURVE  A C0-continuous chain of NURBSCurve segments.
%
%   Used to represent closed fuselage cross-sections (assembled from
%   conic arcs), multi-segment leading edges, or any curve that
%   requires different degrees or knot structures in different regions.
%
%   Construction:
%     CC = geom.CompositeCurve(curves)
%     CC = geom.CompositeCurve(curves, 'Closed', true)
%
%   curves  - cell array of geom.NURBSCurve, joined end-to-end.
%             Each curve's start point should equal the previous curve's
%             end point within tolerance.
%
%   Key methods:
%     pt  = CC.evaluate(t)          % t in [0,1] global parameter
%     D   = CC.derivative(t, k)     % k-th derivative
%     T   = CC.tangent(t)           % unit tangent
%     L   = CC.arcLength()          % total arc length
%     t_u = CC.arcLengthParam(n)    % n uniform arc-length samples
%     C   = CC.toNURBS(p)           % merge into single NURBSCurve
%     CC.plot(n_pts, varargin)       % plot composite curve
%
%   Static factories:
%     CC = geom.CompositeCurve.fromConicSection(P_top,P_bot,P_fwd,P_aft,rhos)
%         Assemble a closed conic fuselage cross-section from 4 quadrant arcs.

    properties
        segments    % {1 x nseg} cell of NURBSCurve
        closed      % logical: last segment end = first segment start
    end

    properties (Dependent)
        nseg        % number of segments
        breakParams % [1 x nseg+1] global parameter breaks (0..1)
    end

    % ------------------------------------------------------------------ %
    methods

        function obj = CompositeCurve(curves, varargin)
            pa = inputParser;
            addParameter(pa, 'Closed',    false);
            addParameter(pa, 'Tolerance', 1e-6);
            parse(pa, varargin{:});

            if ~iscell(curves)
                curves = {curves};
            end
            obj.segments = curves(:)';
            obj.closed   = pa.Results.Closed;

            % Optionally warn on gap
            tol = pa.Results.Tolerance;
            for k = 1:numel(curves)-1
                p_end   = curves{k}.evaluate(curves{k}.U(end));
                p_start = curves{k+1}.evaluate(curves{k+1}.U(1));
                gap = norm(p_end - p_start);
                if gap > tol
                    warning('CompositeCurve: gap of %.2e between segments %d and %d.', gap, k, k+1);
                end
            end
            if obj.closed && numel(curves) > 1
                p_end   = curves{end}.evaluate(curves{end}.U(end));
                p_start = curves{1}.evaluate(curves{1}.U(1));
                gap = norm(p_end - p_start);
                if gap > tol
                    warning('CompositeCurve: closure gap of %.2e.', gap);
                end
            end
        end

    end

    % ------------------------------------------------------------------ %
    %  Dependent properties
    % ------------------------------------------------------------------ %
    methods
        function v = get.nseg(obj)
            v = numel(obj.segments);
        end

        function bp = get.breakParams(obj)
            % Arc-length weighted global parameter breaks
            lens = zeros(1, obj.nseg);
            for k = 1:obj.nseg
                lens(k) = obj.segments{k}.arcLength();
            end
            tot = sum(lens);
            if tot < eps, tot = 1; end
            bp  = [0, cumsum(lens)/tot];
        end
    end

    % ------------------------------------------------------------------ %
    %  Evaluation
    % ------------------------------------------------------------------ %
    methods

        function C = evaluate(obj, t)
        % EVALUATE  Point(s) on composite curve.
        %
        %   t  - global parameter in [0,1], scalar or vector

            t  = t(:)';
            C  = zeros(numel(t), 3);
            bp = obj.breakParams;

            for idx = 1:numel(t)
                tk = max(0, min(1, t(idx)));
                % Find segment
                seg = obj.nseg;
                for k = 1:obj.nseg
                    if tk <= bp(k+1) + eps
                        seg = k; break;
                    end
                end
                % Map to local parameter
                t0 = bp(seg); t1 = bp(seg+1);
                if t1 > t0
                    s = (tk - t0)/(t1 - t0);
                else
                    s = 0;
                end
                U = obj.segments{seg}.U;
                u_local = U(1) + s*(U(end)-U(1));
                C(idx,:) = obj.segments{seg}.evaluate(u_local);
            end
        end

        function D = derivative(obj, t, k)
        % DERIVATIVE  k-th derivative at global parameter t.

            if nargin < 3, k = 1; end
            t  = t(:)';
            D  = zeros(numel(t), 3);
            bp = obj.breakParams;

            for idx = 1:numel(t)
                tk  = max(0, min(1, t(idx)));
                seg = obj.nseg;
                for s = 1:obj.nseg
                    if tk <= bp(s+1) + eps, seg = s; break; end
                end
                t0 = bp(seg); t1 = bp(seg+1);
                s_local = (tk-t0)/(t1-t0+eps);
                U = obj.segments{seg}.U;
                u_local = U(1) + s_local*(U(end)-U(1));
                % Chain rule: d/dt = d/du * du/dt,  du/dt = (U_end-U_1)/(t1-t0)
                scale = (U(end)-U(1)) / (t1-t0+eps);
                D(idx,:) = obj.segments{seg}.derivative(u_local, k) * scale^k;
            end
        end

        function T = tangent(obj, t)
            D = obj.derivative(t, 1);
            nrm = sqrt(sum(D.^2, 2));
            nrm(nrm < eps) = 1;
            T = D ./ nrm;
        end

        function L = arcLength(obj)
            L = 0;
            for k = 1:obj.nseg
                L = L + obj.segments{k}.arcLength();
            end
        end

        function t_u = arcLengthParam(obj, n)
        % ARCLENGTHPARAM  n parameter values uniformly spaced in arc length.

            % Sample densely then re-interpolate
            n_dense = max(500, 10*n);
            t_dense = linspace(0, 1, n_dense);
            pts     = obj.evaluate(t_dense);
            d = [0; cumsum(sqrt(sum(diff(pts).^2, 2)))];
            d = d / d(end);
            s_unif  = linspace(0, 1, n);
            t_u     = interp1(d, t_dense, s_unif, 'pchip');
        end

    end  % evaluation

    % ------------------------------------------------------------------ %
    %  Conversion
    % ------------------------------------------------------------------ %
    methods

        function C = toNURBS(obj, p, n_cp)
        % TONURBS  Merge all segments into a single NURBSCurve.
        %
        %   Resamples at arc-length uniform points and fits a global
        %   NURBS curve.  Useful for lofting operations that need a
        %   single consistent parameterization.
        %
        %   p     - output degree (default 3)
        %   n_cp  - output control points (default 60)

            if nargin < 2, p    = 3;  end
            if nargin < 3, n_cp = 60; end

            t_samp = obj.arcLengthParam(max(4*n_cp, 200));
            pts    = obj.evaluate(t_samp);
            C      = geom.LoftedSurface.globalCurveInterp(pts, p, n_cp);
        end

    end

    % ------------------------------------------------------------------ %
    %  Visualization
    % ------------------------------------------------------------------ %
    methods

        function plot(obj, n_pts, varargin)
            if nargin < 2, n_pts = 300; end
            pa = inputParser;
            addParameter(pa, 'Color',    [0.1 0.4 0.9]);
            addParameter(pa, 'LineWidth', 1.5);
            addParameter(pa, 'ShowCP',   false);
            addParameter(pa, 'ShowBreaks', false);
            parse(pa, varargin{:});
            opts = pa.Results;

            hold on;
            t   = linspace(0, 1, n_pts);
            pts = obj.evaluate(t);
            plot3(pts(:,1), pts(:,2), pts(:,3), ...
                  'Color', opts.Color, 'LineWidth', opts.LineWidth);

            if opts.ShowCP
                for k = 1:obj.nseg
                    P = obj.segments{k}.P;
                    plot3(P(:,1), P(:,2), P(:,3), 'o--', ...
                          'Color', [0.6 0.6 0.6], 'MarkerSize', 4);
                end
            end

            if opts.ShowBreaks
                bp   = obj.breakParams;
                bpts = obj.evaluate(bp);
                plot3(bpts(:,1), bpts(:,2), bpts(:,3), 'ks', ...
                      'MarkerSize', 7, 'MarkerFaceColor', 'k');
            end

            axis equal; grid on;
            xlabel('X'); ylabel('Y'); zlabel('Z');
        end

    end

    % ------------------------------------------------------------------ %
    %  Static factories
    % ------------------------------------------------------------------ %
    methods (Static)

        function CC = fromConicSection(P_top, P_bot, P_fwd, P_aft, rhos)
        % FROMCONICSECTION  Closed cross-section from 4 conic quadrant arcs.
        %
        %   CC = fromConicSection(P_top, P_bot, P_fwd, P_aft, rhos)
        %
        %   P_top  - topmost point    (e.g. [x, 0, r_top])
        %   P_bot  - bottommost point (e.g. [x, 0, -r_bot])
        %   P_fwd  - forwardmost (or widest) point (e.g. [x, r_side, 0])
        %   P_aft  - aftmost (or other side) point  (e.g. [x, -r_side, 0])
        %
        %   rhos   - scalar or [4x1] shoulder params per quadrant (default 0.5)
        %
        %   Returns a closed CompositeCurve (4 rational quadratic arcs).
        %
        %   For a circular cross-section:  rhos = cos(pi/4) / (1+cos(pi/4)) ≈ 0.414

            if nargin < 5, rhos = 0.5; end
            if isscalar(rhos), rhos = repmat(rhos, 4, 1); end

            P_top = P_top(:)'; P_bot = P_bot(:)';
            P_fwd = P_fwd(:)'; P_aft = P_aft(:)';

            % Tangent directions: perpendicular to radius at each cardinal point
            % These are in the YZ cross-section plane (X is axial)
            T_top = (P_fwd - P_top); T_top(1) = 0; T_top = T_top/norm(T_top);
            T_fwd = (P_bot - P_fwd); T_fwd(1) = 0; T_fwd = T_fwd/norm(T_fwd);
            T_bot = (P_aft - P_bot); T_bot(1) = 0; T_bot = T_bot/norm(T_bot);
            T_aft = (P_top - P_aft); T_aft(1) = 0; T_aft = T_aft/norm(T_aft);

            arcs{1} = geom.Aircraft.conicArc(P_top, P_fwd, T_top, T_fwd, rhos(1));
            arcs{2} = geom.Aircraft.conicArc(P_fwd, P_bot, T_fwd, T_bot, rhos(2));
            arcs{3} = geom.Aircraft.conicArc(P_bot, P_aft, T_bot, T_aft, rhos(3));
            arcs{4} = geom.Aircraft.conicArc(P_aft, P_top, T_aft, T_top, rhos(4));

            CC = geom.CompositeCurve(arcs, 'Closed', true);
        end

        function CC = circle(center, radius, normal, n_arcs)
        % CIRCLE  Exact NURBS circle as CompositeCurve.
        %
        %   CC = circle(center, radius, normal, n_arcs)
        %
        %   center  - [1x3]
        %   radius  - scalar
        %   normal  - [1x3] axis normal to circle plane (default [0,0,1])
        %   n_arcs  - number of quadrant arcs (default 4, must be 3 or 4)

            if nargin < 3 || isempty(normal), normal = [0 0 1]; end
            if nargin < 4, n_arcs = 4; end

            center = center(:)'; normal = normal(:)'/norm(normal);

            % Build local frame
            ref = [1 0 0];
            if abs(dot(normal, ref)) > 0.9, ref = [0 1 0]; end
            ex  = cross(normal, ref); ex = ex/norm(ex);
            ey  = cross(normal, ex);

            dtheta = 2*pi / n_arcs;
            w1     = cos(dtheta/2);

            arcs = cell(1, n_arcs);
            for k = 1:n_arcs
                th0 = (k-1)*dtheta;
                th1 = k*dtheta;
                thm = (th0+th1)/2;

                P0 = center + radius*(cos(th0)*ex + sin(th0)*ey);
                P2 = center + radius*(cos(th1)*ex + sin(th1)*ey);
                Pm = center + radius*(cos(thm)*ex + sin(thm)*ey);

                % Shoulder point via tangent intersection
                T0 = -sin(th0)*ex + cos(th0)*ey;
                T2 = -sin(th1)*ex + cos(th1)*ey;
                A  = [T0(:), -T2(:)];
                dp = (P2-P0)';
                st = A \ dp;
                P1 = P0 + st(1)*T0;

                arcs{k} = geom.NURBSCurve([P0; P1; P2], 2, [0 0 0 1 1 1], [1; w1; 1]);
            end
            CC = geom.CompositeCurve(arcs, 'Closed', true);
        end

    end  % static factories

end  % classdef CompositeCurve
