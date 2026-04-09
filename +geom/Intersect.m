classdef Intersect
% INTERSECT  Intersection helpers for NURBS curve/surface workflows.
%
% This file focuses on the first practical milestone for the panair geometry
% kernel: robust-ish curve/surface intersection suitable for wing-body work,
% section cuts, and boundary extraction.
%
% Main entry points:
%   hits = geom.Intersect.curveSurface(C, S)
%   hits = geom.Intersect.curveSurface(C, S, 'Name', value, ...)
%
% Each hit contains:
%   .uCurve      parameter on the curve
%   .uSurface    u parameter on the surface
%   .vSurface    v parameter on the surface
%   .point       3D point
%   .residual    norm(C(u)-S(us,vs))
%
% Strategy:
%   1) coarse sampling to find seeds
%   2) damped Gauss-Newton on residual R = C(u)-S(us,vs)
%   3) deduplicate nearby roots in parameter and physical space

    methods (Static)
        function hits = curveSurface(C, S, varargin)
        % CURVESURFACE  Find intersections between a NURBS curve and surface.
            pa = inputParser;
            addParameter(pa, 'nCurve', 81);
            addParameter(pa, 'nU', 31);
            addParameter(pa, 'nV', 31);
            addParameter(pa, 'tolXYZ', 1e-8);
            addParameter(pa, 'tolParam', 1e-8);
            addParameter(pa, 'maxIter', 40);
            addParameter(pa, 'seedDistanceFactor', 3.0);
            parse(pa, varargin{:});
            opt = pa.Results;

            uc = linspace(C.domain(1), C.domain(2), opt.nCurve);
            us = linspace(S.domainU(1), S.domainU(2), opt.nU);
            vs = linspace(S.domainV(1), S.domainV(2), opt.nV);

            cpts = C.evaluate(uc);
            spts = zeros(opt.nU, opt.nV, 3);
            for i = 1:opt.nU
                for j = 1:opt.nV
                    spts(i,j,:) = S.evaluate(us(i), vs(j));
                end
            end

            % characteristic coarse spacing for candidate acceptance
            dSurf = inf;
            for i = 1:opt.nU-1
                for j = 1:opt.nV-1
                    p00 = squeeze(spts(i,j,:))';
                    p10 = squeeze(spts(i+1,j,:))';
                    p01 = squeeze(spts(i,j+1,:))';
                    dSurf = min(dSurf, norm(p10-p00) + norm(p01-p00));
                end
            end
            if ~isfinite(dSurf)
                dSurf = 1.0;
            end
            seedTol = opt.seedDistanceFactor * dSurf;

            seeds = [];
            for ic = 1:opt.nCurve
                Pc = cpts(ic,:);
                best = inf;
                bestIJ = [1,1];
                for i = 1:opt.nU
                    for j = 1:opt.nV
                        Ps = squeeze(spts(i,j,:))';
                        dij = norm(Pc - Ps);
                        if dij < best
                            best = dij;
                            bestIJ = [i,j];
                        end
                    end
                end
                if best <= seedTol
                    seeds(end+1,:) = [uc(ic), us(bestIJ(1)), vs(bestIJ(2))]; %#ok<AGROW>
                end
            end

            % Also add a few surface-grid guided seeds using curve projection.
            for i = 1:opt.nU
                for j = 1:opt.nV
                    Ps = squeeze(spts(i,j,:))';
                    try
                        uc0 = geom.Projection.toCurve(C, Ps);
                    catch
                        uc0 = geom.Intersect.projectPointToCurve(C, Ps, uc);
                    end
                    if numel(uc0) > 1
                        uc0 = uc0(1);
                    end
                    Pc = C.evaluate(uc0);
                    if norm(Pc - Ps) <= seedTol
                        seeds(end+1,:) = [uc0, us(i), vs(j)]; %#ok<AGROW>
                    end
                end
            end

            if isempty(seeds)
                hits = struct('uCurve', {}, 'uSurface', {}, 'vSurface', {}, ...
                              'point', {}, 'residual', {});
                return;
            end

            roots = struct('uCurve', {}, 'uSurface', {}, 'vSurface', {}, ...
                           'point', {}, 'residual', {});
            for k = 1:size(seeds,1)
                [ok, hit] = geom.Intersect.refineCurveSurface(C, S, seeds(k,:), opt);
                if ok
                    roots = geom.Intersect.appendUnique(roots, hit, opt.tolXYZ, opt.tolParam);
                end
            end

            hits = roots;
        end

        function [ok, hit] = refineCurveSurface(C, S, x0, opt)
        % REFINECURVESURFACE  Damped Gauss-Newton solve for C(uc)=S(us,vs).
            uc = x0(1);
            us = x0(2);
            vs = x0(3);

            ok = false;
            hit = struct('uCurve', NaN, 'uSurface', NaN, 'vSurface', NaN, ...
                         'point', [NaN NaN NaN], 'residual', inf);

            for iter = 1:opt.maxIter
                Pc = C.evaluate(uc);
                Ct = C.derivative(uc, 1);

                SKL = geom.SurfaceEval.derivatives(S, us, vs, 1);
                Ps = SKL{1,1};
                Su = SKL{2,1};
                Sv = SKL{1,2};

                R = Pc - Ps;
                res = norm(R);
                J = [Ct(:), -Su(:), -Sv(:)];

                % Gauss-Newton / minimum-norm LM step
                A = J' * J + 1e-12 * eye(3);
                g = J' * R(:);
                dx = -A \ g;

                alpha = 1.0;
                improved = false;
                for ls = 1:8
                    uc_try = geom.Intersect.clamp(uc + alpha*dx(1), C.domain);
                    us_try = geom.Intersect.clamp(us + alpha*dx(2), S.domainU);
                    vs_try = geom.Intersect.clamp(vs + alpha*dx(3), S.domainV);

                    Pc_try = C.evaluate(uc_try);
                    Ps_try = S.evaluate(us_try, vs_try);
                    res_try = norm(Pc_try - Ps_try);

                    if res_try <= res
                        uc = uc_try;
                        us = us_try;
                        vs = vs_try;
                        improved = true;
                        break;
                    end
                    alpha = 0.5 * alpha;
                end

                if ~improved
                    break;
                end

                if norm(alpha*dx) < opt.tolParam && res < opt.tolXYZ
                    ok = true;
                    break;
                end
            end

            Pc = C.evaluate(uc);
            Ps = S.evaluate(us, vs);
            res = norm(Pc - Ps);
            if res < 10*opt.tolXYZ
                ok = true;
            end

            hit.uCurve   = uc;
            hit.uSurface = us;
            hit.vSurface = vs;
            hit.point    = 0.5 * (Pc + Ps);
            hit.residual = res;
        end
    end

    methods (Static, Access = private)
        function roots = appendUnique(roots, hit, tolXYZ, tolParam)
            for i = 1:numel(roots)
                dp = norm(hit.point - roots(i).point);
                du = abs(hit.uCurve   - roots(i).uCurve);
                ds = abs(hit.uSurface - roots(i).uSurface) + abs(hit.vSurface - roots(i).vSurface);
                if dp < tolXYZ || (du < tolParam && ds < 2*tolParam)
                    if hit.residual < roots(i).residual
                        roots(i) = hit;
                    end
                    return;
                end
            end
            roots(end+1) = hit; %#ok<AGROW>
        end

        function u = projectPointToCurve(C, P, uc)
            pts = C.evaluate(uc);
            d2 = sum((pts - P).^2, 2);
            [~, idx] = min(d2);
            u = uc(idx);
        end

        function val = clamp(val, dom)
            val = max(dom(1), min(dom(2), val));
        end
    end
end
