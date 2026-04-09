classdef SurfaceSurfaceIntersect
% SURFACESURFACEINTERSECT  Practical surface-surface intersection utilities.
%
%   This is a first-pass aerodynamic geometry tool aimed at wing/body work.
%   It is not a full CAD-kernel marching solver yet. Instead it:
%     1) samples iso-curves on S1,
%     2) intersects those curves with S2 using geom.Intersect.curveSurface,
%     3) deduplicates and orders the hits,
%     4) optionally refines each hit as a simultaneous surface/surface solve.
%
%   Main entry points:
%     hits = geom.SurfaceSurfaceIntersect.intersect(S1,S2)
%     C    = geom.SurfaceSurfaceIntersect.toCompositeCurve(hits)
%
%   Returned hit struct fields:
%     .pt      [1x3] Cartesian point
%     .u1 .v1  parameters on S1
%     .u2 .v2  parameters on S2
%     .residual
%     .seedType    'isoU' or 'isoV'
%     .seedValue   iso-parameter value on S1

    methods (Static)

        function hits = intersect(S1, S2, varargin)
        % INTERSECT  Sample-based surface-surface intersection.
        %
        %   hits = intersect(S1,S2,'Nu',14,'Nv',12)
        %
        % Parameters:
        %   Nu, Nv       - number of iso-curves taken from S1
        %   TolPoint     - merge tolerance in xyz
        %   TolParam     - duplicate tolerance in parameters
        %   Refine       - logical, refine with simultaneous Newton solve
        %   CurveSeeds   - number of starting seeds passed to curve-surface solver
        %
        % Returns:
        %   hits - struct array with intersection points and parameters

            pa = inputParser;
            addParameter(pa, 'Nu', 16);
            addParameter(pa, 'Nv', 12);
            addParameter(pa, 'TolPoint', 1e-6);
            addParameter(pa, 'TolParam', 1e-6);
            addParameter(pa, 'Refine', true);
            addParameter(pa, 'CurveSeeds', 8);
            parse(pa, varargin{:});
            opts = pa.Results;

            rawHits = struct('pt', {}, 'u1', {}, 'v1', {}, 'u2', {}, 'v2', {}, ...
                             'residual', {}, 'seedType', {}, 'seedValue', {});

            % Interior iso-U curves on S1
            uu = linspace(S1.U(1), S1.U(end), opts.Nu + 2);
            uu = uu(2:end-1);
            for k = 1:numel(uu)
                C = S1.isoCurveV(uu(k));   % u = const, curve in v
                chits = geom.SurfaceSurfaceIntersect.intersectCurveWithSurface( ...
                    C, S2, 'CurveSeeds', opts.CurveSeeds);
                for j = 1:numel(chits)
                    rawHits(end+1).pt = chits(j).pt; %#ok<AGROW>
                    rawHits(end).u1 = uu(k);
                    rawHits(end).v1 = chits(j).uc;
                    rawHits(end).u2 = chits(j).us;
                    rawHits(end).v2 = chits(j).vs;
                    rawHits(end).residual = chits(j).residual;
                    rawHits(end).seedType = 'isoU';
                    rawHits(end).seedValue = uu(k);
                end
            end

            % Interior iso-V curves on S1
            vv = linspace(S1.V(1), S1.V(end), opts.Nv + 2);
            vv = vv(2:end-1);
            for k = 1:numel(vv)
                C = S1.isoCurveU(vv(k));   % v = const, curve in u
                chits = geom.SurfaceSurfaceIntersect.intersectCurveWithSurface( ...
                    C, S2, 'CurveSeeds', opts.CurveSeeds);
                for j = 1:numel(chits)
                    rawHits(end+1).pt = chits(j).pt; %#ok<AGROW>
                    rawHits(end).u1 = chits(j).uc;
                    rawHits(end).v1 = vv(k);
                    rawHits(end).u2 = chits(j).us;
                    rawHits(end).v2 = chits(j).vs;
                    rawHits(end).residual = chits(j).residual;
                    rawHits(end).seedType = 'isoV';
                    rawHits(end).seedValue = vv(k);
                end
            end

            if isempty(rawHits)
                hits = rawHits;
                return;
            end

            % Optional simultaneous refinement
            if opts.Refine
                for k = 1:numel(rawHits)
                    [u1,v1,u2,v2,pt,res] = geom.SurfaceSurfaceIntersect.refinePoint( ...
                        S1, S2, rawHits(k).u1, rawHits(k).v1, rawHits(k).u2, rawHits(k).v2);
                    rawHits(k).u1 = u1;
                    rawHits(k).v1 = v1;
                    rawHits(k).u2 = u2;
                    rawHits(k).v2 = v2;
                    rawHits(k).pt = pt;
                    rawHits(k).residual = res;
                end
            end

            hits = geom.SurfaceSurfaceIntersect.uniqueHits(rawHits, opts.TolPoint, opts.TolParam);
            hits = geom.SurfaceSurfaceIntersect.orderHits(hits);
        end

        function hits = intersectCurveWithSurface(C, S, varargin)
        % INTERSECTCURVEWITHSURFACE  Wrapper around geom.Intersect.curveSurface if present.
            pa = inputParser;
            addParameter(pa, 'CurveSeeds', 8);
            parse(pa, varargin{:});
            opts = pa.Results;

            % Preferred path: use addon Intersect class if present
            if exist('geom.Intersect', 'class') == 8
                try
                    hits = geom.Intersect.curveSurface(C, S, 'NumSeeds', opts.CurveSeeds);
                    if isempty(hits)
                        hits = struct('uc', {}, 'us', {}, 'vs', {}, 'pt', {}, 'residual', {});
                    end
                    return;
                catch
                    % Fall through to local fallback
                end
            end

            hits = geom.SurfaceSurfaceIntersect.curveSurfaceFallback(C, S, opts.CurveSeeds);
        end

        function hits = curveSurfaceFallback(C, S, nSeeds)
        % CURVESURFACEFALLBACK  Simple seed/refine curve-surface intersection fallback.
            domC = C.domain;
            us = linspace(domC(1), domC(2), max(12, 4*nSeeds));
            pts = C.evaluate(us);

            candidates = [];
            for k = 1:numel(us)
                try
                    [u0, v0, proj, d] = geom.Projection.toSurface(S, pts(k,:));
                catch
                    continue;
                end
                if d < 0.1
                    candidates(end+1,:) = [us(k), u0, v0]; %#ok<AGROW>
                end
            end

            hits = struct('uc', {}, 'us', {}, 'vs', {}, 'pt', {}, 'residual', {});
            for k = 1:size(candidates,1)
                [uc, u2, v2, pt, res, ok] = geom.SurfaceSurfaceIntersect.refineCurveSurfacePoint( ...
                    C, S, candidates(k,1), candidates(k,2), candidates(k,3));
                if ok
                    hits(end+1).uc = uc; %#ok<AGROW>
                    hits(end).us = u2;
                    hits(end).vs = v2;
                    hits(end).pt = pt;
                    hits(end).residual = res;
                end
            end
            hits = geom.SurfaceSurfaceIntersect.uniqueCurveSurfaceHits(hits, 1e-6);
        end

        function [uc, us, vs, pt, res, ok] = refineCurveSurfacePoint(C, S, uc, us, vs)
        % REFINECURVESURFACEPOINT  3x3 Newton solve for C(uc) = S(us,vs).
            ok = false;
            maxIter = 30;
            tol = 1e-10;

            for it = 1:maxIter
                Pc = C.evaluate(uc);
                dC = C.derivative(uc, 1);

                Ps = S.evaluate(us, vs);
                [Su, Sv] = S.partialDerivatives(us, vs);

                R = Pc - Ps;
                J = [dC(:), -Su(:), -Sv(:)];

                if rcond(J) < 1e-12
                    break;
                end

                delta = -J \ R(:);

                uc = max(C.domain(1), min(C.domain(2), uc + delta(1)));
                us = max(S.domainU(1), min(S.domainU(2), us + delta(2)));
                vs = max(S.domainV(1), min(S.domainV(2), vs + delta(3)));

                if norm(delta) < tol && norm(R) < 1e-8
                    ok = true;
                    break;
                end
            end

            pt = 0.5 * (C.evaluate(uc) + S.evaluate(us,vs));
            res = norm(C.evaluate(uc) - S.evaluate(us,vs));
            ok = ok || (res < 1e-7);
        end

        function [u1,v1,u2,v2,pt,res] = refinePoint(S1, S2, u1, v1, u2, v2)
        % REFINEPOINT  Simultaneous surface-surface intersection refinement.
        %
        % Solves:
        %   S1(u1,v1) - S2(u2,v2) = 0
        % with a minimum-norm correction in 4 variables.

            maxIter = 30;
            tol = 1e-10;

            for it = 1:maxIter
                P1 = S1.evaluate(u1, v1);
                P2 = S2.evaluate(u2, v2);
                [S1u, S1v] = S1.partialDerivatives(u1, v1);
                [S2u, S2v] = S2.partialDerivatives(u2, v2);

                R = P1 - P2;
                J = [S1u(:), S1v(:), -S2u(:), -S2v(:)];  % 3x4

                % Minimum-norm least-squares step
                delta = -pinv(J) * R(:);

                u1 = max(S1.domainU(1), min(S1.domainU(2), u1 + delta(1)));
                v1 = max(S1.domainV(1), min(S1.domainV(2), v1 + delta(2)));
                u2 = max(S2.domainU(1), min(S2.domainU(2), u2 + delta(3)));
                v2 = max(S2.domainV(1), min(S2.domainV(2), v2 + delta(4)));

                if norm(delta) < tol && norm(R) < 1e-8
                    break;
                end
            end

            P1 = S1.evaluate(u1, v1);
            P2 = S2.evaluate(u2, v2);
            pt = 0.5 * (P1 + P2);
            res = norm(P1 - P2);
        end

        function hits = uniqueHits(rawHits, tolPt, tolParam)
        % UNIQUEHITS  Remove duplicate surface-surface hits.
            keep = true(1, numel(rawHits));
            for i = 1:numel(rawHits)
                if ~keep(i), continue; end
                for j = i+1:numel(rawHits)
                    if ~keep(j), continue; end
                    samePt = norm(rawHits(i).pt - rawHits(j).pt) < tolPt;
                    sameP1 = norm([rawHits(i).u1 - rawHits(j).u1, rawHits(i).v1 - rawHits(j).v1]) < tolParam;
                    sameP2 = norm([rawHits(i).u2 - rawHits(j).u2, rawHits(i).v2 - rawHits(j).v2]) < tolParam;
                    if samePt || (sameP1 && sameP2)
                        if rawHits(j).residual < rawHits(i).residual
                            keep(i) = false;
                            break;
                        else
                            keep(j) = false;
                        end
                    end
                end
            end
            hits = rawHits(keep);
        end

        function hits = orderHits(hits)
        % ORDERHITS  Order hits along the dominant parameter direction on S1.
            if numel(hits) <= 2
                return;
            end

            uv = [[hits.u1].', [hits.v1].'];
            if std(uv(:,1)) >= std(uv(:,2))
                [~, idx] = sort(uv(:,1));
            else
                [~, idx] = sort(uv(:,2));
            end
            hits = hits(idx);
        end

        function CC = toCompositeCurve(hits)
        % TOCOMPOSITECURVE  Convert ordered hits to a smooth interpolated curve.
            if isempty(hits)
                CC = [];
                return;
            end
            pts = reshape([hits.pt], 3, []).';
            if size(pts,1) == 1
                C = geom.NURBSCurve(pts, 1, [0 0 1 1], [1;1]);
                CC = geom.CompositeCurve({C});
                return;
            end

            if exist('geom.LoftedSurface', 'class') == 8 && ismethod('geom.LoftedSurface', 'globalCurveInterp')
                C = geom.LoftedSurface.globalCurveInterp(pts, 3, min(max(size(pts,1),4), size(pts,1)));
                CC = geom.CompositeCurve({C});
            else
                C = geom.NURBSCurve(pts, 1);
                CC = geom.CompositeCurve({C});
            end
        end

        function hits = uniqueCurveSurfaceHits(hits, tolPt)
            keep = true(1, numel(hits));
            for i = 1:numel(hits)
                if ~keep(i), continue; end
                for j = i+1:numel(hits)
                    if ~keep(j), continue; end
                    if norm(hits(i).pt - hits(j).pt) < tolPt
                        if hits(j).residual < hits(i).residual
                            keep(i) = false;
                            break;
                        else
                            keep(j) = false;
                        end
                    end
                end
            end
            hits = hits(keep);
        end
    end
end
