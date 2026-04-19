
classdef Loft
    %LOFT Higher-level geometric construction helpers built on top of
    % geom.NURBSCurve and geom.NURBSSurface.
    %
    % This class is intended to hold lofting / construction techniques that
    % are not core generic NURBS primitives.  The methods return ordinary
    % geom.NURBSCurve / geom.NURBSSurface objects and plain MATLAB structs.
    %
    % Initial capabilities:
    %   - sample / intersect 3D curves at a station plane
    %   - simple local plane-frame utilities
    %   - line/line intersection in 2D
    %   - Liming-style quadratic rational conic from:
    %       start point, end point, tangent intersection, shoulder point
    %   - fuselage section construction from longitudinal guide curves
    %   - section loft wrapper
    %
    % Notes:
    %   - Liming conic uses a rational quadratic Bezier / NURBS segment.
    %   - The middle weight and shoulder parameter can be solved directly
    %     from the geometry (start, end, tangent intersection, shoulder).

    methods (Static)
        function [u, pt, info] = sampleCurveAtStation(C, stationValue, varargin)
            pa = inputParser;
            addParameter(pa, 'Axis', 'x');
            addParameter(pa, 'Normal', []);
            addParameter(pa, 'BracketN', 400);
            addParameter(pa, 'Tol', 1e-12);
            addParameter(pa, 'MaxIter', 60);
            addParameter(pa, 'Occurrence', 1);
            parse(pa, varargin{:});
            opt = pa.Results;

            nrm = geom.Loft.axisToNormal(opt.Axis);
            if ~isempty(opt.Normal)
                nrm = opt.Normal(:).';
                if norm(nrm) < eps
                    error('Loft:sampleCurveAtStation', 'Normal vector must be nonzero.');
                end
                nrm = nrm / norm(nrm);
            end

            dom = C.domain;
            us = linspace(dom(1), dom(2), max(20, opt.BracketN));
            pts = C.evaluate(us);
            f = pts * nrm(:) - stationValue;

            [fminAbs, idx0] = min(abs(f));
            if fminAbs < opt.Tol
                u = us(idx0);
                pt = pts(idx0,:);
                info = struct('normal', nrm, 'usedBracket', false, ...
                    'occurrence', 1, 'residual', f(idx0), 'iterations', 0);
                return;
            end

            brackets = zeros(0,2);
            for k = 1:(numel(us)-1)
                if f(k) == 0
                    brackets(end+1,:) = [us(k), us(k)]; %#ok<AGROW>
                elseif f(k) * f(k+1) < 0
                    brackets(end+1,:) = [us(k), us(k+1)]; %#ok<AGROW>
                end
            end

            if isempty(brackets)
                error('Loft:sampleCurveAtStation', ...
                    ['Could not bracket intersection with station plane. ' ...
                     'Check that the curve crosses the requested station value.']);
            end

            occ = opt.Occurrence;
            if occ < 1 || occ ~= floor(occ) || occ > size(brackets,1)
                error('Loft:sampleCurveAtStation', ...
                    'Requested occurrence %d is invalid; only %d bracket(s) found.', ...
                    occ, size(brackets,1));
            end

            a = brackets(occ,1);
            b = brackets(occ,2);
            if a == b
                u = a;
                pt = C.evaluate(u);
                info = struct('normal', nrm, 'usedBracket', true, ...
                    'occurrence', occ, 'residual', dot(nrm, pt) - stationValue, ...
                    'iterations', 0);
                return;
            end

            fa = dot(C.evaluate(a), nrm) - stationValue;
            fb = dot(C.evaluate(b), nrm) - stationValue;
            u = 0.5 * (a + b);
            iters = 0;

            for iter = 1:opt.MaxIter
                iters = iter;
                p = C.evaluate(u);
                d1 = C.derivative(u, 1);
                fu = dot(p, nrm) - stationValue;
                df = dot(d1, nrm);

                if abs(fu) < opt.Tol
                    break;
                end

                usedNewton = false;
                if abs(df) > 1e-14
                    un = u - fu / df;
                    if un > a && un < b
                        uTrial = un;
                        usedNewton = true;
                    end
                end
                if ~usedNewton
                    uTrial = 0.5 * (a + b);
                end

                pTrial = C.evaluate(uTrial);
                fTrial = dot(pTrial, nrm) - stationValue;

                if fa == 0
                    u = a;
                    break;
                elseif fb == 0
                    u = b;
                    break;
                elseif fa * fTrial <= 0
                    b = uTrial;
                    fb = fTrial;
                else
                    a = uTrial;
                    fa = fTrial;
                end

                uNew = 0.5 * (a + b);
                if abs(uNew - u) < opt.Tol
                    u = uNew;
                    break;
                end
                u = uNew;
            end

            pt = C.evaluate(u);
            info = struct('normal', nrm, 'usedBracket', true, ...
                'occurrence', occ, 'residual', dot(nrm, pt) - stationValue, ...
                'iterations', iters);
        end

        function [u, pt, info] = intersectCurveWithPlane(C, planeNormal, planeOffset, varargin)
            [u, pt, info] = geom.Loft.sampleCurveAtStation(C, planeOffset, ...
                'Normal', planeNormal, varargin{:});
        end

        function frame = makePlaneFrame(origin, normal, xHint)
            origin = origin(:).';
            normal = normal(:).';
            if nargin < 3 || isempty(xHint)
                xHint = [1 0 0];
                if abs(dot(normal / norm(normal), xHint)) > 0.95
                    xHint = [0 1 0];
                end
            end
            xHint = xHint(:).';

            nhat = normal / norm(normal);
            xhat = xHint - dot(xHint, nhat) * nhat;
            if norm(xhat) < 1e-12
                alt = [1 0 0];
                if abs(dot(alt, nhat)) > 0.95
                    alt = [0 1 0];
                end
                xhat = alt - dot(alt, nhat) * nhat;
            end
            xhat = xhat / norm(xhat);
            yhat = cross(nhat, xhat);
            yhat = yhat / norm(yhat);

            frame = struct();
            frame.origin = origin;
            frame.normal = nhat;
            frame.xhat = xhat;
            frame.yhat = yhat;
            frame.R = [xhat(:), yhat(:), nhat(:)];
        end

        function uv = toPlane2D(P, frame)
            P = geom.Loft.ensure3DPoints(P);
            d = P - frame.origin;
            uv = [d * frame.xhat(:), d * frame.yhat(:)];
        end

        function P = fromPlane2D(uv, frame)
            if size(uv,2) ~= 2
                error('Loft:fromPlane2D', 'uv must be Nx2.');
            end
            P = frame.origin + uv(:,1) * frame.xhat + uv(:,2) * frame.yhat;
        end

        function [pInt, t1, t2] = lineIntersection2D(P1, d1, P2, d2, tol)
            if nargin < 5 || isempty(tol), tol = 1e-12; end
            P1 = P1(:).';
            P2 = P2(:).';
            d1 = d1(:).';
            d2 = d2(:).';
            A = [d1(:), -d2(:)];
            rhs = (P2 - P1).';
            if abs(det(A)) < tol
                error('Loft:lineIntersection2D', '2D lines are parallel or nearly parallel.');
            end
            t = A \ rhs;
            t1 = t(1);
            t2 = t(2);
            pInt = P1 + t1 * d1;
        end

        function [C, meta] = limingConic(P0, P1, T, S, varargin)
            pa = inputParser;
            addParameter(pa, 'ShoulderParameter', []);
            addParameter(pa, 'Tolerance', 1e-8);
            addParameter(pa, 'FitTolerance', []);
            parse(pa, varargin{:});
            opt = pa.Results;

            P0 = geom.Loft.ensure3DRow(P0);
            P1 = geom.Loft.ensure3DRow(P1);
            T  = geom.Loft.ensure3DRow(T);
            S  = geom.Loft.ensure3DRow(S);

            frame = geom.Loft.makeLocalPlaneFrameFromPoints(P0, P1, T, S, opt.Tolerance);
            p02 = geom.Loft.toPlane2D(P0, frame);
            p12 = geom.Loft.toPlane2D(P1, frame);
            t2  = geom.Loft.toPlane2D(T,  frame);
            s2  = geom.Loft.toPlane2D(S,  frame);

            if isempty(opt.ShoulderParameter) || ...
               (ischar(opt.ShoulderParameter) && strcmpi(strtrim(opt.ShoulderParameter), 'auto')) || ...
               (isstring(opt.ShoulderParameter) && strcmpi(strtrim(opt.ShoulderParameter), "auto"))
                [us, w, solveInfo] = geom.Loft.solveLimingShoulderParameter(p02, p12, t2, s2, opt.Tolerance);
            else
                us = opt.ShoulderParameter;
                if ~isscalar(us) || us <= 0 || us >= 1
                    error('Loft:limingConic', 'ShoulderParameter must lie strictly between 0 and 1, or be empty/''auto''.');
                end
                [w, solveInfo] = geom.Loft.solveLimingWeightAtParameter(p02, p12, t2, s2, us, opt.Tolerance);
            end

            P = [P0; T; P1];
            U = [0 0 0 1 1 1];
            W = [1; w; 1];
            C = geom.NURBSCurve(P, 2, U, W);

            pChk = C.evaluate(us);
            err = norm(pChk - S);
            scale = max([1, norm(P1-P0), norm(T-P0), norm(S-P0)]);
            fitTol = opt.FitTolerance;
            if isempty(fitTol)
                fitTol = max(1e-5 * scale, 100 * opt.Tolerance);
            end
            if err > fitTol
                relaxedFitTol = max(10 * fitTol, 1e-4 * scale);
                if err > relaxedFitTol
                    error('Loft:limingConic', ...
                        ['Input geometry is not well represented by a single rational quadratic ' ...
                         'for the solved shoulder parameter. Error = %.3e, FitTolerance = %.3e'], err, relaxedFitTol);
                end
                fitAcceptedWithRelaxation = true;
            else
                fitAcceptedWithRelaxation = false;
            end

            meta = struct();
            meta.parameter = us;
            meta.weight = w;
            meta.pointError = err;
            meta.fitToleranceUsed = fitTol;
            meta.fitAcceptedWithRelaxation = fitAcceptedWithRelaxation;
            meta.solveInfo = solveInfo;
            meta.frame = frame;
        end

        function sec = buildFuselageSectionAtStation(stationValue, varargin)
            pa = inputParser;
            addParameter(pa, 'UpperProfile', []);
            addParameter(pa, 'LowerProfile', []);
            addParameter(pa, 'MaxBreadth', []);
            addParameter(pa, 'UpperShoulder', []);
            addParameter(pa, 'LowerShoulder', []);
            addParameter(pa, 'Axis', 'x');
            addParameter(pa, 'PlaneNormal', []);
            addParameter(pa, 'BreadthDirection', []);
            addParameter(pa, 'VerticalDirection', []);
            addParameter(pa, 'ShoulderParameter', []);
            addParameter(pa, 'Method', 'liming');
            parse(pa, varargin{:});
            opt = pa.Results;

            if isempty(opt.UpperProfile) || isempty(opt.LowerProfile) || ...
               isempty(opt.MaxBreadth) || isempty(opt.UpperShoulder) || ...
               isempty(opt.LowerShoulder)
                error('Loft:buildFuselageSectionAtStation', ...
                    'UpperProfile, LowerProfile, MaxBreadth, UpperShoulder, and LowerShoulder are required.');
            end

            nrm = geom.Loft.axisToNormal(opt.Axis);
            if ~isempty(opt.PlaneNormal)
                nrm = opt.PlaneNormal(:).';
                nrm = nrm / norm(nrm);
            end

            [~, Ptop] = geom.Loft.sampleCurveAtStation(opt.UpperProfile, stationValue, ...
                'Normal', nrm);
            [~, Pbot] = geom.Loft.sampleCurveAtStation(opt.LowerProfile, stationValue, ...
                'Normal', nrm);
            [~, Pmax] = geom.Loft.sampleCurveAtStation(opt.MaxBreadth, stationValue, ...
                'Normal', nrm);
            [~, Sup]  = geom.Loft.sampleCurveAtStation(opt.UpperShoulder, stationValue, ...
                'Normal', nrm);
            [~, Slo]  = geom.Loft.sampleCurveAtStation(opt.LowerShoulder, stationValue, ...
                'Normal', nrm);

            origin = Ptop;
            if strcmpi(strtrim(opt.Axis), 'x')
                xHint = [0 1 0];
            elseif strcmpi(strtrim(opt.Axis), 'y')
                xHint = [1 0 0];
            else
                xHint = [1 0 0];
            end
            frame = geom.Loft.makePlaneFrame(origin, nrm, xHint);

            pTop2 = geom.Loft.toPlane2D(Ptop, frame);
            pBot2 = geom.Loft.toPlane2D(Pbot, frame);
            pMax2 = geom.Loft.toPlane2D(Pmax, frame);
            sUp2  = geom.Loft.toPlane2D(Sup,  frame);
            sLo2  = geom.Loft.toPlane2D(Slo,  frame);

            if isempty(opt.BreadthDirection)
                dBreadth2 = [1 0];
            else
                tmp = geom.Loft.toPlane2D(frame.origin + opt.BreadthDirection(:).', frame) ...
                    - geom.Loft.toPlane2D(frame.origin, frame);
                dBreadth2 = tmp / norm(tmp);
            end
            if isempty(opt.VerticalDirection)
                dVert2 = [0 1];
            else
                tmp = geom.Loft.toPlane2D(frame.origin + opt.VerticalDirection(:).', frame) ...
                    - geom.Loft.toPlane2D(frame.origin, frame);
                dVert2 = tmp / norm(tmp);
            end

            [tUp2, ~, ~] = geom.Loft.lineIntersection2D(pTop2, dVert2, pMax2, dBreadth2);
            [tLo2, ~, ~] = geom.Loft.lineIntersection2D(pBot2, dVert2, pMax2, dBreadth2);

            Tup = geom.Loft.fromPlane2D(tUp2, frame);
            Tlo = geom.Loft.fromPlane2D(tLo2, frame);

            switch lower(strtrim(opt.Method))
                case 'liming'
                    [Cup, upperMeta] = geom.Loft.limingConic(Ptop, Pmax, Tup, Sup, ...
                        'ShoulderParameter', opt.ShoulderParameter);
                    [Clo, lowerMeta] = geom.Loft.limingConic(Pmax, Pbot, Tlo, Slo, ...
                        'ShoulderParameter', opt.ShoulderParameter);
                otherwise
                    error('Loft:buildFuselageSectionAtStation', ...
                        'Unknown section construction method: %s', opt.Method);
            end

            sec = struct();
            sec.stationValue = stationValue;
            sec.frame = frame;
            sec.points = struct('top', Ptop, 'bottom', Pbot, 'maxBreadth', Pmax, ...
                'upperShoulder', Sup, 'lowerShoulder', Slo);
            sec.points2D = struct('top', pTop2, 'bottom', pBot2, 'maxBreadth', pMax2, ...
                'upperShoulder', sUp2, 'lowerShoulder', sLo2);
            sec.tangentIntersectionUpper = Tup;
            sec.tangentIntersectionLower = Tlo;
            sec.upperCurve = Cup;
            sec.lowerCurve = Clo;
            if exist('upperMeta','var'), sec.upperMeta = upperMeta; end
            if exist('lowerMeta','var'), sec.lowerMeta = lowerMeta; end
        end

        function sections = buildSectionFamily(stations, varargin)
            stations = stations(:).';
            sections = cell(size(stations));
            for k = 1:numel(stations)
                sections{k} = geom.Loft.buildFuselageSectionAtStation(stations(k), varargin{:});
            end
        end

        function S = loftSections(curves, q, method, sectionParams)
            if nargin < 2 || isempty(q), q = 3; end
            if nargin < 3 || isempty(method), method = 'centripetal'; end
            if nargin < 4
                S = geom.NURBSSurface.loft(curves, q, method);
            else
                S = geom.NURBSSurface.loft(curves, q, method, sectionParams);
            end
        end
    end

    methods (Static, Access = private)
        function [us, w, info] = solveLimingShoulderParameter(P0, P1, T, S, tol)
            if nargin < 5 || isempty(tol), tol = 1e-10; end

            P0 = P0(:).';
            P1 = P1(:).';
            T  = T(:).';
            S  = S(:).';

            uMin = max(1e-6, 10*tol);
            uMax = 1 - uMin;
            uGrid = linspace(uMin, uMax, 2001);
            best = struct('res', inf, 'u', NaN, 'w', NaN, 'info', struct());

            for k = 1:numel(uGrid)
                uk = uGrid(k);
                try
                    [wk, wInfo] = geom.Loft.solveLimingWeightAtParameter(P0, P1, T, S, uk, tol);
                catch
                    continue;
                end
                if isfinite(wk) && wk > 0 && wInfo.pointResidual < best.res
                    best.res = wInfo.pointResidual;
                    best.u = uk;
                    best.w = wk;
                    best.info = wInfo;
                end
            end

            if ~(isfinite(best.u))
                error('Loft:solveLimingShoulderParameter', ...
                    'Could not find a valid shoulder parameter in (0,1) with positive middle weight.');
            end

            % Local refinement around the best coarse point.
            du = (uMax - uMin) / (numel(uGrid) - 1);
            uLo = max(uMin, best.u - 2*du);
            uHi = min(uMax, best.u + 2*du);
            uFine = linspace(uLo, uHi, 2001);
            for k = 1:numel(uFine)
                uk = uFine(k);
                try
                    [wk, wInfo] = geom.Loft.solveLimingWeightAtParameter(P0, P1, T, S, uk, tol);
                catch
                    continue;
                end
                if isfinite(wk) && wk > 0 && wInfo.pointResidual < best.res
                    best.res = wInfo.pointResidual;
                    best.u = uk;
                    best.w = wk;
                    best.info = wInfo;
                end
            end

            us = best.u;
            w = best.w;
            info = best.info;
            info.searchResidual = best.res;
            info.searchGridBest = best.u;
        end

        function [w, info] = solveLimingWeightAtParameter(P0, P1, T, S, us, tol)
            if nargin < 6 || isempty(tol), tol = 1e-10; end

            P0 = P0(:).';
            P1 = P1(:).';
            T  = T(:).';
            S  = S(:).';

            B0 = (1 - us)^2;
            B1 = 2 * us * (1 - us);
            B2 = us^2;

            D = T - S;
            R = B0 * (S - P0) + B2 * (S - P1);

            if norm(D) < tol
                error('Loft:solveLimingWeightAtParameter', ...
                    'Tangent intersection T must not coincide with shoulder point S.');
            end
            if abs(B1) < tol
                error('Loft:solveLimingWeightAtParameter', ...
                    'Shoulder parameter is too close to an endpoint.');
            end

            w = dot(R, D) / (B1 * dot(D, D));
            resVec = B1 * w * D - R;
            ptResidual = norm(resVec);

            if ~(isfinite(w))
                error('Loft:solveLimingWeightAtParameter', ...
                    'Computed middle weight is not finite.');
            end

            info = struct();
            info.pointResidual = ptResidual;
            info.R = R;
            info.D = D;
        end

        function frame = makeLocalPlaneFrameFromPoints(P0, P1, T, S, tol)
            if nargin < 5 || isempty(tol), tol = 1e-10; end
            P0 = P0(:).'; P1 = P1(:).'; T = T(:).'; S = S(:).';
            n = cross(P1 - P0, T - P0);
            if norm(n) < tol, n = cross(P1 - P0, S - P0); end
            if norm(n) < tol, n = cross(T - P0, S - P0); end
            if norm(n) < tol
                error('Loft:makeLocalPlaneFrameFromPoints', ...
                    'Could not determine a stable local plane from P0, P1, T, and S.');
            end
            n = n / norm(n);
            xHint = P1 - P0;
            if norm(xHint) < tol, xHint = [1 0 0]; end
            frame = geom.Loft.makePlaneFrame(P0, n, xHint);
        end

        function nrm = axisToNormal(axisSpec)
            if isnumeric(axisSpec) && numel(axisSpec) == 3
                nrm = axisSpec(:).';
                nrm = nrm / norm(nrm);
                return;
            end
            switch lower(strtrim(char(axisSpec)))
                case 'x'
                    nrm = [1 0 0];
                case 'y'
                    nrm = [0 1 0];
                case 'z'
                    nrm = [0 0 1];
                otherwise
                    error('Loft:axisToNormal', 'Unknown axis specifier: %s', char(axisSpec));
            end
        end

        function P = ensure3DPoints(P)
            if size(P,2) == 2
                P = [P, zeros(size(P,1),1)];
            elseif size(P,2) ~= 3
                error('Loft:ensure3DPoints', 'Points must be Nx2 or Nx3.');
            end
        end

        function p = ensure3DRow(p)
            p = p(:).';
            if numel(p) == 2
                p = [p, 0];
            elseif numel(p) ~= 3
                error('Loft:ensure3DRow', 'Point must have 2 or 3 coordinates.');
            end
        end
    end
end
