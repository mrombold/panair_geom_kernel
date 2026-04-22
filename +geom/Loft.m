classdef Loft
    %LOFT Minimal guide/sampling helpers used by ConicSurface workflows.
    %
    % This trimmed version focuses on the geometry operations currently used
    % by the guide-driven conic surface workflow:
    %
    %   - sample/intersect 3D curves with station planes
    %   - build local construction planes
    %   - project 3D geometry into a local 2D plane
    %   - construct a Liming-style rational quadratic conic from 4 points
    %   - combine XY and XZ planar guides into a sampled/interpolated 3D guide
    %
    % Exact NURBS evaluation, interpolation, fitting, splitting, etc. remain in
    % geom.NURBSCurve / geom.NURBSSurface.

    methods (Static)

        %% -----------------------------------------------------------------
        % Curve / station sampling
        % ------------------------------------------------------------------
        function [u, pt, info] = sampleCurveAtStation(C, stationValue, varargin)
            %SAMPLECURVEATSTATION Intersect a curve with a station plane.
            %
            % [u, pt, info] = sampleCurveAtStation(C, stationValue, ...)
            %
            % By default, stationValue is interpreted along the X axis:
            %   dot([x y z], [1 0 0]) = stationValue
            %
            % Optional name/value pairs:
            %   'Axis'       : 'x', 'y', or 'z' (default 'x')
            %   'Normal'     : explicit plane normal, overrides Axis
            %   'BracketN'   : coarse samples used to bracket intersections
            %   'Tol'        : solve tolerance
            %   'MaxIter'    : max Newton/bisection iterations
            %   'Occurrence' : which crossing to use if multiple exist

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
                    error('Loft:sampleCurveAtStation', ...
                        'Normal vector must be nonzero.');
                end
                nrm = nrm / norm(nrm);
            end

            dom = C.domain;
            us = linspace(dom(1), dom(2), max(20, opt.BracketN));
            pts = C.evaluate(us);
            if size(pts, 2) ~= 3
                error('Loft:sampleCurveAtStation', ...
                    'Curve evaluation must return Nx3 points.');
            end

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

        function nrm = axisToNormal(axisName)
            if isstring(axisName), axisName = char(axisName); end
            switch lower(strtrim(axisName))
                case 'x'
                    nrm = [1 0 0];
                case 'y'
                    nrm = [0 1 0];
                case 'z'
                    nrm = [0 0 1];
                otherwise
                    error('Loft:axisToNormal', ...
                        'Unknown axis "%s".', axisName);
            end
        end

        function idx = axisIndex(axisName)
            if isstring(axisName), axisName = char(axisName); end
            switch lower(strtrim(axisName))
                case 'x'
                    idx = 1;
                case 'y'
                    idx = 2;
                case 'z'
                    idx = 3;
                otherwise
                    error('Loft:axisIndex', ...
                        'Unknown axis "%s".', axisName);
            end
        end

        %% -----------------------------------------------------------------
        % Local plane utilities
        % ------------------------------------------------------------------
        function frame = makePlaneFrame(origin, normal, xHint)
            origin = geom.Loft.ensure3DRow(origin);
            normal = geom.Loft.ensure3DRow(normal);

            if norm(normal) < 1e-14
                error('Loft:makePlaneFrame', ...
                    'Plane normal must be nonzero.');
            end

            if nargin < 3 || isempty(xHint)
                xHint = [1 0 0];
                if abs(dot(normal / norm(normal), xHint)) > 0.95
                    xHint = [0 1 0];
                end
            end

            xHint = geom.Loft.ensure3DRow(xHint);
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

        function frame = makeLocalPlaneFrameFromPoints(varargin)
            %MAKELOCALPLANEFRAMEFROMPOINTS Best-fit local plane through 3D points.
            tol = 1e-10;
            if nargin >= 1 && isnumeric(varargin{end}) && isscalar(varargin{end})
                tol = varargin{end};
                ptsIn = varargin(1:end-1);
            else
                ptsIn = varargin;
            end

            if numel(ptsIn) < 3
                error('Loft:makeLocalPlaneFrameFromPoints', ...
                    'Need at least three points.');
            end

            P = zeros(numel(ptsIn), 3);
            for k = 1:numel(ptsIn)
                P(k,:) = geom.Loft.ensure3DRow(ptsIn{k});
            end

            origin = P(1,:);
            [~,~,V] = svd(P - mean(P,1), 0);
            nhat = V(:,end).';
            if norm(nhat) < tol
                error('Loft:makeLocalPlaneFrameFromPoints', ...
                    'Could not determine a stable local plane normal.');
            end

            for i = 2:size(P,1)-1
                v1 = P(i,:)   - P(1,:);
                v2 = P(i+1,:) - P(1,:);
                c = cross(v1, v2);
                if norm(c) > 100*tol
                    if dot(c, nhat) < 0
                        nhat = -nhat;
                    end
                    break;
                end
            end

            xHint = P(end,:) - P(1,:);
            if norm(xHint) < tol
                xHint = P(2,:) - P(1,:);
            end
            frame = geom.Loft.makePlaneFrame(origin, nhat, xHint);
        end

        function uv = toPlane2D(P, frame)
            P = geom.Loft.ensure3DPoints(P);
            d = P - frame.origin;
            uv = [d * frame.xhat(:), d * frame.yhat(:)];
        end

        %% -----------------------------------------------------------------
        % Liming conic construction
        % ------------------------------------------------------------------
        function [C, meta] = limingConic(P0, P1, T, S, varargin)
            %LIMINGCONIC Build a single rational quadratic conic through 4 points.
            %
            % Inputs:
            %   P0 : start point
            %   P1 : end point
            %   T  : tangent-line intersection point
            %   S  : draw-through / shoulder point
            %
            % Optional name/value pairs:
            %   'ShoulderParameter' : prescribed u in (0,1), or []/'auto'
            %   'Tolerance'         : geometric solve tolerance
            %   'FitTolerance'      : tolerance for shoulder-point fit warning

            pa = inputParser;
            addParameter(pa, 'ShoulderParameter', []);
            addParameter(pa, 'Tolerance', 1e-3);
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
                    error('Loft:limingConic', ...
                        ['ShoulderParameter must lie strictly between 0 and 1, ' ...
                         'or be empty/''auto''.']);
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
                fitTol = max(opt.Tolerance);
            end
            
            if err > fitTol
                warning('Loft:limingConic', ...
                    ['Geometry is not exactly representable by a single rational quadratic ' ...
                     'to the requested tolerance. Error = %.16e, FitTolerance = %.16e.'], ...
                    err, fitTol);
            end
            
            relaxedFitTol = fitTol;
            fitAcceptedWithRelaxation = false;

            meta = struct();
            meta.parameter = us;
            meta.weight = w;
            meta.pointError = err;
            meta.fitToleranceUsed = fitTol;
            meta.relaxedFitToleranceUsed = relaxedFitTol;
            meta.fitAcceptedWithRelaxation = fitAcceptedWithRelaxation;
            meta.fitWarning = '';
            if err > relaxedFitTol
                meta.fitWarning = sprintf(['Input geometry is only approximately represented by a single rational quadratic ' ...
                    '(Error = %.3e, RelaxedFitTolerance = %.3e).'], err, relaxedFitTol);
            end
            meta.solveInfo = solveInfo;
            meta.frame = frame;
        end



        
        function [us, w, info] = solveLimingShoulderParameter(P0, P1, T, S, tol)
            if nargin < 5 || isempty(tol)
                tol = 1e-14;
            end
        
            P0 = P0(:).';
            P1 = P1(:).';
            T  = T(:).';
            S  = S(:).';
        
            if numel(P0) ~= 2 || numel(P1) ~= 2 || numel(T) ~= 2 || numel(S) ~= 2
                error('Loft:solveLimingShoulderParameter', ...
                    'P0, P1, T, and S must all be 2D row vectors.');
            end
        
            % Use:
            %   a = P0 - S
            %   b = P1 - S
            %   c = T  - S
            a = P0 - S;
            b = P1 - S;
            c = T  - S;
        
            % Exact conic condition:
            %   (1-u)^2 * det(a,c) + u^2 * det(b,c) = 0
            A = a(1)*c(2) - a(2)*c(1);
            B = b(1)*c(2) - b(2)*c(1);
        
            if abs(A) < tol || abs(B) < tol || A*B >= 0
                error('Loft:solveLimingShoulderParameter', ...
                    ['No exact interior conic solution exists for this geometry. ' ...
                     'Need A*B < 0.']);
            end
        
            rho = sqrt(-A / B);
            us = rho / (1 + rho);
        
            if ~(isfinite(us) && us > 0 && us < 1)
                error('Loft:solveLimingShoulderParameter', ...
                    'Computed shoulder parameter is invalid: us = %.16g.', us);
            end
        
            B0 = (1-us)^2;
            B1 = 2*us*(1-us);
            B2 = us^2;
        
            % From:
            %   B1*w*c = -(B0*a + B2*b)
            num = -dot(c, B0*a + B2*b);
            den =  B1 * dot(c, c);
        
            if abs(den) < tol
                error('Loft:solveLimingShoulderParameter', ...
                    'Degenerate denominator while solving middle weight.');
            end
        
            w = num / den;
        
            if ~(isfinite(w) && w > 0)
                error('Loft:solveLimingShoulderParameter', ...
                    'Computed middle weight is invalid: w = %.16g.', w);
            end
        
            % Exact residual check in local 2D plane
            Cpt = (B0*P0 + B1*w*T + B2*P1) / (B0 + B1*w + B2);
            err = norm(Cpt - S);
        
            info = struct();
            info.method = 'exact_conic_algebraic';
            info.A = A;
            info.B = B;
            info.rho = rho;
            info.pointError = err;
        end
        
        function [w, info] = solveLimingWeightAtParameter(P0, P1, T, S, u, tol)
            if nargin < 6 || isempty(tol)
                tol = 1e-14;
            end
        
            P0 = P0(:).';
            P1 = P1(:).';
            T  = T(:).';
            S  = S(:).';
        
            if numel(P0) ~= 2 || numel(P1) ~= 2 || numel(T) ~= 2 || numel(S) ~= 2
                error('Loft:solveLimingWeightAtParameter', ...
                    'P0, P1, T, and S must all be 2D row vectors.');
            end
            if ~(isscalar(u) && isfinite(u) && u > 0 && u < 1)
                error('Loft:solveLimingWeightAtParameter', ...
                    'u must be a finite scalar in (0,1).');
            end
        
            a = P0 - S;
            b = P1 - S;
            c = T  - S;
        
            B0 = (1-u)^2;
            B1 = 2*u*(1-u);
            B2 = u^2;
        
            % From:
            %   B1*w*c = -(B0*a + B2*b)
            rhs = -(B0*a + B2*b);
            lhs =  B1*c;
        
            mask = abs(lhs) > tol;
            if ~any(mask)
                error('Loft:solveLimingWeightAtParameter', ...
                    'Could not determine middle weight from supplied geometry.');
            end
        
            cands = rhs(mask) ./ lhs(mask);
        
            % Scalar least-squares estimate
            w = dot(lhs(mask), rhs(mask)) / dot(lhs(mask), lhs(mask));
        
            if ~(isfinite(w) && w > 0)
                error('Loft:solveLimingWeightAtParameter', ...
                    'Computed middle weight is invalid: w = %.16g.', w);
            end
        
            residuals = abs(cands - w);
        
            info = struct();
            info.candidates = cands;
            info.candidateSpread = max(residuals);
            info.usedCoordinates = find(mask);
            info.B0 = B0;
            info.B1 = B1;
            info.B2 = B2;
        end










        %% -----------------------------------------------------------------
        % Combine planar guides into a 3D guide
        % ------------------------------------------------------------------
        function [C3, data] = combinePlanarGuidesTo3D(Cxy, Cxz, varargin)
            %COMBINEPLANARGUIDESTO3D Interpolate a 3D guide from XY and XZ guides.
            %
            % This helper:
            %   1) chooses station samples on a shared axis,
            %   2) intersects both planar guides with those station planes,
            %   3) assembles 3D points from the paired projections,
            %   4) globally interpolates a 3D NURBS curve through those points.
            %
            % Name/value pairs:
            %   'Axis'              : shared station axis ('x','y','z'), default 'x'
            %   'Master'            : 'xy', 'xz', or 'stations', default 'xy'
            %   'NumSamples'        : number of station samples, default 31
            %   'Stations'          : explicit station vector when Master='stations'
            %   'Degree'            : interpolation degree, default 3
            %   'ParameterMethod'   : interpolation parameterization, default 'centripetal'
            %   'OccurrenceXY'      : station-plane occurrence on Cxy, default 1
            %   'OccurrenceXZ'      : station-plane occurrence on Cxz, default 1
            %   'BracketN'          : coarse bracketing count, default 400
            %   'Tol'               : station solve tolerance, default 1e-12
            %   'MaxIter'           : max station solve iterations, default 60
            %   'PlanarityTol'      : warning tolerance for off-plane drift, default 1e-6
            %   'RestrictToOverlap' : clamp to shared axis range, default true

            pa = inputParser;
            addParameter(pa, 'Axis', 'x');
            addParameter(pa, 'Master', 'xy');
            addParameter(pa, 'NumSamples', 31);
            addParameter(pa, 'Stations', []);
            addParameter(pa, 'Degree', 3);
            addParameter(pa, 'ParameterMethod', 'centripetal');
            addParameter(pa, 'OccurrenceXY', 1);
            addParameter(pa, 'OccurrenceXZ', 1);
            addParameter(pa, 'BracketN', 400);
            addParameter(pa, 'Tol', 1e-12);
            addParameter(pa, 'MaxIter', 60);
            addParameter(pa, 'PlanarityTol', 1e-6);
            addParameter(pa, 'RestrictToOverlap', true);
            parse(pa, varargin{:});
            opt = pa.Results;

            master = lower(string(opt.Master));
            axisName = char(lower(string(opt.Axis)));
            if ~ismember(axisName, {'x','y','z'})
                error('Loft:combinePlanarGuidesTo3D', ...
                    'Axis must be ''x'', ''y'', or ''z''.');
            end

            geom.Loft.checkPlanarity(Cxy, 'xy', axisName, opt.PlanarityTol);
            geom.Loft.checkPlanarity(Cxz, 'xz', axisName, opt.PlanarityTol);

            switch master
                case "xy"
                    usMaster = linspace(Cxy.domain(1), Cxy.domain(2), max(2, round(opt.NumSamples)));
                    Pm = Cxy.evaluate(usMaster);
                    stations = Pm(:, geom.Loft.axisIndex(axisName));
                case "xz"
                    usMaster = linspace(Cxz.domain(1), Cxz.domain(2), max(2, round(opt.NumSamples)));
                    Pm = Cxz.evaluate(usMaster);
                    stations = Pm(:, geom.Loft.axisIndex(axisName));
                case "stations"
                    if isempty(opt.Stations)
                        error('Loft:combinePlanarGuidesTo3D', ...
                            'Stations must be supplied when Master is ''stations''.');
                    end
                    stations = opt.Stations(:);
                otherwise
                    error('Loft:combinePlanarGuidesTo3D', ...
                        'Master must be ''xy'', ''xz'', or ''stations''.');
            end

            if opt.RestrictToOverlap && master ~= "stations"
                overlap = geom.Loft.estimateOverlapRange(Cxy, Cxz, axisName);
                keep = stations >= overlap(1) & stations <= overlap(2);
                stations = stations(keep);
                if numel(stations) < 2
                    error('Loft:combinePlanarGuidesTo3D', ...
                        'Too few station samples remain in the overlapping range.');
                end
            end

            nS = numel(stations);
            P3 = zeros(nS, 3);
            uxy = zeros(nS, 1);
            uxz = zeros(nS, 1);
            infoXY = cell(nS,1);
            infoXZ = cell(nS,1);

            for k = 1:nS
                sk = stations(k);

                [uxy(k), pxy, infoXY{k}] = geom.Loft.sampleCurveAtStation( ...
                    Cxy, sk, 'Axis', axisName, ...
                    'Occurrence', opt.OccurrenceXY, ...
                    'BracketN', opt.BracketN, ...
                    'Tol', opt.Tol, ...
                    'MaxIter', opt.MaxIter);

                [uxz(k), pxz, infoXZ{k}] = geom.Loft.sampleCurveAtStation( ...
                    Cxz, sk, 'Axis', axisName, ...
                    'Occurrence', opt.OccurrenceXZ, ...
                    'BracketN', opt.BracketN, ...
                    'Tol', opt.Tol, ...
                    'MaxIter', opt.MaxIter);

                P3(k,:) = geom.Loft.mergeProjectedPoint(pxy, pxz, axisName, sk);
            end

            [P3u, ia] = uniquetol(P3, max(1e-12, 10*opt.Tol), ...
                'ByRows', true, 'DataScale', 1);
            [ia, order] = sort(ia);
            P3u = P3u(order,:);
            stations = stations(ia);
            uxy = uxy(ia);
            uxz = uxz(ia);
            infoXY = infoXY(ia);
            infoXZ = infoXZ(ia);

            if size(P3u,1) < 2
                error('Loft:combinePlanarGuidesTo3D', ...
                    'Need at least two distinct 3D points to build a curve.');
            end

            p = min(max(1, round(opt.Degree)), size(P3u,1)-1);
            C3 = geom.NURBSCurve.globalInterp(P3u, p, opt.ParameterMethod);

            data = struct();
            data.points = P3u;
            data.stations = stations;
            data.degreeUsed = p;
            data.parameterMethod = opt.ParameterMethod;
            data.master = char(master);
            data.axis = axisName;
            data.uXY = uxy;
            data.uXZ = uxz;
            data.infoXY = infoXY;
            data.infoXZ = infoXZ;
            data.boundingBox = [min(P3u,[],1); max(P3u,[],1)];
        end

        function checkPlanarity(C, planeTag, stationAxis, tol)
            pts = C.evaluate(linspace(C.domain(1), C.domain(2), 101));
            planeTag = lower(string(planeTag));
            stationAxis = lower(string(stationAxis));

            switch planeTag
                case "xy"
                    off = max(abs(pts(:,3) - mean(pts(:,3))));
                    if off > tol
                        warning('Loft:combinePlanarGuidesTo3D:PlanarityXY', ...
                            'Cxy is not perfectly planar in XY: max |z-zmean| = %.3e.', off);
                    end
                    if stationAxis == "z"
                        warning('Loft:combinePlanarGuidesTo3D:AxisChoice', ...
                            'Using Axis=''z'' with an XY guide is usually not meaningful.');
                    end
                case "xz"
                    off = max(abs(pts(:,2) - mean(pts(:,2))));
                    if off > tol
                        warning('Loft:combinePlanarGuidesTo3D:PlanarityXZ', ...
                            'Cxz is not perfectly planar in XZ: max |y-ymean| = %.3e.', off);
                    end
                    if stationAxis == "y"
                        warning('Loft:combinePlanarGuidesTo3D:AxisChoice', ...
                            'Using Axis=''y'' with an XZ guide is usually not meaningful.');
                    end
                otherwise
                    error('Loft:combinePlanarGuidesTo3D', ...
                        'Unknown plane tag "%s".', planeTag);
            end
        end

        function range = estimateOverlapRange(C1, C2, axisName)
            pts1 = C1.evaluate(linspace(C1.domain(1), C1.domain(2), 201));
            pts2 = C2.evaluate(linspace(C2.domain(1), C2.domain(2), 201));
            idx = geom.Loft.axisIndex(axisName);
            r1 = [min(pts1(:,idx)), max(pts1(:,idx))];
            r2 = [min(pts2(:,idx)), max(pts2(:,idx))];
            range = [max(r1(1), r2(1)), min(r1(2), r2(2))];
            if range(2) < range(1)
                error('Loft:combinePlanarGuidesTo3D', ...
                    'The two guide curves do not overlap in %s.', upper(axisName));
            end
        end

        function P = mergeProjectedPoint(pxy, pxz, axisName, stationValue)
            switch lower(axisName)
                case 'x'
                    P = [stationValue, pxy(2), pxz(3)];
                case 'y'
                    P = [pxy(1), stationValue, pxz(3)];
                case 'z'
                    P = [pxy(1), pxz(2), stationValue];
                otherwise
                    error('Loft:combinePlanarGuidesTo3D', ...
                        'Unknown axis "%s".', axisName);
            end
        end

        %% -----------------------------------------------------------------
        % Basic shape/type checks
        % ------------------------------------------------------------------
        function P = ensure3DPoints(P)
            if isempty(P)
                P = zeros(0,3);
                return;
            end
            if size(P,2) ~= 3
                error('Loft:ensure3DPoints', ...
                    'Points must be Nx3.');
            end
        end

        function p = ensure3DRow(p)
            p = p(:).';
            if numel(p) ~= 3
                error('Loft:ensure3DRow', ...
                    'Point/vector must have 3 elements.');
            end
        end
    end
end
