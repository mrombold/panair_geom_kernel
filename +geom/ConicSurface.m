classdef ConicSurface < handle
    %CONICSURFACE Guide-driven conic patch generator.
    %
    % Supported marching modes:
    %   1) axis-aligned stationing using SweepOrigin / SweepVector
    %   2) spine-normal stationing using a 3D NURBSCurve spine
    %
    % In spine-normal mode, the section plane at parameter s is:
    %   origin = Spine.evaluate(s)
    %   normal = normalize(Spine.derivative(s,1))
    %
    % The plane frame is completed with PlaneXHint through
    % geom.Loft.makePlaneFrame(...).
    %
    % At each station, the guide curves are intersected with the station
    % plane and the resulting points are used to build one Liming conic.

    properties
        UpperGuide
        LowerGuide
        TangencyGuide
        ShoulderGuide

        MarchMode char = 'axis'

        SweepOrigin (1,3) double = [0 0 0]
        SweepVector (1,3) double = [1 0 0]
        StationRange (1,2) double = [0 1]

        Spine = []
        SpineDomain (1,2) double = [0 1]
        PlaneXHint (1,3) double = [0 0 1]

        uDomain (1,2) double = [0 1]
        vDomain (1,2) double = [0 1]

        BracketN (1,1) double = 400
        Tol (1,1) double = 1e-12
        MaxIter (1,1) double = 60
        OccurrenceUpper (1,1) double = 1
        OccurrenceLower (1,1) double = 1
        OccurrenceTangency (1,1) double = 1
        OccurrenceShoulder (1,1) double = 1

        LimingArgs cell = {}

        EnableSectionCache logical = false
        CacheTol (1,1) double = 1e-10
    end

    properties (Access = private)
        sectionCache
    end

    methods
        function obj = ConicSurface(varargin)
            if rem(nargin,2) ~= 0
                error('ConicSurface:Constructor', ...
                    'Arguments must be supplied as name/value pairs.');
            end
            props = properties(obj);
            for k = 1:2:nargin
                name = varargin{k};
                value = varargin{k+1};
                idx = find(strcmpi(char(name), props), 1, 'first');
                if isempty(idx)
                    error('ConicSurface:Constructor', ...
                        'Unknown property "%s".', char(name));
                end
                obj.(props{idx}) = value;
            end
            obj.sectionCache = containers.Map('KeyType','char','ValueType','any');
            obj.validateSetup();
        end

        function validateSetup(obj)
            req = {'UpperGuide','LowerGuide','TangencyGuide','ShoulderGuide'};
            for i = 1:numel(req)
                if isempty(obj.(req{i}))
                    error('ConicSurface:validateSetup', '%s is required.', req{i});
                end
            end

            if obj.uDomain(2) <= obj.uDomain(1)
                error('ConicSurface:validateSetup', 'uDomain must be increasing.');
            end
            if obj.vDomain(2) <= obj.vDomain(1)
                error('ConicSurface:validateSetup', 'vDomain must be increasing.');
            end

            switch lower(strtrim(obj.MarchMode))
                case 'axis'
                    if norm(obj.SweepVector) < 1e-14
                        error('ConicSurface:validateSetup', ...
                            'SweepVector must be nonzero in axis marching mode.');
                    end
                    if obj.StationRange(2) <= obj.StationRange(1)
                        error('ConicSurface:validateSetup', ...
                            'StationRange must be increasing.');
                    end

                case 'spine_normal'
                    if isempty(obj.Spine)
                        error('ConicSurface:validateSetup', ...
                            'Spine is required in spine_normal mode.');
                    end
                    if obj.SpineDomain(2) <= obj.SpineDomain(1)
                        error('ConicSurface:validateSetup', ...
                            'SpineDomain must be increasing.');
                    end
                    if norm(obj.PlaneXHint) < 1e-14
                        error('ConicSurface:validateSetup', ...
                            'PlaneXHint must be nonzero.');
                    end

                otherwise
                    error('ConicSurface:validateSetup', ...
                        'Unsupported MarchMode "%s".', obj.MarchMode);
            end
        end

        function s = stationAtU(obj, u)
            t = (u - obj.uDomain(1)) / diff(obj.uDomain);
            switch lower(strtrim(obj.MarchMode))
                case 'axis'
                    s = obj.StationRange(1) + t * diff(obj.StationRange);
                case 'spine_normal'
                    s = obj.SpineDomain(1) + t * diff(obj.SpineDomain);
                otherwise
                    error('ConicSurface:stationAtU', ...
                        'Unsupported MarchMode "%s".', obj.MarchMode);
            end
        end

        function u = uAtStation(obj, s)
            switch lower(strtrim(obj.MarchMode))
                case 'axis'
                    t = (s - obj.StationRange(1)) / diff(obj.StationRange);
                case 'spine_normal'
                    t = (s - obj.SpineDomain(1)) / diff(obj.SpineDomain);
                otherwise
                    error('ConicSurface:uAtStation', ...
                        'Unsupported MarchMode "%s".', obj.MarchMode);
            end
            u = obj.uDomain(1) + t * diff(obj.uDomain);
        end

        function frame = frameAtStation(obj, s)
            switch lower(strtrim(obj.MarchMode))
                case 'axis'
                    frame = obj.axisFrameAtStation(s);
                case 'spine_normal'
                    frame = obj.spineNormalFrameAtStation(s);
                otherwise
                    error('ConicSurface:frameAtStation', ...
                        'Unsupported MarchMode "%s".', obj.MarchMode);
            end
        end

        function frame = axisFrameAtStation(obj, s)
            nhat = obj.SweepVector / norm(obj.SweepVector);
            origin = obj.SweepOrigin + s * nhat;
            frame = geom.Loft.makePlaneFrame(origin, nhat, obj.PlaneXHint);
        end

        function frame = spineNormalFrameAtStation(obj, s)
            origin = obj.Spine.evaluate(s);
            d1 = obj.Spine.derivative(s, 1);
            if norm(d1) < 1e-14
                error('ConicSurface:spineNormalFrameAtStation', ...
                    'Spine tangent is degenerate at s = %.16g.', s);
            end
            nhat = d1 / norm(d1);
            frame = geom.Loft.makePlaneFrame(origin, nhat, obj.PlaneXHint);
            frame.spineParameter = s;
            frame.spinePoint = origin;
            frame.spineTangent = nhat;
        end

        function offset = planeOffsetFromFrame(~, frame)
            offset = dot(frame.origin, frame.normal);
        end

        function [P0,P1,T,S,info] = sampleGuidesAtStation(obj, s)
            frame = obj.frameAtStation(s);
            offset = obj.planeOffsetFromFrame(frame);

            [u0, P0, info0] = geom.Loft.sampleCurveAtStation( ...
                obj.UpperGuide, offset, 'Normal', frame.normal, ...
                'Occurrence', obj.OccurrenceUpper, ...
                'BracketN', obj.BracketN, 'Tol', obj.Tol, 'MaxIter', obj.MaxIter);

            [u1, P1, info1] = geom.Loft.sampleCurveAtStation( ...
                obj.LowerGuide, offset, 'Normal', frame.normal, ...
                'Occurrence', obj.OccurrenceLower, ...
                'BracketN', obj.BracketN, 'Tol', obj.Tol, 'MaxIter', obj.MaxIter);

            [ut, T, infot] = geom.Loft.sampleCurveAtStation( ...
                obj.TangencyGuide, offset, 'Normal', frame.normal, ...
                'Occurrence', obj.OccurrenceTangency, ...
                'BracketN', obj.BracketN, 'Tol', obj.Tol, 'MaxIter', obj.MaxIter);

            [us, S, infos] = geom.Loft.sampleCurveAtStation( ...
                obj.ShoulderGuide, offset, 'Normal', frame.normal, ...
                'Occurrence', obj.OccurrenceShoulder, ...
                'BracketN', obj.BracketN, 'Tol', obj.Tol, 'MaxIter', obj.MaxIter);

            info = struct();
            info.frame = frame;
            info.offset = offset;
            info.uUpper = u0;
            info.uLower = u1;
            info.uTangency = ut;
            info.uShoulder = us;
            info.upper = info0;
            info.lower = info1;
            info.tangency = infot;
            info.shoulder = infos;
        end

        function [C, meta] = sectionAtStation(obj, s)
            if obj.EnableSectionCache
                key = sprintf('%.15g', round(s / obj.CacheTol) * obj.CacheTol);
                if isKey(obj.sectionCache, key)
                    data = obj.sectionCache(key);
                    C = data.curve;
                    meta = data.meta;
                    return;
                end
            else
                key = '';
            end

            [P0,P1,T,S,guideInfo] = obj.sampleGuidesAtStation(s);
            [C, meta] = geom.Loft.limingConic(P0, P1, T, S, obj.LimingArgs{:});
            meta.station = s;
            meta.guideInfo = guideInfo;

            if obj.EnableSectionCache
                obj.sectionCache(key) = struct('curve', C, 'meta', meta);
            end
        end

        function C = sectionAt(obj, u)
            s = obj.stationAtU(u);
            [C, ~] = obj.sectionAtStation(s);
        end

        function P = evaluate(obj, u, v)
            C = obj.sectionAt(u);
            vv = obj.mapVToCurveParam(C, v);
            P = C.evaluate(vv);
        end

        function [u_iso, v_iso, pts] = isoGrid(obj, nu, nv)
            u_iso = linspace(obj.uDomain(1), obj.uDomain(2), nu);
            v_iso = linspace(obj.vDomain(1), obj.vDomain(2), nv);
            pts = zeros(nu, nv, 3);

            for i = 1:nu
                C = obj.sectionAt(u_iso(i));
                vv = obj.mapVToCurveParam(C, v_iso(:));
                Pi = C.evaluate(vv);
                pts(i,:,:) = reshape(Pi, [1, nv, 3]);
            end
        end

        function mesh = isoMesh(obj, nu, nv)
            [u_iso, v_iso, pts] = obj.isoGrid(nu, nv);

            mesh.X = squeeze(pts(:,:,1));
            mesh.Y = squeeze(pts(:,:,2));
            mesh.Z = squeeze(pts(:,:,3));
            mesh.nu = nu;
            mesh.nv = nv;
            [mesh.u, mesh.v] = meshgrid(u_iso, v_iso);
            mesh.u = mesh.u.';
            mesh.v = mesh.v.';
            mesh.trimMask = true(size(mesh.u));

            nquads = (nu - 1) * (nv - 1);
            conn = zeros(nquads, 4);
            qid = 1;
            for i = 1:nu-1
                for j = 1:nv-1
                    n1 = (i-1) * nv + j;
                    n2 = i * nv + j;
                    n3 = i * nv + j + 1;
                    n4 = (i-1) * nv + j + 1;
                    conn(qid,:) = [n1 n2 n3 n4];
                    qid = qid + 1;
                end
            end
            mesh.connectivity = conn;
            mesh.normals = [];
        end

        function clearCache(obj)
            obj.sectionCache = containers.Map('KeyType','char','ValueType','any');
        end

        function [Sfit, report, sections] = toNURBS(obj, varargin)
            %TONURBS Approximate a procedural ConicSurface with a NURBSSurface.
            %
            %   [Sfit, report, sections] = obj.toNURBS()
            %   [Sfit, report, sections] = obj.toNURBS('Name', value, ...)
            %
            % This routine samples the ConicSurface on a structured station/section grid,
            % interpolates a NURBS section curve at each sampled station, then lofts those
            % section curves into a NURBS surface. A validation report is returned by
            % comparing the fitted NURBS surface against the analytic/procedural
            % ConicSurface on a denser grid.
            %
            % Name/value options
            % ------------------
            % 'Nu'            : number of station sections to build         (default 41)
            % 'Nv'            : samples per section for curve interpolation (default 41)
            % 'SectionDegree' : degree of each fitted section curve         (default 3)
            % 'LoftDegree'    : loft degree in station direction            (default 3)
            % 'ParameterMethod': section interpolation parameterization     (default 'centripetal')
            % 'ValidationNu'  : validation samples in station direction     (default 81)
            % 'ValidationNv'  : validation samples in section direction     (default 81)
            %
            % Outputs
            % -------
            % Sfit     : geom.NURBSSurface approximate fitted surface
            % report   : struct with error statistics and settings
            % sections : cell array of fitted station-section NURBS curves
            
                pa = inputParser;
                addParameter(pa, 'Nu', 41);
                addParameter(pa, 'Nv', 41);
                addParameter(pa, 'SectionDegree', 3);
                addParameter(pa, 'LoftDegree', 3);
                addParameter(pa, 'ParameterMethod', 'centripetal');
                addParameter(pa, 'ValidationNu', 81);
                addParameter(pa, 'ValidationNv', 81);
                parse(pa, varargin{:});
                opt = pa.Results;
            
                nu = max(2, round(opt.Nu));
                nv = max(2, round(opt.Nv));
                pSec = max(1, round(opt.SectionDegree));
                qLoft = max(1, round(opt.LoftDegree));
            
                uIso = linspace(obj.uDomain(1), obj.uDomain(2), nu);
                vIso = linspace(obj.vDomain(1), obj.vDomain(2), nv);
            
                sections = cell(1, nu);
                secPts = cell(1, nu);
            
                for i = 1:nu
                    P = zeros(nv, 3);
                    for j = 1:nv
                        P(j,:) = obj.evaluate(uIso(i), vIso(j));
                    end
                    secPts{i} = P;
            
                    pUse = min(pSec, size(P,1) - 1);
                    sections{i} = geom.NURBSCurve.globalInterp(P, pUse, opt.ParameterMethod);
                end
            
                qUse = min(qLoft, numel(sections) - 1);
                Sfit = geom.NURBSSurface.loft(sections, qUse);
            
                % Validation against procedural surface
                uVal = linspace(obj.uDomain(1), obj.uDomain(2), max(2, round(opt.ValidationNu)));
                vVal = linspace(obj.vDomain(1), obj.vDomain(2), max(2, round(opt.ValidationNv)));
            
                errs = zeros(numel(uVal), numel(vVal));
                for i = 1:numel(uVal)
                    for j = 1:numel(vVal)
                        Ptrue = obj.evaluate(uVal(i), vVal(j));
                        Pfit  = Sfit.evaluate(uVal(i), vVal(j));
                        errs(i,j) = norm(Pfit - Ptrue);
                    end
                end
            
                report = struct();
                report.nu = nu;
                report.nv = nv;
                report.sectionDegree = pUse;
                report.loftDegree = qUse;
                report.parameterMethod = opt.ParameterMethod;
                report.validationNu = numel(uVal);
                report.validationNv = numel(vVal);
                report.maxError = max(errs(:));
                report.rmsError = sqrt(mean(errs(:).^2));
                report.meanError = mean(errs(:));
                report.errorGrid = errs;
                report.uValidation = uVal;
                report.vValidation = vVal;
                report.sectionPoints = secPts;
            end

    end

    methods (Access = private)
        function vv = mapVToCurveParam(obj, C, v)
            t = (v(:) - obj.vDomain(1)) / diff(obj.vDomain);
            t = max(0, min(1, t));
            vv = C.domain(1) + t * diff(C.domain);
        end
    end

    methods (Static)
        function obj = fromGuides(varargin)
            obj = geom.ConicSurface(varargin{:});
        end

        function obj = fromSpineNormal(varargin)
            args = [{'MarchMode','spine_normal'}, varargin];
            obj = geom.ConicSurface(args{:});
        end
    end
end
