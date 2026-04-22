classdef ConicSurface < handle
    %CONICSURFACE Minimal guide-driven conic patch generator.
    %
    % Inputs are four longitudinal guide curves:
    %   UpperGuide    -> P0(s)
    %   LowerGuide    -> P1(s)
    %   TangencyGuide -> T(s)
    %   ShoulderGuide -> S(s)
    %
    % At each station plane, the four guides are sampled with
    % geom.Loft.sampleCurveAtStation(...), then a single Liming conic is built
    % through geom.Loft.limingConic(P0,P1,T,S,...).
    %
    % This class intentionally only focuses on generating an isomesh.

    properties
        UpperGuide
        LowerGuide
        TangencyGuide
        ShoulderGuide

        SweepOrigin (1,3) double = [0 0 0]
        SweepVector (1,3) double = [1 0 0]
        StationRange (1,2) double = [0 1]

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
        EnableSectionCache logical = true
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
            if norm(obj.SweepVector) < 1e-14
                error('ConicSurface:validateSetup', 'SweepVector must be nonzero.');
            end
            if obj.StationRange(2) <= obj.StationRange(1)
                error('ConicSurface:validateSetup', 'StationRange must be increasing.');
            end
            if obj.uDomain(2) <= obj.uDomain(1)
                error('ConicSurface:validateSetup', 'uDomain must be increasing.');
            end
            if obj.vDomain(2) <= obj.vDomain(1)
                error('ConicSurface:validateSetup', 'vDomain must be increasing.');
            end
        end

        function s = stationAtU(obj, u)
            t = (u - obj.uDomain(1)) / diff(obj.uDomain);
            s = obj.StationRange(1) + t * diff(obj.StationRange);
        end

        function u = uAtStation(obj, s)
            t = (s - obj.StationRange(1)) / diff(obj.StationRange);
            u = obj.uDomain(1) + t * diff(obj.uDomain);
        end

        function offset = stationPlaneOffset(obj, s)
            nhat = obj.SweepVector / norm(obj.SweepVector);
            offset = dot(obj.SweepOrigin, nhat) + s;
        end

        function [P0,P1,T,S,info] = sampleGuidesAtStation(obj, s)
            nhat = obj.SweepVector / norm(obj.SweepVector);
            off = obj.stationPlaneOffset(s);

            [u0, P0, info0] = geom.Loft.sampleCurveAtStation( ...
                obj.UpperGuide, off, 'Normal', nhat, ...
                'Occurrence', obj.OccurrenceUpper, ...
                'BracketN', obj.BracketN, 'Tol', obj.Tol, 'MaxIter', obj.MaxIter);

            [u1, P1, info1] = geom.Loft.sampleCurveAtStation( ...
                obj.LowerGuide, off, 'Normal', nhat, ...
                'Occurrence', obj.OccurrenceLower, ...
                'BracketN', obj.BracketN, 'Tol', obj.Tol, 'MaxIter', obj.MaxIter);

            [ut, T, infot] = geom.Loft.sampleCurveAtStation( ...
                obj.TangencyGuide, off, 'Normal', nhat, ...
                'Occurrence', obj.OccurrenceTangency, ...
                'BracketN', obj.BracketN, 'Tol', obj.Tol, 'MaxIter', obj.MaxIter);

            [us, S, infos] = geom.Loft.sampleCurveAtStation( ...
                obj.ShoulderGuide, off, 'Normal', nhat, ...
                'Occurrence', obj.OccurrenceShoulder, ...
                'BracketN', obj.BracketN, 'Tol', obj.Tol, 'MaxIter', obj.MaxIter);

            info = struct('uUpper',u0,'uLower',u1,'uTangency',ut,'uShoulder',us, ...
                          'upper',info0,'lower',info1,'tangency',infot,'shoulder',infos);
        end

        function C = sectionAtStation(obj, s)
            if obj.EnableSectionCache
                key = sprintf('%.15g', round(s / obj.CacheTol) * obj.CacheTol);
                if isKey(obj.sectionCache, key)
                    C = obj.sectionCache(key);
                    return;
                end
            else
                key = '';
            end

            [P0,P1,T,S] = obj.sampleGuidesAtStation(s);
            C = geom.Loft.limingConic(P0, P1, T, S, obj.LimingArgs{:});

            if obj.EnableSectionCache
                obj.sectionCache(key) = C;
            end
        end

        function C = sectionAt(obj, u)
            s = obj.stationAtU(u);
            C = obj.sectionAtStation(s);
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
    end
end
