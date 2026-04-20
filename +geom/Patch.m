classdef Patch < handle
    %PATCH Continuity-aware wrapper around geom.NURBSSurface.
    %
    % Focuses on boundary edge extraction and G0/G1 control for:
    %   - symmetry-plane tangency on a patch edge
    %   - tangent continuity between two adjacent patches
    %
    % Edge ids:
    %   'u0' : first control row in u   (u = domainU(1))
    %   'u1' : last  control row in u   (u = domainU(2))
    %   'v0' : first control row in v   (v = domainV(1))
    %   'v1' : last  control row in v   (v = domainV(2))

    properties
        name char
        S    % geom.NURBSSurface
        topoSurface = []  % optional geom.TopologySurface mirror
    end

    methods
        function obj = Patch(name, S)
            if nargin < 2
                error('geom.Patch:Constructor', 'Requires name and geom.NURBSSurface.');
            end
            obj.name = char(name);
            obj.S = S;
            if exist('geom.TopologySurface', 'class') == 8
                obj.topoSurface = geom.TopologySurface(obj.name, obj.S);
            end
        end

        function info = edgeInfo(obj, edgeId)
            info = geom.Patch.edgeInfoFromSurface(obj.S, edgeId);
        end

        function report = enforceSymmetryPlane(obj, edgeId, planePoint, planeNormal, varargin)
            % Enforce mirror-compatible tangency at a symmetry plane.
            %
            % Geometry enforced:
            %   1) edge control row lies on the plane
            %   2) first interior row differs from the edge only along the
            %      plane normal, so the cross-edge derivative is normal to
            %      the plane and the mirrored mate is G1-compatible.
            p = inputParser;
            addParameter(p, 'ProjectEdgeToPlane', true);
            addParameter(p, 'KeepNormalOffset', true);
            addParameter(p, 'NormalSign', []); % optional preferred sign
            parse(p, varargin{:});
            opt = p.Results;

            n = planeNormal(:).';
            nrm = norm(n);
            if nrm <= eps
                error('geom.Patch:enforceSymmetryPlane', 'planeNormal must be nonzero.');
            end
            n = n / nrm;
            p0 = planePoint(:).';

            ei = obj.edgeInfo(edgeId);
            Hedge = ei.edgeH;
            Hadj  = ei.adjH;

            % Project edge onto plane in Euclidean space while preserving weights.
            if opt.ProjectEdgeToPlane
                for k = 1:size(Hedge,1)
                    wk = Hedge(k,4);
                    Pk = Hedge(k,1:3) / wk;
                    Pk = geom.Patch.projectPointToPlane(Pk, p0, n);
                    Hedge(k,1:3) = wk * Pk;
                end
            end

            % Make first interior row differ only along plane normal.
            for k = 1:size(Hadj,1)
                we = Hedge(k,4);
                Pe = Hedge(k,1:3) / we;

                wa = Hadj(k,4);
                Pa = Hadj(k,1:3) / wa;

                alpha = dot(Pa - Pe, n);
                if ~opt.KeepNormalOffset
                    alpha = 0.0;
                end
                if ~isempty(opt.NormalSign) && alpha ~= 0
                    alpha = sign(opt.NormalSign) * abs(alpha);
                end
                Pa_new = Pe + alpha * n;
                Hadj(k,1:3) = wa * Pa_new;
            end

            obj = obj.setEdgeRows(edgeId, Hedge, Hadj);
            report = obj.symmetryReport(edgeId, p0, n);
        end

        function report = symmetryReport(obj, edgeId, planePoint, planeNormal, varargin)
            p = inputParser;
            addParameter(p, 'NumSamples', 11);
            parse(p, varargin{:});
            nSamples = p.Results.NumSamples;

            n = planeNormal(:).';
            n = n / norm(n);
            p0 = planePoint(:).';

            uv = geom.Patch.edgeSampleUV(obj.S, edgeId, nSamples);
            posErr = zeros(nSamples,1);
            tangentPlaneErr = zeros(nSamples,1);
            for i = 1:nSamples
                [Spt, Su, Sv] = obj.S.firstPartials(uv(i,1), uv(i,2));
                posErr(i) = abs(dot(Spt - p0, n));
                dCross = geom.Patch.crossDerivativeForEdge(edgeId, Su, Sv);
                dCross_n = dCross / max(norm(dCross), eps);
                tangentPlaneErr(i) = norm(cross(dCross_n, n));
            end
            report = struct();
            report.edge = char(edgeId);
            report.maxEdgePlaneDistance = max(posErr);
            report.maxCrossDerivativeOffNormal = max(tangentPlaneErr);
        end

        function report = enforceG1WithPatch(obj, thisEdge, otherPatch, otherEdge, varargin)
            % Edit THIS patch so that it is G0/G1-compatible with otherPatch
            % along the shared edge.
            %
            % Uses homogeneous control rows so rational surfaces are handled
            % consistently at the control-net level.
            p = inputParser;
            addParameter(p, 'Lambda', 1.0);
            addParameter(p, 'AverageEdge', false);
            addParameter(p, 'Tolerance', 1e-8);
            parse(p, varargin{:});
            opt = p.Results;

            eiThis  = obj.edgeInfo(thisEdge);
            eiOther = otherPatch.edgeInfo(otherEdge);

            % Reorient the other rows so edge parameter runs the same way.
            [edgeOther, adjOther, reversed] = geom.Patch.alignRowsToReference( ...
                eiThis.edgeXYZ, eiOther.edgeXYZ, eiOther.edgeH, eiOther.adjH);

            edgeThis = eiThis.edgeH;
            adjThis  = eiThis.adjH;

            if opt.AverageEdge
                edgeTarget = 0.5 * (edgeThis + edgeOther);
            else
                edgeTarget = edgeOther;
            end

            % Enforce opposite first interior rows in homogeneous space.
            adjNew = edgeTarget - opt.Lambda * (adjOther - edgeTarget);

            obj = obj.setEdgeRows(thisEdge, edgeTarget, adjNew);

            report = obj.g1ReportWithPatch(thisEdge, otherPatch, otherEdge, ...
                'NumSamples', 11, 'Tolerance', opt.Tolerance);
            report.reversedParameterization = reversed;
        end

        function report = g1ReportWithPatch(obj, thisEdge, otherPatch, otherEdge, varargin)
            p = inputParser;
            addParameter(p, 'NumSamples', 11);
            addParameter(p, 'Tolerance', 1e-8);
            parse(p, varargin{:});
            nSamples = p.Results.NumSamples;
            tol = p.Results.Tolerance;

            uvThis  = geom.Patch.edgeSampleUV(obj.S, thisEdge, nSamples);
            uvOther = geom.Patch.edgeSampleUV(otherPatch.S, otherEdge, nSamples);

            % Detect parameter reversal from endpoint positions.
            A0 = obj.S.evaluate(uvThis(1,1), uvThis(1,2));
            A1 = obj.S.evaluate(uvThis(end,1), uvThis(end,2));
            B0 = otherPatch.S.evaluate(uvOther(1,1), uvOther(1,2));
            B1 = otherPatch.S.evaluate(uvOther(end,1), uvOther(end,2));
            sameErr = norm(A0 - B0) + norm(A1 - B1);
            revErr  = norm(A0 - B1) + norm(A1 - B0);
            if revErr < sameErr
                uvOther = flipud(uvOther);
            end

            posErr = zeros(nSamples,1);
            tanErr = zeros(nSamples,1);
            planeErr = zeros(nSamples,1);
            for i = 1:nSamples
                [PA, SuA, SvA] = obj.S.firstPartials(uvThis(i,1), uvThis(i,2));
                [PB, SuB, SvB] = otherPatch.S.firstPartials(uvOther(i,1), uvOther(i,2));
                posErr(i) = norm(PA - PB);

                tA = geom.Patch.edgeTangentForEdge(thisEdge, SuA, SvA);
                cA = geom.Patch.crossDerivativeForEdge(thisEdge, SuA, SvA);
                tB = geom.Patch.edgeTangentForEdge(otherEdge, SuB, SvB);
                cB = geom.Patch.crossDerivativeForEdge(otherEdge, SuB, SvB);

                if dot(tA, tB) < 0
                    tB = -tB;
                end
                planeErr(i) = norm(cross(cA, cB)) / max(norm(cA)*norm(cB), eps);
                cA = cA / max(norm(cA), eps);
                cB = cB / max(norm(cB), eps);
                tanErr(i) = min(norm(cA - cB), norm(cA + cB));
            end

            report = struct();
            report.edgeA = char(thisEdge);
            report.edgeB = char(otherEdge);
            report.maxPositionError = max(posErr);
            report.maxCrossDirectionMismatch = max(tanErr);
            report.maxCrossPlaneMismatch = max(planeErr);
            report.isG0 = report.maxPositionError <= tol;
            report.isG1 = report.isG0 && report.maxCrossPlaneMismatch <= sqrt(tol);
        end

        function obj = setEdgeRows(obj, edgeId, edgeH, adjH)
            S = obj.S;
            switch lower(strtrim(edgeId))
                case 'u0'
                    S = geom.Patch.writeRowU(S, 1, edgeH);
                    S = geom.Patch.writeRowU(S, 2, adjH);
                case 'u1'
                    S = geom.Patch.writeRowU(S, size(S.P,1), edgeH);
                    S = geom.Patch.writeRowU(S, size(S.P,1)-1, adjH);
                case 'v0'
                    S = geom.Patch.writeRowV(S, 1, edgeH);
                    S = geom.Patch.writeRowV(S, 2, adjH);
                case 'v1'
                    S = geom.Patch.writeRowV(S, size(S.P,2), edgeH);
                    S = geom.Patch.writeRowV(S, size(S.P,2)-1, adjH);
                otherwise
                    error('geom.Patch:setEdgeRows', 'Unknown edge id "%s".', edgeId);
            end
            obj.S = S;
            if ~isempty(obj.topoSurface)
                obj.topoSurface.S = S;
            end
        end
    end

    methods (Static)
        function info = edgeInfoFromSurface(S, edgeId)
            switch lower(strtrim(edgeId))
                case 'u0'
                    edgeH = geom.Patch.readRowU(S, 1);
                    adjH  = geom.Patch.readRowU(S, 2);
                    varying = 'v';
                case 'u1'
                    edgeH = geom.Patch.readRowU(S, size(S.P,1));
                    adjH  = geom.Patch.readRowU(S, size(S.P,1)-1);
                    varying = 'v';
                case 'v0'
                    edgeH = geom.Patch.readRowV(S, 1);
                    adjH  = geom.Patch.readRowV(S, 2);
                    varying = 'u';
                case 'v1'
                    edgeH = geom.Patch.readRowV(S, size(S.P,2));
                    adjH  = geom.Patch.readRowV(S, size(S.P,2)-1);
                    varying = 'u';
                otherwise
                    error('geom.Patch:edgeInfoFromSurface', 'Unknown edge id "%s".', edgeId);
            end
            info = struct();
            info.edge = char(edgeId);
            info.varyingParameter = varying;
            info.edgeH = edgeH;
            info.adjH = adjH;
            info.edgeXYZ = edgeH(:,1:3) ./ edgeH(:,4);
            info.adjXYZ = adjH(:,1:3) ./ adjH(:,4);
        end

        function H = readRowU(S, idx)
            nV = size(S.P,2);
            H = zeros(nV, 4);
            for j = 1:nV
                w = S.W(idx,j);
                p = reshape(S.P(idx,j,:),1,3);
                H(j,:) = [w*p, w];
            end
        end

        function H = readRowV(S, idx)
            nU = size(S.P,1);
            H = zeros(nU, 4);
            for i = 1:nU
                w = S.W(i,idx);
                p = reshape(S.P(i,idx,:),1,3);
                H(i,:) = [w*p, w];
            end
        end

        function S = writeRowU(S, idx, H)
            if size(H,1) ~= size(S.P,2) || size(H,2) ~= 4
                error('geom.Patch:writeRowU', 'Row size mismatch.');
            end
            for j = 1:size(S.P,2)
                w = H(j,4);
                if w <= 0
                    error('geom.Patch:writeRowU', 'Weights must remain positive.');
                end
                S.W(idx,j) = w;
                S.P(idx,j,:) = reshape(H(j,1:3)/w,1,1,3);
            end
        end

        function S = writeRowV(S, idx, H)
            if size(H,1) ~= size(S.P,1) || size(H,2) ~= 4
                error('geom.Patch:writeRowV', 'Row size mismatch.');
            end
            for i = 1:size(S.P,1)
                w = H(i,4);
                if w <= 0
                    error('geom.Patch:writeRowV', 'Weights must remain positive.');
                end
                S.W(i,idx) = w;
                S.P(i,idx,:) = reshape(H(i,1:3)/w,1,1,3);
            end
        end

        function [rowB, adjB, reversed] = alignRowsToReference(refXYZ, testXYZ, rowBH, adjBH)
            sameErr = norm(refXYZ(1,:) - testXYZ(1,:)) + norm(refXYZ(end,:) - testXYZ(end,:));
            revErr  = norm(refXYZ(1,:) - testXYZ(end,:)) + norm(refXYZ(end,:) - testXYZ(1,:));
            reversed = revErr < sameErr;
            if reversed
                rowB = flipud(rowBH);
                adjB = flipud(adjBH);
            else
                rowB = rowBH;
                adjB = adjBH;
            end
        end

        function uv = edgeSampleUV(S, edgeId, n)
            if nargin < 3 || isempty(n)
                n = 11;
            end
            switch lower(strtrim(edgeId))
                case 'u0'
                    u = S.domainU(1) * ones(n,1);
                    v = linspace(S.domainV(1), S.domainV(2), n).';
                case 'u1'
                    u = S.domainU(2) * ones(n,1);
                    v = linspace(S.domainV(1), S.domainV(2), n).';
                case 'v0'
                    u = linspace(S.domainU(1), S.domainU(2), n).';
                    v = S.domainV(1) * ones(n,1);
                case 'v1'
                    u = linspace(S.domainU(1), S.domainU(2), n).';
                    v = S.domainV(2) * ones(n,1);
                otherwise
                    error('geom.Patch:edgeSampleUV', 'Unknown edge id "%s".', edgeId);
            end
            uv = [u, v];
        end

        function d = edgeTangentForEdge(edgeId, Su, Sv)
            switch lower(strtrim(edgeId))
                case {'u0','u1'}
                    d = Sv;
                case {'v0','v1'}
                    d = Su;
                otherwise
                    error('geom.Patch:edgeTangentForEdge', 'Unknown edge id "%s".', edgeId);
            end
        end

        function d = crossDerivativeForEdge(edgeId, Su, Sv)
            switch lower(strtrim(edgeId))
                case {'u0','u1'}
                    d = Su;
                case {'v0','v1'}
                    d = Sv;
                otherwise
                    error('geom.Patch:crossDerivativeForEdge', 'Unknown edge id "%s".', edgeId);
            end
        end

        function Pproj = projectPointToPlane(P, planePoint, planeNormal)
            Pproj = P - dot(P - planePoint, planeNormal) * planeNormal;
        end
    end
end
