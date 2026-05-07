classdef SurfaceTools
%SURFACETOOLS  CATIA-style surface utilities for geom.NURBSSurface.
%
% Drop this file in +geom/SurfaceTools.m
%
% First-pass feature set:
%   - edge extraction from a NURBSSurface side
%   - continuity reporting between surface edges
%   - two-edge fill / blend surface with G0 or G1 boundary options
%
% Edge names:
%   'u0' : S(domainU(1), v)
%   'u1' : S(domainU(2), v)
%   'v0' : S(u, domainV(1))
%   'v1' : S(u, domainV(2))

    methods (Static)

        function E = edge(S, side)
            side = lower(char(string(side)));

            if ~isa(S, 'geom.NURBSSurface')
                error('SurfaceTools:edge', 'S must be a geom.NURBSSurface.');
            end

            switch side
                case 'u0'
                    C = S.isoCurveV(S.domainU(1));
                case 'u1'
                    C = S.isoCurveV(S.domainU(2));
                case 'v0'
                    C = S.isoCurveU(S.domainV(1));
                case 'v1'
                    C = S.isoCurveU(S.domainV(2));
                otherwise
                    error('SurfaceTools:edge', ...
                        'side must be one of: u0, u1, v0, v1.');
            end

            E = struct();
            E.surface = S;
            E.side = side;
            E.curve = C;
        end


        function report = edgeReport(Sa, edgeA, Sb, edgeB, nSample)
            if nargin < 5 || isempty(nSample)
                nSample = 41;
            end

            ea = geom.SurfaceTools.edge(Sa, edgeA);
            eb = geom.SurfaceTools.edge(Sb, edgeB);

            [Ca, Cb, flipB] = geom.SurfaceTools.compatibleEdgeCurves(ea.curve, eb.curve);

            t = linspace(0, 1, nSample).';
            ua = geom.SurfaceTools.map01ToDomain(t, Ca.domain);
            ub = geom.SurfaceTools.map01ToDomain(t, Cb.domain);

            Pa = Ca.evaluate(ua);
            Pb = Cb.evaluate(ub);

            gaps = vecnorm(Pa - Pb, 2, 2);

            ang = zeros(nSample,1);
            for k = 1:nSample
                uvA = geom.SurfaceTools.edgeParam01(Sa, edgeA, t(k));

                if flipB
                    tb = 1 - t(k);
                else
                    tb = t(k);
                end
                uvB = geom.SurfaceTools.edgeParam01(Sb, edgeB, tb);

                Na = Sa.normal(uvA(1), uvA(2));
                Nb = Sb.normal(uvB(1), uvB(2));

                if norm(Na) < eps || norm(Nb) < eps
                    ang(k) = NaN;
                else
                    c = abs(dot(Na, Nb) / (norm(Na)*norm(Nb)));
                    c = min(1,max(-1,c));
                    ang(k) = acosd(c);
                end
            end

            report = struct();
            report.maxGap = max(gaps);
            report.rmsGap = sqrt(mean(gaps.^2));
            report.maxNormalAngleDeg = max(ang, [], 'omitnan');
            report.rmsNormalAngleDeg = sqrt(mean(ang.^2, 'omitnan'));
            report.nSample = nSample;
        end

        function F = twoEdgeFill(E1, E2, varargin)
        %TWOEDGEFILL Create a CATIA-style two-boundary fill/blend surface.
        %
        % F = geom.SurfaceTools.twoEdgeFill(E1,E2, ...
        %     'Continuity1','G1', ...
        %     'Continuity2','G1', ...
        %     'DegreeV',3, ...
        %     'Scale1',1.0, ...
        %     'Scale2',1.0)
        %
        % The fill direction is v:
        %   F(u,0) = E1
        %   F(u,1) = E2
        %
        % For G1, tangent signs are chosen automatically so each parent
        % cross-boundary derivative points into the gap.
        
            pa = inputParser;
            addParameter(pa, 'Continuity1', 'G1');
            addParameter(pa, 'Continuity2', 'G1');
            addParameter(pa, 'DegreeV', 3);
            addParameter(pa, 'Scale1', 1.0);
            addParameter(pa, 'Scale2', 1.0);
            addParameter(pa, 'MinTangent', 1e-12);
            parse(pa, varargin{:});
            opts = pa.Results;
        
            cont1 = upper(string(opts.Continuity1));
            cont2 = upper(string(opts.Continuity2));
            q = opts.DegreeV;
        
            if q < 1 || q ~= floor(q)
                error('SurfaceTools:twoEdgeFill', ...
                    'DegreeV must be a positive integer.');
            end
        
            if q < 3 && (cont1 == "G1" || cont2 == "G1")
                error('SurfaceTools:twoEdgeFill', ...
                    'DegreeV must be at least 3 for G1 fill.');
            end
        
            [C1, C2, flip2] = geom.SurfaceTools.compatibleEdgeCurves(E1.curve, E2.curve);
        
            p = C1.p;
            U = C1.U;
            nCtrl = size(C1.P,1);
        
            nvCtrl = q + 1;
            V = [zeros(1,q+1), ones(1,q+1)];
        
            P = zeros(nCtrl, nvCtrl, 3);
            W = ones(nCtrl, nvCtrl);
        
            P1 = C1.P;
            P2 = C2.P;
            W1 = C1.W(:);
            W2 = C2.W(:);
        
            P(:,1,:) = reshape(P1, nCtrl, 1, 3);
            P(:,end,:) = reshape(P2, nCtrl, 1, 3);
            W(:,1) = W1;
            W(:,end) = W2;
        
            % Default interior is linear interpolation in Cartesian space.
            for j = 2:nvCtrl-1
                a = (j-1)/(nvCtrl-1);
                P(:,j,:) = reshape((1-a)*P1 + a*P2, nCtrl, 1, 3);
                W(:,j) = (1-a)*W1 + a*W2;
            end
        
            if q >= 3
                t1 = geom.SurfaceTools.grevilleParamsForCurve(C1);
        
                t2 = t1;
                if flip2
                    t2 = 1 - t1;
                end
        
                % Raw parent cross-boundary derivatives.
                D1 = geom.SurfaceTools.crossBoundaryDerivativeAlongEdge( ...
                    E1.surface, E1.side, t1, +1);
        
                D2 = geom.SurfaceTools.crossBoundaryDerivativeAlongEdge( ...
                    E2.surface, E2.side, t2, +1);
        
                if flip2
                    D2 = flipud(D2);
                end
        
                % Gap direction from boundary 1 to boundary 2.
                G = P2 - P1;
        
                % Robust tangent orientation:
                %   edge 1 tangent should point toward edge 2:  dot(D1,  G) > 0
                %   edge 2 tangent should point toward edge 1:  dot(D2, -G) > 0
                for i = 1:nCtrl
                    gi = G(i,:);
        
                    if norm(gi) > opts.MinTangent
                        if dot(D1(i,:), gi) < 0
                            D1(i,:) = -D1(i,:);
                        end
        
                        if dot(D2(i,:), -gi) < 0
                            D2(i,:) = -D2(i,:);
                        end
                    end
                end
        
                D1 = opts.Scale1 * D1;
                D2 = opts.Scale2 * D2;
        
                if cont1 == "G1"
                    for i = 1:nCtrl
                        P(i,2,:) = reshape(P1(i,:) + D1(i,:)/q, 1, 1, 3);
                    end
                    W(:,2) = W(:,1);
                end
        
                if cont2 == "G1"
                    for i = 1:nCtrl
                        P(i,end-1,:) = reshape(P2(i,:) + D2(i,:)/q, 1, 1, 3);
                    end
                    W(:,end-1) = W(:,end);
                end
            end
        
            F = geom.NURBSSurface(P, p, q, U, V, W);
        end

        function F = twoEdgeFillFair(E1, E2, varargin)
        %TWOEDGEFILLFAIR Create a two-edge fill with extra interior fairness rows.
        %
        % This keeps the same CATIA-like G0/G1 boundary logic as twoEdgeFill(),
        % but uses a higher degree in v and relaxes interior rows using cubic Hermite
        % interpolation between the two boundary curves/tangent fields.
        %
        % Example:
        %   F = geom.SurfaceTools.twoEdgeFillFair(e1,e2, ...
        %       'Continuity1','G1', ...
        %       'Continuity2','G1', ...
        %       'DegreeV',5, ...
        %       'Scale1',0.65, ...
        %       'Scale2',0.65);
        
            pa = inputParser;
            addParameter(pa, 'Continuity1', 'G1');
            addParameter(pa, 'Continuity2', 'G1');
            addParameter(pa, 'DegreeV', 5);
            addParameter(pa, 'Scale1', 1.0);
            addParameter(pa, 'Scale2', 1.0);
            addParameter(pa, 'Fairness', 0.5);
            addParameter(pa, 'MinTangent', 1e-12);
            parse(pa, varargin{:});
            opts = pa.Results;
        
            cont1 = upper(string(opts.Continuity1));
            cont2 = upper(string(opts.Continuity2));
            q = opts.DegreeV;
        
            if q < 3
                error('SurfaceTools:twoEdgeFillFair', ...
                    'DegreeV must be at least 3.');
            end
        
            [C1, C2, flip2] = geom.SurfaceTools.compatibleEdgeCurves(E1.curve, E2.curve);
        
            p = C1.p;
            U = C1.U;
            nCtrl = size(C1.P,1);
        
            nvCtrl = q + 1;
            V = [zeros(1,q+1), ones(1,q+1)];
        
            P = zeros(nCtrl, nvCtrl, 3);
            W = ones(nCtrl, nvCtrl);
        
            P1 = C1.P;
            P2 = C2.P;
            W1 = C1.W(:);
            W2 = C2.W(:);
        
            t1 = geom.SurfaceTools.grevilleParamsForCurve(C1);
            t2 = t1;
            if flip2
                t2 = 1 - t1;
            end
        
            D1 = geom.SurfaceTools.crossBoundaryDerivativeAlongEdge( ...
                E1.surface, E1.side, t1, +1);
        
            D2 = geom.SurfaceTools.crossBoundaryDerivativeAlongEdge( ...
                E2.surface, E2.side, t2, +1);
        
            if flip2
                D2 = flipud(D2);
            end
        
            G = P2 - P1;
        
            for i = 1:nCtrl
                gi = G(i,:);
        
                if norm(gi) > opts.MinTangent
                    if dot(D1(i,:), gi) < 0
                        D1(i,:) = -D1(i,:);
                    end
        
                    if dot(D2(i,:), -gi) < 0
                        D2(i,:) = -D2(i,:);
                    end
                end
            end
        
            if cont1 ~= "G1"
                D1 = G;
            else
                D1 = opts.Scale1 * D1;
            end
        
            if cont2 ~= "G1"
                D2 = -G;
            else
                D2 = opts.Scale2 * D2;
            end
        
            % Clamp fairness blend.
            fair = max(0, min(1, opts.Fairness));
        
            for j = 1:nvCtrl
                s = (j-1) / (nvCtrl-1);
        
                h00 =  2*s^3 - 3*s^2 + 1;
                h10 =      s^3 - 2*s^2 + s;
                h01 = -2*s^3 + 3*s^2;
                h11 =      s^3 -   s^2;
        
                PHermite = h00*P1 + h10*D1 + h01*P2 + h11*(-D2);
                PLinear  = (1-s)*P1 + s*P2;
        
                Pj = (1-fair)*PHermite + fair*PLinear;
        
                P(:,j,:) = reshape(Pj, nCtrl, 1, 3);
                W(:,j) = (1-s)*W1 + s*W2;
            end
        
            % Restore exact boundary rows.
            P(:,1,:) = reshape(P1, nCtrl, 1, 3);
            P(:,end,:) = reshape(P2, nCtrl, 1, 3);
            W(:,1) = W1;
            W(:,end) = W2;
        
            % Restore exact G1 rows if requested.
            if cont1 == "G1"
                for i = 1:nCtrl
                    P(i,2,:) = reshape(P1(i,:) + D1(i,:)/q, 1, 1, 3);
                end
                W(:,2) = W(:,1);
            end
        
            if cont2 == "G1"
                for i = 1:nCtrl
                    P(i,end-1,:) = reshape(P2(i,:) + D2(i,:)/q, 1, 1, 3);
                end
                W(:,end-1) = W(:,end);
            end
        
            F = geom.NURBSSurface(P, p, q, U, V, W);
        end



        function [S2, report] = matchG1Edge(S, edgeName, Sref, refEdgeName, varargin)
            pa = inputParser;
            addParameter(pa, 'Scale', 1.0);
            addParameter(pa, 'Samples', 41);
            parse(pa, varargin{:});
            opts = pa.Results;

            side = lower(char(string(edgeName)));
            refSide = lower(char(string(refEdgeName)));

            E = geom.SurfaceTools.edge(S, side);
            Eref = geom.SurfaceTools.edge(Sref, refSide);

            [C, Cref, flipRef] = geom.SurfaceTools.compatibleEdgeCurves(E.curve, Eref.curve);

            P = S.P;
            W = S.W;

            nCtrl = size(C.P,1);
            t = geom.SurfaceTools.grevilleParamsForCurve(C);

            tref = t;
            if flipRef
                tref = 1 - t;
            end

            Dref = geom.SurfaceTools.crossBoundaryDerivativeAlongEdge( ...
                Sref, refSide, tref, -1);

            if flipRef
                Dref = flipud(Dref);
            end

            Dref = opts.Scale * Dref;

            switch side
                case 'v0'
                    q = S.q;
                    P(:,1,:) = reshape(Cref.P, nCtrl, 1, 3);
                    for i = 1:nCtrl
                        P(i,2,:) = reshape(Cref.P(i,:) + Dref(i,:)/q, 1,1,3);
                    end

                case 'v1'
                    q = S.q;
                    P(:,end,:) = reshape(Cref.P, nCtrl, 1, 3);
                    for i = 1:nCtrl
                        P(i,end-1,:) = reshape(Cref.P(i,:) - Dref(i,:)/q, 1,1,3);
                    end

                case 'u0'
                    p = S.p;
                    P(1,:,:) = reshape(Cref.P, 1, nCtrl, 3);
                    for j = 1:nCtrl
                        P(2,j,:) = reshape(Cref.P(j,:) + Dref(j,:)/p, 1,1,3);
                    end

                case 'u1'
                    p = S.p;
                    P(end,:,:) = reshape(Cref.P, 1, nCtrl, 3);
                    for j = 1:nCtrl
                        P(end-1,j,:) = reshape(Cref.P(j,:) - Dref(j,:)/p, 1,1,3);
                    end

                otherwise
                    error('SurfaceTools:matchG1Edge', ...
                        'edgeName must be u0,u1,v0,v1.');
            end

            S2 = geom.NURBSSurface(P, S.p, S.q, S.U, S.V, W);
            report = geom.SurfaceTools.edgeReport(S2, edgeName, Sref, refEdgeName, opts.Samples);
        end


        function F = twoEdgeFillLaw(E1, E2, varargin)
        %TWOEDGEFILLLAW Two-edge fill with variable tangent scale laws.
        %
        % Example:
        %   F = geom.SurfaceTools.twoEdgeFillLaw(e1,e2, ...
        %       'Continuity1','G1', ...
        %       'Continuity2','G1', ...
        %       'DegreeV',5, ...
        %       'ScaleLaw1',@(t) 0.4 + 0.4*sin(pi*t), ...
        %       'ScaleLaw2',@(t) 0.7*ones(size(t)), ...
        %       'Fairness',0.25);
        
            pa = inputParser;
            addParameter(pa, 'Continuity1', 'G1');
            addParameter(pa, 'Continuity2', 'G1');
            addParameter(pa, 'DegreeV', 5);
            addParameter(pa, 'ScaleLaw1', @(t) ones(size(t)));
            addParameter(pa, 'ScaleLaw2', @(t) ones(size(t)));
            addParameter(pa, 'Fairness', 0.35);
            addParameter(pa, 'MinTangent', 1e-12);
            parse(pa, varargin{:});
            opts = pa.Results;
        
            cont1 = upper(string(opts.Continuity1));
            cont2 = upper(string(opts.Continuity2));
            q = opts.DegreeV;
        
            if q < 3
                error('SurfaceTools:twoEdgeFillLaw', ...
                    'DegreeV must be at least 3.');
            end
        
            [C1, C2, flip2] = geom.SurfaceTools.compatibleEdgeCurves(E1.curve, E2.curve);
        
            p = C1.p;
            U = C1.U;
            nCtrl = size(C1.P,1);
        
            nvCtrl = q + 1;
            V = [zeros(1,q+1), ones(1,q+1)];
        
            P = zeros(nCtrl, nvCtrl, 3);
            W = ones(nCtrl, nvCtrl);
        
            P1 = C1.P;
            P2 = C2.P;
            W1 = C1.W(:);
            W2 = C2.W(:);
        
            t1 = geom.SurfaceTools.grevilleParamsForCurve(C1);
            t2 = t1;
            if flip2
                t2 = 1 - t1;
            end
        
            D1 = geom.SurfaceTools.crossBoundaryDerivativeAlongEdge( ...
                E1.surface, E1.side, t1, +1);
        
            D2 = geom.SurfaceTools.crossBoundaryDerivativeAlongEdge( ...
                E2.surface, E2.side, t2, +1);
        
            if flip2
                D2 = flipud(D2);
            end
        
            G = P2 - P1;
        
            for i = 1:nCtrl
                gi = G(i,:);
                if norm(gi) > opts.MinTangent
                    if dot(D1(i,:), gi) < 0
                        D1(i,:) = -D1(i,:);
                    end
                    if dot(D2(i,:), -gi) < 0
                        D2(i,:) = -D2(i,:);
                    end
                end
            end
        
            s1 = opts.ScaleLaw1(t1);
            s2 = opts.ScaleLaw2(t1);
        
            if isscalar(s1), s1 = repmat(s1, nCtrl, 1); end
            if isscalar(s2), s2 = repmat(s2, nCtrl, 1); end
        
            s1 = s1(:);
            s2 = s2(:);
        
            if numel(s1) ~= nCtrl || numel(s2) ~= nCtrl
                error('SurfaceTools:twoEdgeFillLaw', ...
                    'ScaleLaw1 and ScaleLaw2 must return scalar or nCtrl values.');
            end
        
            if cont1 ~= "G1"
                D1 = G;
            else
                D1 = D1 .* s1;
            end
        
            if cont2 ~= "G1"
                D2 = -G;
            else
                D2 = D2 .* s2;
            end
        
            fair = max(0, min(1, opts.Fairness));
        
            for j = 1:nvCtrl
                s = (j-1) / (nvCtrl-1);
        
                h00 =  2*s^3 - 3*s^2 + 1;
                h10 =      s^3 - 2*s^2 + s;
                h01 = -2*s^3 + 3*s^2;
                h11 =      s^3 -   s^2;
        
                PHermite = h00*P1 + h10*D1 + h01*P2 + h11*(-D2);
                PLinear  = (1-s)*P1 + s*P2;
        
                Pj = (1-fair)*PHermite + fair*PLinear;
        
                P(:,j,:) = reshape(Pj, nCtrl, 1, 3);
                W(:,j) = (1-s)*W1 + s*W2;
            end
        
            P(:,1,:) = reshape(P1, nCtrl, 1, 3);
            P(:,end,:) = reshape(P2, nCtrl, 1, 3);
            W(:,1) = W1;
            W(:,end) = W2;
        
            if cont1 == "G1"
                for i = 1:nCtrl
                    P(i,2,:) = reshape(P1(i,:) + D1(i,:)/q, 1, 1, 3);
                end
                W(:,2) = W(:,1);
            end
        
            if cont2 == "G1"
                for i = 1:nCtrl
                    P(i,end-1,:) = reshape(P2(i,:) + D2(i,:)/q, 1, 1, 3);
                end
                W(:,end-1) = W(:,end);
            end
        
            F = geom.NURBSSurface(P, p, q, U, V, W);
        end


        function S2 = relaxInterior(S, varargin)
        %RELAXINTERIOR Laplacian relax only unconstrained interior control points.
        %
        % This preserves boundary rows/columns and, optionally, first tangent rows.
        %
        % Example:
        %   F2 = geom.SurfaceTools.relaxInterior(F, ...
        %       'Iterations',20, ...
        %       'Lambda',0.35, ...
        %       'PreserveTangencyRows',true);
        
            pa = inputParser;
            addParameter(pa, 'Iterations', 10);
            addParameter(pa, 'Lambda', 0.25);
            addParameter(pa, 'PreserveTangencyRows', true);
            parse(pa, varargin{:});
            opts = pa.Results;
        
            P = S.P;
            W = S.W;
        
            ni = size(P,1);
            nj = size(P,2);
        
            lam = max(0, min(1, opts.Lambda));
            nIter = max(0, floor(opts.Iterations));
        
            if opts.PreserveTangencyRows
                iRange = 3:ni-2;
                jRange = 3:nj-2;
            else
                iRange = 2:ni-1;
                jRange = 2:nj-1;
            end
        
            if isempty(iRange) || isempty(jRange)
                S2 = S;
                return;
            end
        
            for it = 1:nIter
                Pold = P;
                for i = iRange
                    for j = jRange
                        nbr = 0.25 * ( ...
                            reshape(Pold(i-1,j,:),1,3) + ...
                            reshape(Pold(i+1,j,:),1,3) + ...
                            reshape(Pold(i,j-1,:),1,3) + ...
                            reshape(Pold(i,j+1,:),1,3));
        
                        cur = reshape(Pold(i,j,:),1,3);
                        P(i,j,:) = reshape((1-lam)*cur + lam*nbr, 1, 1, 3);
                    end
                end
            end
        
            S2 = geom.NURBSSurface(P, S.p, S.q, S.U, S.V, W);
        end

        function [S2, report] = matchG1EdgeRobust(S, edgeName, Sref, refEdgeName, varargin)
        %MATCHG1EDGEROBUST Force one edge of S to be G0/G1 to a reference edge.
        %
        % This is intended for CATIA-style "Match Surface" behavior.
        %
        % It:
        %   1) matches the boundary control row/column to the reference edge
        %   2) chooses the reference tangent sign so it points into the existing S patch
        %   3) sets the first interior row/column to enforce tangent-plane continuity
        %
        % Example:
        %   [S2, rpt] = geom.SurfaceTools.matchG1EdgeRobust( ...
        %       Sfill, 'v0', Sparent, 'v1', ...
        %       'MagnitudeMode','target', ...
        %       'Scale',1.0);
        
            pa = inputParser;
            addParameter(pa, 'Scale', 1.0);
            addParameter(pa, 'Samples', 51);
            addParameter(pa, 'MagnitudeMode', 'target'); % 'target', 'reference', 'blend'
            addParameter(pa, 'Blend', 0.5);
            addParameter(pa, 'MinTangent', 1e-12);
            parse(pa, varargin{:});
            opts = pa.Results;
        
            side    = lower(char(string(edgeName)));
            refSide = lower(char(string(refEdgeName)));
        
            E    = geom.SurfaceTools.edge(S, side);
            Eref = geom.SurfaceTools.edge(Sref, refSide);
        
            [C, Cref, flipRef] = geom.SurfaceTools.compatibleEdgeCurves(E.curve, Eref.curve);
        
            % This direct control-row method requires the matched edge curve to have
            % the same number of control points as the actual target edge.
            nCtrl = size(C.P,1);
        
            switch side
                case {'v0','v1'}
                    if nCtrl ~= size(S.P,1)
                        error('SurfaceTools:matchG1EdgeRobust', ...
                            ['Target surface edge is not control-compatible after curve compatibility. ', ...
                             'Refine/elevate the target surface first.']);
                    end
                    degCross = S.q;
        
                case {'u0','u1'}
                    if nCtrl ~= size(S.P,2)
                        error('SurfaceTools:matchG1EdgeRobust', ...
                            ['Target surface edge is not control-compatible after curve compatibility. ', ...
                             'Refine/elevate the target surface first.']);
                    end
                    degCross = S.p;
        
                otherwise
                    error('SurfaceTools:matchG1EdgeRobust', ...
                        'edgeName must be u0,u1,v0,v1.');
            end
        
            if degCross < 1
                error('SurfaceTools:matchG1EdgeRobust', ...
                    'Target surface degree normal to edge must be at least 1.');
            end
        
            P = S.P;
            W = S.W;
        
            % Use Greville parameters of the target edge control curve.
            t = geom.SurfaceTools.grevilleParamsForCurve(C);
        
            tref = t;
            if flipRef
                tref = 1 - t;
            end
        
            % Current target inward derivative. This tells us which side of the seam
            % the target patch currently occupies.
            Dcur = geom.SurfaceTools.crossBoundaryDerivativeAlongEdge( ...
                S, side, t, +1);
        
            % Reference inward derivative.
            Dref = geom.SurfaceTools.crossBoundaryDerivativeAlongEdge( ...
                Sref, refSide, tref, +1);
        
            if flipRef
                Dref = flipud(Dref);
            end
        
            % Choose reference tangent sign so it points into the target patch.
            for i = 1:nCtrl
                if norm(Dcur(i,:)) > opts.MinTangent && norm(Dref(i,:)) > opts.MinTangent
                    if dot(Dref(i,:), Dcur(i,:)) < 0
                        Dref(i,:) = -Dref(i,:);
                    end
                end
            end
        
            % Magnitude policy. G1 only needs direction/tangent plane, not equal
            % derivative magnitude.
            mode = lower(char(string(opts.MagnitudeMode)));
        
            for i = 1:nCtrl
                nr = norm(Dref(i,:));
                nt = norm(Dcur(i,:));
        
                if nr < opts.MinTangent
                    continue;
                end
        
                dHat = Dref(i,:) / nr;
        
                switch mode
                    case 'reference'
                        mag = nr;
        
                    case 'target'
                        if nt > opts.MinTangent
                            mag = nt;
                        else
                            mag = nr;
                        end
        
                    case 'blend'
                        b = max(0,min(1,opts.Blend));
                        if nt > opts.MinTangent
                            mag = (1-b)*nr + b*nt;
                        else
                            mag = nr;
                        end
        
                    otherwise
                        error('SurfaceTools:matchG1EdgeRobust', ...
                            'MagnitudeMode must be reference, target, or blend.');
                end
        
                Dref(i,:) = opts.Scale * mag * dHat;
            end
        
            Pedge = Cref.P;
        
            switch side
                case 'v0'
                    P(:,1,:) = reshape(Pedge, nCtrl, 1, 3);
                    for i = 1:nCtrl
                        P(i,2,:) = reshape(Pedge(i,:) + Dref(i,:)/degCross, 1, 1, 3);
                    end
        
                case 'v1'
                    P(:,end,:) = reshape(Pedge, nCtrl, 1, 3);
                    for i = 1:nCtrl
                        P(i,end-1,:) = reshape(Pedge(i,:) + Dref(i,:)/degCross, 1, 1, 3);
                    end
        
                case 'u0'
                    P(1,:,:) = reshape(Pedge, 1, nCtrl, 3);
                    for j = 1:nCtrl
                        P(2,j,:) = reshape(Pedge(j,:) + Dref(j,:)/degCross, 1, 1, 3);
                    end
        
                case 'u1'
                    P(end,:,:) = reshape(Pedge, 1, nCtrl, 3);
                    for j = 1:nCtrl
                        P(end-1,j,:) = reshape(Pedge(j,:) + Dref(j,:)/degCross, 1, 1, 3);
                    end
            end
        
            S2 = geom.NURBSSurface(P, S.p, S.q, S.U, S.V, W);
        
            report = geom.SurfaceTools.edgeReport( ...
                S2, edgeName, Sref, refEdgeName, opts.Samples);
        end

        
        function report = edgeG1Report(Sa, edgeA, Sb, edgeB, nSample)
        %EDGEG1REPORT More detailed seam continuity report.
        %
        % Reports:
        %   G0 gap
        %   normal-angle mismatch
        %   cross-boundary tangent parallelism mismatch
        %
        % For true G1, normal angle should be near zero. The cross-boundary
        % derivative may point opposite across a seam, so this reports parallelism:
        % min(angle, 180-angle).
        
            if nargin < 5 || isempty(nSample)
                nSample = 51;
            end
        
            base = geom.SurfaceTools.edgeReport(Sa, edgeA, Sb, edgeB, nSample);
        
            ea = geom.SurfaceTools.edge(Sa, edgeA);
            eb = geom.SurfaceTools.edge(Sb, edgeB);
            [Ca, ~, flipB] = geom.SurfaceTools.compatibleEdgeCurves(ea.curve, eb.curve);
        
            t = linspace(0,1,nSample).';
        
            crossAng = zeros(nSample,1);
        
            for k = 1:nSample
                ta = t(k);
                tb = t(k);
                if flipB
                    tb = 1 - tb;
                end
        
                Da = geom.SurfaceTools.crossBoundaryDerivativeAlongEdge( ...
                    Sa, edgeA, ta, +1);
        
                Db = geom.SurfaceTools.crossBoundaryDerivativeAlongEdge( ...
                    Sb, edgeB, tb, +1);
        
                if norm(Da) < eps || norm(Db) < eps
                    crossAng(k) = NaN;
                else
                    c = dot(Da, Db) / (norm(Da)*norm(Db));
                    c = min(1,max(-1,c));
                    a = acosd(c);
                    crossAng(k) = min(a, 180-a);
                end
            end
        
            report = base;
            report.maxCrossTangentParallelAngleDeg = max(crossAng, [], 'omitnan');
            report.rmsCrossTangentParallelAngleDeg = sqrt(mean(crossAng.^2, 'omitnan'));
            report.nSample = nSample;
            report.edgeA = edgeA;
            report.edgeB = edgeB;
            report.curveDegree = Ca.p;
        end

        function F = fourEdgeFill(Ev0, Ev1, Eu0, Eu1, varargin)
        %FOUREDGEFILL Four-sided CATIA-style fill patch with G0/G1 edge options.
        %
        % Boundary convention:
        %   Ev0 -> F(u,0)  bottom edge
        %   Ev1 -> F(u,1)  top edge
        %   Eu0 -> F(0,v)  left edge
        %   Eu1 -> F(1,v)  right edge
        %
        % This version:
        %   - degree-elevates the fill space by default
        %   - restores exact G0 boundaries
        %   - converts sampled derivative fields into derivative-control rows
        %     before setting G1 rows/columns.
        
            pa = inputParser;
            addParameter(pa, 'ContinuityV0', 'G1');
            addParameter(pa, 'ContinuityV1', 'G1');
            addParameter(pa, 'ContinuityU0', 'G1');
            addParameter(pa, 'ContinuityU1', 'G1');
        
            addParameter(pa, 'DegreeU', 5);
            addParameter(pa, 'DegreeV', 5);
        
            addParameter(pa, 'ScaleV0', 1.0);
            addParameter(pa, 'ScaleV1', 1.0);
            addParameter(pa, 'ScaleU0', 1.0);
            addParameter(pa, 'ScaleU1', 1.0);
        
            addParameter(pa, 'CornerTolerance', 1e-7);
            addParameter(pa, 'MinTangent', 1e-12);
            parse(pa, varargin{:});
            opts = pa.Results;
        
            contV0 = upper(string(opts.ContinuityV0));
            contV1 = upper(string(opts.ContinuityV1));
            contU0 = upper(string(opts.ContinuityU0));
            contU1 = upper(string(opts.ContinuityU1));
        
            % ------------------------------------------------------------------
            % 1. Auto-orient boundary curves.
            % ------------------------------------------------------------------
            Cb0 = Ev0.curve;
            Ct0 = Ev1.curve;
            Cl0 = Eu0.curve;
            Cr0 = Eu1.curve;
        
            bestErr = inf;
            best = [];
        
            for fb = 0:1
                CbTry = maybeReverse(Cb0, fb);
                for ft = 0:1
                    CtTry = maybeReverse(Ct0, ft);
                    for fl = 0:1
                        ClTry = maybeReverse(Cl0, fl);
                        for fr = 0:1
                            CrTry = maybeReverse(Cr0, fr);
        
                            P00a = CbTry.evaluate(CbTry.domain(1));
                            P10a = CbTry.evaluate(CbTry.domain(2));
                            P01a = CtTry.evaluate(CtTry.domain(1));
                            P11a = CtTry.evaluate(CtTry.domain(2));
        
                            P00b = ClTry.evaluate(ClTry.domain(1));
                            P01b = ClTry.evaluate(ClTry.domain(2));
                            P10b = CrTry.evaluate(CrTry.domain(1));
                            P11b = CrTry.evaluate(CrTry.domain(2));
        
                            err = norm(P00a-P00b) + norm(P10a-P10b) + ...
                                  norm(P01a-P01b) + norm(P11a-P11b);
        
                            if err < bestErr
                                bestErr = err;
                                best.Cb = CbTry;
                                best.Ct = CtTry;
                                best.Cl = ClTry;
                                best.Cr = CrTry;
                                best.flipB = logical(fb);
                                best.flipT = logical(ft);
                                best.flipL = logical(fl);
                                best.flipR = logical(fr);
                            end
                        end
                    end
                end
            end
        
            if bestErr > opts.CornerTolerance
                warning('SurfaceTools:fourEdgeFill:CornerMismatch', ...
                    'Best four-edge corner closure error = %.3e.', bestErr);
            end
        
            Cb = best.Cb;
            Ct = best.Ct;
            Cl = best.Cl;
            Cr = best.Cr;
        
            % ------------------------------------------------------------------
            % 2. Degree elevate first, then make opposite pairs compatible.
            % ------------------------------------------------------------------
            pFill = max([opts.DegreeU, Cb.p, Ct.p, 3]);
            qFill = max([opts.DegreeV, Cl.p, Cr.p, 3]);
        
            if Cb.p < pFill, Cb = Cb.elevate(pFill - Cb.p); end
            if Ct.p < pFill, Ct = Ct.elevate(pFill - Ct.p); end
            if Cl.p < qFill, Cl = Cl.elevate(qFill - Cl.p); end
            if Cr.p < qFill, Cr = Cr.elevate(qFill - Cr.p); end
        
            [uvCurves, p, U] = geom.NURBSSurface.makeCompatibleCurves({Cb, Ct});
            Cb = uvCurves{1};
            Ct = uvCurves{2};
        
            [lrCurves, q, V] = geom.NURBSSurface.makeCompatibleCurves({Cl, Cr});
            Cl = lrCurves{1};
            Cr = lrCurves{2};
        
            nu = size(Cb.P,1);
            nv = size(Cl.P,1);
        
            % ------------------------------------------------------------------
            % 3. Initial Coons-like control net.
            % ------------------------------------------------------------------
            P = zeros(nu, nv, 3);
            W = ones(nu, nv);
        
            P00 = 0.5 * (Cb.evaluate(Cb.domain(1)) + Cl.evaluate(Cl.domain(1)));
            P10 = 0.5 * (Cb.evaluate(Cb.domain(2)) + Cr.evaluate(Cr.domain(1)));
            P01 = 0.5 * (Ct.evaluate(Ct.domain(1)) + Cl.evaluate(Cl.domain(2)));
            P11 = 0.5 * (Ct.evaluate(Ct.domain(2)) + Cr.evaluate(Cr.domain(2)));
        
            for i = 1:nu
                uhat = greville01(i, p, U);
                Pb = Cb.P(i,:);
                Pt = Ct.P(i,:);
        
                for j = 1:nv
                    vhat = greville01(j, q, V);
                    Pl = Cl.P(j,:);
                    Pr = Cr.P(j,:);
        
                    Su = (1-vhat)*Pb + vhat*Pt;
                    Sv = (1-uhat)*Pl + uhat*Pr;
        
                    Sb = (1-uhat)*(1-vhat)*P00 + ...
                          uhat   *(1-vhat)*P10 + ...
                         (1-uhat)*vhat   *P01 + ...
                          uhat   *vhat   *P11;
        
                    P(i,j,:) = reshape(Su + Sv - Sb, 1, 1, 3);
                end
            end
        
            restoreBoundaries();
        
            % ------------------------------------------------------------------
            % 4. G1 rows/columns.
            % ------------------------------------------------------------------
            if contV0 == "G1"
                t = grevilleVector01(Cb);
                D = edgeDerivativeSamples(Ev0, t, best.flipB);
                G = squeeze(P(:,end,:) - P(:,1,:));
                D = orientIntoGap(D, G, opts.MinTangent);
                D = opts.ScaleV0 * D;
        
                % Convert pointwise derivative samples into spline control values
                % along the u direction.
                D = derivativeSamplesToControl(D, p, U);
        
                for i = 1:nu
                    P(i,2,:) = reshape(squeeze(P(i,1,:)).' + D(i,:)/q, 1, 1, 3);
                end
                W(:,2) = W(:,1);
            end
        
            if contV1 == "G1"
                t = grevilleVector01(Ct);
                D = edgeDerivativeSamples(Ev1, t, best.flipT);
                G = squeeze(P(:,1,:) - P(:,end,:));
                D = orientIntoGap(D, G, opts.MinTangent);
                D = opts.ScaleV1 * D;
        
                % Convert pointwise derivative samples into spline control values
                % along the u direction.
                D = derivativeSamplesToControl(D, p, U);
        
                for i = 1:nu
                    P(i,end-1,:) = reshape(squeeze(P(i,end,:)).' + D(i,:)/q, 1, 1, 3);
                end
                W(:,end-1) = W(:,end);
            end
        
            if contU0 == "G1"
                t = grevilleVector01(Cl);
                D = edgeDerivativeSamples(Eu0, t, best.flipL);
                G = squeeze(P(end,:,:) - P(1,:,:));
                D = orientIntoGap(D, G, opts.MinTangent);
                D = opts.ScaleU0 * D;
        
                % Convert pointwise derivative samples into spline control values
                % along the v direction.
                D = derivativeSamplesToControl(D, q, V);
        
                for j = 1:nv
                    P(2,j,:) = reshape(squeeze(P(1,j,:)).' + D(j,:)/p, 1, 1, 3);
                end
                W(2,:) = W(1,:);
            end
        
            if contU1 == "G1"
                t = grevilleVector01(Cr);
                D = edgeDerivativeSamples(Eu1, t, best.flipR);
                G = squeeze(P(1,:,:) - P(end,:,:));
                D = orientIntoGap(D, G, opts.MinTangent);
                D = opts.ScaleU1 * D;
        
                % Convert pointwise derivative samples into spline control values
                % along the v direction.
                D = derivativeSamplesToControl(D, q, V);
        
                for j = 1:nv
                    P(end-1,j,:) = reshape(squeeze(P(end,j,:)).' + D(j,:)/p, 1, 1, 3);
                end
                W(end-1,:) = W(end,:);
            end
        
            % Corner-adjacent compromise. This is intentionally mild; the main
            % improvement comes from using derivative-control rows.
            blendCorner(2,    2);
            blendCorner(nu-1, 2);
            blendCorner(2,    nv-1);
            blendCorner(nu-1, nv-1);
        
            restoreBoundaries();
        
            F = geom.NURBSSurface(P, p, q, U, V, W);
        
            % ------------------------------------------------------------------
            % Local helpers
            % ------------------------------------------------------------------
            function C = maybeReverse(Cin, doFlip)
                C = Cin;
                if doFlip
                    C = C.reverse();
                end
            end
        
            function restoreBoundaries()
                P(:,1,:)   = reshape(Cb.P, nu, 1, 3);
                P(:,end,:) = reshape(Ct.P, nu, 1, 3);
                P(1,:,:)   = reshape(Cl.P, 1, nv, 3);
                P(end,:,:) = reshape(Cr.P, 1, nv, 3);
        
                W(:,1)     = Cb.W(:);
                W(:,end)   = Ct.W(:);
                W(1,:)     = Cl.W(:).';
                W(end,:)   = Cr.W(:).';
        
                W(1,1)       = 0.5*(Cb.W(1)   + Cl.W(1));
                W(end,1)     = 0.5*(Cb.W(end) + Cr.W(1));
                W(1,end)     = 0.5*(Ct.W(1)   + Cl.W(end));
                W(end,end)   = 0.5*(Ct.W(end) + Cr.W(end));
            end
        
            function blendCorner(ii,jj)
                if ii < 2 || ii > nu-1 || jj < 2 || jj > nv-1
                    return;
                end
        
                vals = zeros(0,3);
        
                if jj == 2 && contV0 == "G1"
                    vals(end+1,:) = squeeze(P(ii,2,:)).'; %#ok<AGROW>
                end
                if jj == nv-1 && contV1 == "G1"
                    vals(end+1,:) = squeeze(P(ii,nv-1,:)).'; %#ok<AGROW>
                end
                if ii == 2 && contU0 == "G1"
                    vals(end+1,:) = squeeze(P(2,jj,:)).'; %#ok<AGROW>
                end
                if ii == nu-1 && contU1 == "G1"
                    vals(end+1,:) = squeeze(P(nu-1,jj,:)).'; %#ok<AGROW>
                end
        
                if size(vals,1) > 1
                    P(ii,jj,:) = reshape(mean(vals,1), 1, 1, 3);
                end
            end
        
            function Dctrl = derivativeSamplesToControl(Dsamp, deg, K)
                n = size(Dsamp,1);
                A = zeros(n,n);
            
                dom = [K(deg+1), K(end-deg)];
            
                for rr = 1:n
                    tt = greville01(rr, deg, K);
                    uAbs = dom(1) + tt*(dom(2)-dom(1));
            
                    span = geom.BasisFunctions.FindSpan(n-1, deg, uAbs, K);
                    N = geom.BasisFunctions.BasisFuns(span, uAbs, deg, K);
            
                    % FindSpan is 1-based. BasisFuns corresponds to basis indices:
                    % span-deg, ..., span
                    cols = (span-deg):span;
            
                    A(rr, cols) = N;
                end
            
                Dctrl = A \ Dsamp;
            end
                    
            function g = greville01(ii, deg, K)
                if deg == 0
                    g = 0.0;
                else
                    gAbs = sum(K(ii+1:ii+deg)) / deg;
                    dom = [K(deg+1), K(end-deg)];
                    if abs(dom(2)-dom(1)) < eps
                        g = 0.0;
                    else
                        g = (gAbs - dom(1)) / (dom(2)-dom(1));
                    end
                end
                g = max(0,min(1,g));
            end
        
            function t = grevilleVector01(C)
                n = size(C.P,1);
                t = zeros(n,1);
                for kk = 1:n
                    t(kk) = greville01(kk, C.p, C.U);
                end
            end
        
            function D = edgeDerivativeSamples(E, t01, wasFlipped)
                if wasFlipped
                    tp = 1 - t01;
                else
                    tp = t01;
                end
        
                D = geom.SurfaceTools.crossBoundaryDerivativeAlongEdge( ...
                    E.surface, E.side, tp, +1);
        
                if wasFlipped
                    D = flipud(D);
                end
            end
        
            function D = orientIntoGap(D, G, minTan)
                for kk = 1:size(D,1)
                    gk = G(kk,:);
                    if norm(gk) > minTan && dot(D(kk,:), gk) < 0
                        D(kk,:) = -D(kk,:);
                    end
                end
            end
        end



        function report = edgeG1ReportInterior(Sa, edgeA, Sb, edgeB, varargin)
        %EDGEG1REPORTINTERIOR G1 report excluding corner/end effects.
        %
        % Example:
        %   r = geom.SurfaceTools.edgeG1ReportInterior(F4g1,'v0',Sb,'v1', ...
        %       'NSample',101, ...
        %       'TrimFraction',0.05);
        
            pa = inputParser;
            addParameter(pa, 'NSample', 101);
            addParameter(pa, 'TrimFraction', 0.05);
            parse(pa, varargin{:});
            opts = pa.Results;
        
            nSample = opts.NSample;
            trim = max(0, min(0.49, opts.TrimFraction));
        
            t = linspace(trim, 1-trim, nSample).';
        
            ea = geom.SurfaceTools.edge(Sa, edgeA);
            eb = geom.SurfaceTools.edge(Sb, edgeB);
            [Ca, Cb, flipB] = geom.SurfaceTools.compatibleEdgeCurves(ea.curve, eb.curve);
        
            ua = geom.SurfaceTools.map01ToDomain(t, Ca.domain);
        
            if flipB
                tb = 1 - t;
            else
                tb = t;
            end
            ub = geom.SurfaceTools.map01ToDomain(tb, Cb.domain);
        
            Pa = Ca.evaluate(ua);
            Pb = Cb.evaluate(ub);
        
            gaps = vecnorm(Pa - Pb, 2, 2);
        
            normalAng = zeros(nSample,1);
            crossAng  = zeros(nSample,1);
        
            for k = 1:nSample
                ta = t(k);
                tbk = tb(k);
        
                uvA = geom.SurfaceTools.edgeParam01(Sa, edgeA, ta);
                uvB = geom.SurfaceTools.edgeParam01(Sb, edgeB, tbk);
        
                Na = Sa.normal(uvA(1), uvA(2));
                Nb = Sb.normal(uvB(1), uvB(2));
        
                if norm(Na) < eps || norm(Nb) < eps
                    normalAng(k) = NaN;
                else
                    c = abs(dot(Na,Nb)/(norm(Na)*norm(Nb)));
                    c = min(1,max(-1,c));
                    normalAng(k) = acosd(c);
                end
        
                Da = geom.SurfaceTools.crossBoundaryDerivativeAlongEdge( ...
                    Sa, edgeA, ta, +1);
        
                Db = geom.SurfaceTools.crossBoundaryDerivativeAlongEdge( ...
                    Sb, edgeB, tbk, +1);
        
                if norm(Da) < eps || norm(Db) < eps
                    crossAng(k) = NaN;
                else
                    c = dot(Da,Db)/(norm(Da)*norm(Db));
                    c = min(1,max(-1,c));
                    a = acosd(c);
                    crossAng(k) = min(a, 180-a);
                end
            end
        
            report = struct();
            report.maxGap = max(gaps);
            report.rmsGap = sqrt(mean(gaps.^2));
            report.maxNormalAngleDeg = max(normalAng, [], 'omitnan');
            report.rmsNormalAngleDeg = sqrt(mean(normalAng.^2, 'omitnan'));
            report.maxCrossTangentParallelAngleDeg = max(crossAng, [], 'omitnan');
            report.rmsCrossTangentParallelAngleDeg = sqrt(mean(crossAng.^2, 'omitnan'));
            report.nSample = nSample;
            report.trimFraction = trim;
        end




    end

    methods (Static, Access = private)

        function [C1, C2, flip2] = compatibleEdgeCurves(C1in, C2in)
            C1 = C1in;
            C2 = C2in;

            a1 = C1.evaluate(C1.domain(1));
            b1 = C1.evaluate(C1.domain(2));
            a2 = C2.evaluate(C2.domain(1));
            b2 = C2.evaluate(C2.domain(2));

            dSame = norm(a1-a2) + norm(b1-b2);
            dFlip = norm(a1-b2) + norm(b1-a2);

            flip2 = dFlip < dSame;
            if flip2
                C2 = C2.reverse();
            end

            [curves, ~, ~] = geom.NURBSSurface.makeCompatibleCurves({C1,C2});
            C1 = curves{1};
            C2 = curves{2};
        end


        function uv = edgeParam01(S, side, t)
            side = lower(char(string(side)));
            du = S.domainU;
            dv = S.domainV;

            switch side
                case 'u0'
                    uv = [du(1), dv(1) + t*(dv(2)-dv(1))];
                case 'u1'
                    uv = [du(2), dv(1) + t*(dv(2)-dv(1))];
                case 'v0'
                    uv = [du(1) + t*(du(2)-du(1)), dv(1)];
                case 'v1'
                    uv = [du(1) + t*(du(2)-du(1)), dv(2)];
                otherwise
                    error('SurfaceTools:edgeParam01', ...
                        'side must be u0,u1,v0,v1.');
            end
        end


        function D = crossBoundaryDerivativeAlongEdge(S, side, t, outwardSign)
            side = lower(char(string(side)));
            t = t(:);
            D = zeros(numel(t),3);

            for k = 1:numel(t)
                uv = geom.SurfaceTools.edgeParam01(S, side, t(k));
                [Su, Sv] = S.partialDerivatives(uv(1), uv(2));

                switch side
                    case 'v0'
                        d = Sv;
                    case 'v1'
                        d = -Sv;
                    case 'u0'
                        d = Su;
                    case 'u1'
                        d = -Su;
                    otherwise
                        error('SurfaceTools:crossBoundaryDerivativeAlongEdge', ...
                            'side must be u0,u1,v0,v1.');
                end

                D(k,:) = outwardSign * d;
            end
        end


        function u = map01ToDomain(t, dom)
            u = dom(1) + t(:) * (dom(2)-dom(1));
        end


        function t = grevilleParamsForCurve(C)
            n = size(C.P,1) - 1;
            p = C.p;
            U = C.U;

            t = zeros(n+1,1);
            if p == 0
                t(:) = 0;
                return;
            end

            for i = 1:n+1
                t(i) = sum(U(i+1:i+p)) / p;
            end

            dom = C.domain;
            if abs(dom(2)-dom(1)) > eps
                t = (t - dom(1)) / (dom(2)-dom(1));
            else
                t(:) = 0;
            end

            t = max(0,min(1,t));
        end

    end
end
