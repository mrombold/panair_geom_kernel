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
