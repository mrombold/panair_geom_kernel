classdef NURBSSurface < handle
% NURBSSURFACE  Tensor-product rational B-spline surface kernel.
%
% Features
%   - exact rational evaluation
%   - exact rational derivative tensor
%   - analytic first/second partials
%   - normals and curvatures
%   - closest-point projection
%   - refinement, exact splitting, iso-curves, subpatch extraction
%   - Bezier patch decomposition
%   - degree elevation / reduction in u and v
%   - knot removal in u and v
%   - rectangular-net interpolation / least-squares fitting
%   - ruled and bilinear constructors
%   - Coons patch
%   - Gordon surface
%   - loft / skinning with compatible knot merging
%   - surface of revolution
%   - swept surface
%   - trimmed-surface support in UV parameter space
%   - network compatibility helpers

    properties (Access = private)
        Pw_   % [n+1 x m+1 x 4] homogeneous control net [wx wy wz w]
        U_    % [1 x ...] knot vector in u
        V_    % [1 x ...] knot vector in v
        p_    % degree in u
        q_    % degree in v
    end
    
    properties
        trimOuterLoops = {}
        trimInnerLoops = {}
        trimTolerance = 1e-8
    end

    properties (Dependent)
        P
        W
        U
        V
        p
        q
        n
        m
        domainU
        domainV
        Pw
    end

    methods
        function obj = NURBSSurface(P, p, q, U, V, W)
            if nargin < 3
                error('NURBSSurface:Constructor', 'Requires at least P, p, q.');
            end

            if isstruct(P)
                P = cat(3, P.x, P.y, P.z);
            end

            sz = size(P);
            if numel(sz) ~= 3 || sz(3) ~= 3
                error('NURBSSurface:Constructor', 'P must be an (n+1)x(m+1)x3 array.');
            end

            if ~isscalar(p) || p < 0 || p ~= floor(p)
                error('NURBSSurface:Constructor', 'p must be a nonnegative integer.');
            end
            if ~isscalar(q) || q < 0 || q ~= floor(q)
                error('NURBSSurface:Constructor', 'q must be a nonnegative integer.');
            end

            n = sz(1) - 1;
            m = sz(2) - 1;
            if p > n
                error('NURBSSurface:Constructor', 'Degree p=%d exceeds n=%d.', p, n);
            end
            if q > m
                error('NURBSSurface:Constructor', 'Degree q=%d exceeds m=%d.', q, m);
            end

            if nargin < 6 || isempty(W)
                W = ones(n+1, m+1);
            else
                if ~isequal(size(W), [n+1, m+1])
                    error('NURBSSurface:Constructor', 'W must be size (%d x %d).', n+1, m+1);
                end
                if any(W(:) <= 0)
                    error('NURBSSurface:Constructor', 'Weights must be strictly positive.');
                end
            end

            if nargin < 4 || isempty(U)
                U = geom.BasisFunctions.MakeUniformKnotVector(n, p);
            else
                U = U(:).';
            end
            if nargin < 5 || isempty(V)
                V = geom.BasisFunctions.MakeUniformKnotVector(m, q);
            else
                V = V(:).';
            end

            expectedU = n + p + 2;
            expectedV = m + q + 2;
            if numel(U) ~= expectedU
                error('NURBSSurface:Constructor', ...
                    'U length must be n+p+2 = %d, got %d.', expectedU, numel(U));
            end
            if numel(V) ~= expectedV
                error('NURBSSurface:Constructor', ...
                    'V length must be m+q+2 = %d, got %d.', expectedV, numel(V));
            end
            if any(diff(U) < 0)
                error('NURBSSurface:Constructor', 'U must be nondecreasing.');
            end
            if any(diff(V) < 0)
                error('NURBSSurface:Constructor', 'V must be nondecreasing.');
            end

            obj.Pw_ = cat(3, ...
                P(:,:,1).*W, ...
                P(:,:,2).*W, ...
                P(:,:,3).*W, ...
                W);
            obj.U_ = U(:).';
            obj.V_ = V(:).';
            obj.p_ = p;
            obj.q_ = q;
            obj.validate();
        end
        
        function v = get.P(obj)
            w = obj.Pw_(:,:,4);
            v = obj.Pw_(:,:,1:3) ./ w;
        end
        
        function set.P(obj, P)
            if ndims(P) ~= 3 || size(P,3) ~= 3
                error('NURBSSurface:set.P', 'P must be [n+1 x m+1 x 3].');
            end
            if size(P,1) ~= size(obj.Pw_,1) || size(P,2) ~= size(obj.Pw_,2)
                error('NURBSSurface:set.P', 'P size must match existing control net.');
            end
            w = obj.Pw_(:,:,4);
            obj.Pw_ = cat(3, ...
                P(:,:,1).*w, ...
                P(:,:,2).*w, ...
                P(:,:,3).*w, ...
                w);
        end
        
        function v = get.W(obj)
            v = obj.Pw_(:,:,4);
        end
        
        function set.W(obj, W)
            if ~isequal(size(W), [size(obj.Pw_,1), size(obj.Pw_,2)])
                error('NURBSSurface:set.W', 'W must match control-net size.');
            end
            if any(W(:) <= 0)
                error('NURBSSurface:set.W', 'Weights must be strictly positive.');
            end
            P = obj.P;
            obj.Pw_ = cat(3, ...
                P(:,:,1).*W, ...
                P(:,:,2).*W, ...
                P(:,:,3).*W, ...
                W);
        end
        
        function v = get.U(obj)
            v = obj.U_;
        end
        
        function set.U(obj, U)
            U = U(:).';
            if any(diff(U) < 0)
                error('NURBSSurface:set.U', 'U must be nondecreasing.');
            end
            obj.U_ = U;
        end
        
        function v = get.V(obj)
            v = obj.V_;
        end
        
        function set.V(obj, V)
            V = V(:).';
            if any(diff(V) < 0)
                error('NURBSSurface:set.V', 'V must be nondecreasing.');
            end
            obj.V_ = V;
        end
        
        function v = get.p(obj)
            v = obj.p_;
        end
        
        function set.p(obj, p)
            if ~isscalar(p) || p < 0 || p ~= floor(p)
                error('NURBSSurface:set.p', 'p must be a nonnegative integer.');
            end
            obj.p_ = p;
        end
        
        function v = get.q(obj)
            v = obj.q_;
        end
        
        function set.q(obj, q)
            if ~isscalar(q) || q < 0 || q ~= floor(q)
                error('NURBSSurface:set.q', 'q must be a nonnegative integer.');
            end
            obj.q_ = q;
        end
        
        function v = get.n(obj)
            v = size(obj.Pw_,1) - 1;
        end
        
        function v = get.m(obj)
            v = size(obj.Pw_,2) - 1;
        end
        
        function v = get.domainU(obj)
            v = [obj.U_(obj.p_+1), obj.U_(end-obj.p_)];
        end
        
        function v = get.domainV(obj)
            v = [obj.V_(obj.q_+1), obj.V_(end-obj.q_)];
        end
        
        function v = get.Pw(obj)
            v = obj.Pw_;
        end
        
      
        function ok = validate(obj)
            if ndims(obj.Pw_) ~= 3 || size(obj.Pw_,3) ~= 4
                error('NURBSSurface:Validate', 'Pw_ must be [n+1 x m+1 x 4].');
            end
            if numel(obj.U_) ~= obj.n + obj.p_ + 2
                error('NURBSSurface:Validate', 'U length inconsistent with n and p.');
            end
            if numel(obj.V_) ~= obj.m + obj.q_ + 2
                error('NURBSSurface:Validate', 'V length inconsistent with m and q.');
            end
            if any(diff(obj.U_) < 0)
                error('NURBSSurface:Validate', 'U must be nondecreasing.');
            end
            if any(diff(obj.V_) < 0)
                error('NURBSSurface:Validate', 'V must be nondecreasing.');
            end
            if any(obj.Pw_(:,:,4) <= 0, 'all')
                error('NURBSSurface:Validate', 'Weights must be strictly positive.');
            end
            ok = true;
        end

        function C = evaluate(obj, u, v)
            u = u(:);
            v = v(:);
        
            if isscalar(u), u = repmat(u, numel(v), 1); end
            if isscalar(v), v = repmat(v, numel(u), 1); end
        
            assert(numel(u) == numel(v), ...
                'NURBSSurface:evaluate', 'u and v must match in length.');
        
            np = numel(u);
            C = zeros(np, 3);
        
            for k = 1:np
                uk = obj.clampU(u(k));
                vk = obj.clampV(v(k));
        
                uspan = geom.BasisFunctions.FindSpan(obj.n, obj.p, uk, obj.U);
                vspan = geom.BasisFunctions.FindSpan(obj.m, obj.q, vk, obj.V);
        
                Nu = geom.BasisFunctions.BasisFuns(uspan, uk, obj.p, obj.U);
                Nv = geom.BasisFunctions.BasisFuns(vspan, vk, obj.q, obj.V);
        
                Sw = zeros(1, 4);
        
                for l = 0:obj.q
                    temp = zeros(1, 4);
                    jj = vspan - obj.q + l;
        
                    for i = 0:obj.p
                        ii = uspan - obj.p + i;
                        Pw_ij = reshape(obj.Pw_(ii, jj, :), 1, 4);
                        temp = temp + Nu(i+1) * Pw_ij;
                    end
        
                    Sw = Sw + Nv(l+1) * temp;
                end
        
                C(k,:) = Sw(1:3) / Sw(4);
            end
        end

        function C = evaluateTrimmed(obj, u, v)
        % Evaluate only if (u,v) lies inside the active trimmed region.
            if ~obj.isInsideTrim(u, v)
                error('NURBSSurface:evaluateTrimmed', ...
                    'Requested parameter point lies outside trimmed region.');
            end
            C = obj.evaluate(u, v);
        end

        function tf = isTrimmed(obj)
            tf = ~isempty(obj.trimOuterLoops) || ~isempty(obj.trimInnerLoops);
        end

        function obj2 = clearTrims(obj)
            obj2 = obj.copySurfaceOnly();
            obj2.trimOuterLoops = {};
            obj2.trimInnerLoops = {};
            obj2.trimTolerance = obj.trimTolerance;
        end

        function obj2 = setTrims(obj, outerLoops, innerLoops)
            if nargin < 2 || isempty(outerLoops), outerLoops = {}; end
            if nargin < 3 || isempty(innerLoops), innerLoops = {}; end

            obj2 = obj.copySurfaceOnly();
            obj2.trimOuterLoops = geom.NURBSSurface.normalizeTrimLoopSet(outerLoops);
            obj2.trimInnerLoops = geom.NURBSSurface.normalizeTrimLoopSet(innerLoops);
            obj2.trimTolerance = obj.trimTolerance;
        end

        function obj2 = addOuterTrimLoop(obj, loop)
            obj2 = obj.copySurfaceOnly();
            obj2.trimOuterLoops = [obj.trimOuterLoops, ...
                geom.NURBSSurface.normalizeTrimLoopSet({loop})];
            obj2.trimInnerLoops = obj.trimInnerLoops;
            obj2.trimTolerance = obj.trimTolerance;
        end

        function obj2 = addInnerTrimLoop(obj, loop)
            obj2 = obj.copySurfaceOnly();
            obj2.trimOuterLoops = obj.trimOuterLoops;
            obj2.trimInnerLoops = [obj.trimInnerLoops, ...
                geom.NURBSSurface.normalizeTrimLoopSet({loop})];
            obj2.trimTolerance = obj.trimTolerance;
        end

        function tf = isInsideTrim(obj, u, v)
        % True if (u,v) lies in the active kept region of the trimmed surface.
            u = obj.clampU(u);
            v = obj.clampV(v);

            tol = obj.trimTolerance;
            if ~obj.isTrimmed()
                tf = true;
                return;
            end

            % Outer logic
            if isempty(obj.trimOuterLoops)
                insideOuter = true;
            else
                insideOuter = false;
                for k = 1:numel(obj.trimOuterLoops)
                    loop = obj.trimOuterLoops{k};
                    if geom.NURBSSurface.pointInPolygon2D(loop, [u v], tol)
                        insideOuter = true;
                        break;
                    end
                end
            end

            if ~insideOuter
                tf = false;
                return;
            end

            % Inner holes
            insideInner = false;
            for k = 1:numel(obj.trimInnerLoops)
                loop = obj.trimInnerLoops{k};
                if geom.NURBSSurface.pointInPolygon2D(loop, [u v], tol)
                    insideInner = true;
                    break;
                end
            end

            tf = ~insideInner;
        end

        function mask = trimMask(obj, uGrid, vGrid)
            if ~isequal(size(uGrid), size(vGrid))
                error('NURBSSurface:trimMask', 'uGrid and vGrid must have same size.');
            end
            mask = false(size(uGrid));
            for i = 1:size(uGrid,1)
                for j = 1:size(uGrid,2)
                    mask(i,j) = obj.isInsideTrim(uGrid(i,j), vGrid(i,j));
                end
            end
        end

        function UV = sampleTrimLoops(obj, nPerLoop)
            if nargin < 2 || isempty(nPerLoop), nPerLoop = 120; end
            UV.outer = cell(size(obj.trimOuterLoops));
            UV.inner = cell(size(obj.trimInnerLoops));

            for k = 1:numel(obj.trimOuterLoops)
                UV.outer{k} = geom.NURBSSurface.resampleTrimLoop(obj.trimOuterLoops{k}, nPerLoop);
            end
            for k = 1:numel(obj.trimInnerLoops)
                UV.inner{k} = geom.NURBSSurface.resampleTrimLoop(obj.trimInnerLoops{k}, nPerLoop);
            end
        end

        function XYZ = sampleTrimLoopsXYZ(obj, nPerLoop)
            if nargin < 2 || isempty(nPerLoop), nPerLoop = 120; end
            UV = obj.sampleTrimLoops(nPerLoop);
            XYZ.outer = cell(size(UV.outer));
            XYZ.inner = cell(size(UV.inner));

            for k = 1:numel(UV.outer)
                uv = UV.outer{k};
                XYZ.outer{k} = obj.evaluate(uv(:,1), uv(:,2));
            end
            for k = 1:numel(UV.inner)
                uv = UV.inner{k};
                XYZ.inner{k} = obj.evaluate(uv(:,1), uv(:,2));
            end
        end

        function SKL = derivatives(obj, u, v, d)
            if nargin < 4 || isempty(d)
                d = 1;
            end
            if ~isscalar(u) || ~isscalar(v)
                error('NURBSSurface:derivatives', 'Use scalar u,v for derivatives().');
            end

            u = obj.clampU(u);
            v = obj.clampV(v);
            d = max(0, floor(d));

            du = min(d, obj.p);
            dv = min(d, obj.q);

            uspan = geom.BasisFunctions.FindSpan(obj.n, obj.p, u, obj.U);
            vspan = geom.BasisFunctions.FindSpan(obj.m, obj.q, v, obj.V);

            Nu = geom.BasisFunctions.DersBasisFuns(uspan, u, obj.p, du, obj.U);
            Nv = geom.BasisFunctions.DersBasisFuns(vspan, v, obj.q, dv, obj.V);

            Aders = cell(d+1, d+1);
            wders = zeros(d+1, d+1);
            for k = 0:d
                for l = 0:d
                    Aders{k+1,l+1} = zeros(1,3);
                end
            end

            for k = 0:du
                for l = 0:dv
                    Akl = zeros(1,3);
                    wkl = 0.0;
                    for i = 0:obj.p
                        ii = uspan - obj.p + i;
                        for j = 0:obj.q
                            jj = vspan - obj.q + j;
                            B = Nu(k+1, i+1) * Nv(l+1, j+1);
                            Pw_ij = reshape(obj.Pw_(ii, jj, :), 1, 4);
                            Akl = Akl + B * Pw_ij(1:3);
                            wkl = wkl + B * Pw_ij(4);
                        end
                    end
                    Aders{k+1,l+1} = Akl;
                    wders(k+1,l+1) = wkl;
                end
            end
            SKL = geom.NURBSSurface.rationalSurfaceDerivativesFromHomogeneous( ...
                      Aders, wders, d);
        end

    
        function PKL = derivativeControlPoints(obj, d, r1, r2, s1, s2)
            if nargin < 2 || isempty(d), d = 1; end
            if nargin < 3 || isempty(r1), r1 = 1; end
            if nargin < 4 || isempty(r2), r2 = obj.n + 1; end
            if nargin < 5 || isempty(s1), s1 = 1; end
            if nargin < 6 || isempty(s2), s2 = obj.m + 1; end
        
            d = floor(d);
            du = min(d,obj.p);
            dv = min(d,obj.q);
        
            r = r2-r1;
            s = s2-s1;
        
            PKL = cell(d+1,d+1);
            for k = 0:d
                for l = 0:d
                    PKL{k+1,l+1} = zeros(max(r-k+1,0), max(s-l+1,0), 4);
                end
            end
        
            for j = 0:s
                Pcol = squeeze(obj.Pw(r1:r2,s1+j,:));
                PKu = geom.NURBSSurface.curveDerivCpts1D(obj.p,obj.U,Pcol,du,r1,r2);
        
                for k = 0:du
                    PKL{k+1,1}(:,j+1,:) = PKu{k+1};
                end
            end
        
            for k = 0:du
                dd = min(d-k,dv);
                for i = 0:(r-k)
                    Prow = squeeze(PKL{k+1,1}(i+1,1:s+1,:));
                    PKv = geom.NURBSSurface.curveDerivCpts1D(obj.q,obj.V,Prow,dd,s1,s2);
        
                    for l = 1:dd
                        PKL{k+1,l+1}(i+1,:,:) = PKv{l+1};
                    end
                end
            end
        end



            
            
        function Sder = derivativeSurfaceHomogeneous(obj, ku, kv)
        %DERIVATIVESURFACEHOMOGENEOUS Return homogeneous derivative B-spline surface.
        %
        % Sder.Pw    homogeneous derivative control net
        % Sder.U     derivative u knot vector
        % Sder.V     derivative v knot vector
        % Sder.p     derivative u degree
        % Sder.q     derivative v degree
        %
        % This is not a Cartesian NURBS surface for rational derivatives.
        % It is the homogeneous derivative surface.
        
            if nargin < 2, ku = 1; end
            if nargin < 3, kv = 0; end
        
            ku = floor(ku);
            kv = floor(kv);
        
            if ku < 0 || kv < 0 || ku > obj.p || kv > obj.q
                error('NURBSSurface:derivativeSurfaceHomogeneous', ...
                    'Requested derivative order exceeds surface degree.');
            end
        
            d = ku + kv;
            PKL = obj.derivativeControlPoints(d);
        
            Sder = struct();
            Sder.Pw = PKL{ku+1,kv+1};
            Sder.p  = obj.p - ku;
            Sder.q  = obj.q - kv;
            Sder.U  = obj.U(ku+1:end-ku);
            Sder.V  = obj.V(kv+1:end-kv);
        end






        function SKL = derivativesViaControlPoints(obj,u,v,d)
            if nargin < 4 || isempty(d), d = 1; end
        
            u = obj.clampU(u);
            v = obj.clampV(v);
            d = floor(d);
        
            du = min(d,obj.p);
            dv = min(d,obj.q);
        
            uspan = geom.BasisFunctions.FindSpan(obj.n,obj.p,u,obj.U);
            vspan = geom.BasisFunctions.FindSpan(obj.m,obj.q,v,obj.V);
        
            Nu = geom.NURBSSurface.allBasisFuns(uspan,u,obj.p,obj.U);
            Nv = geom.NURBSSurface.allBasisFuns(vspan,v,obj.q,obj.V);
        
            r1 = uspan - obj.p;
            r2 = uspan;
            s1 = vspan - obj.q;
            s2 = vspan;
        
            PKL = obj.derivativeControlPoints(d,r1,r2,s1,s2);
        
            Aders = cell(d+1,d+1);
            wders = zeros(d+1,d+1);
        
            for k = 0:d
                for l = 0:d
                    Aders{k+1,l+1} = zeros(1,3);
                end
            end
        
            for k = 0:du
                dd = min(d-k,dv);
                for l = 0:dd
                    temp = zeros(obj.q-l+1,4);
        
                    for s = 0:(obj.q-l)
                        acc = zeros(1,4);
                        for r = 0:(obj.p-k)
                            acc = acc + Nu(obj.p-k+1,r+1) * ...
                                reshape(PKL{k+1,l+1}(r+1,s+1,:),1,4);
                        end
                        temp(s+1,:) = acc;
                    end
        
                    Sw = zeros(1,4);
                    for s = 0:(obj.q-l)
                        Sw = Sw + Nv(obj.q-l+1,s+1) * temp(s+1,:);
                    end
        
                    Aders{k+1,l+1} = Sw(1:3);
                    wders(k+1,l+1) = Sw(4);
                end
            end
        
            SKL = geom.NURBSSurface.rationalSurfaceDerivativesFromHomogeneous(Aders, wders, d);
        end































        function S2 = reverseU(obj)
            %REVERSEU Reverse the surface parameterization in the u-direction.
            %
            % S2(u,v) = S(a+b-u, v), where [a,b] is the active u knot span.
            aK = obj.U(1);
            bK = obj.U(end);
        
            P2 = flip(obj.P, 1);
            W2 = flip(obj.W, 1);
            U2 = aK + bK - fliplr(obj.U);
        
            S2 = geom.NURBSSurface(P2, obj.p, obj.q, U2, obj.V, W2);
        
            % Preserve trim metadata by reflecting UV trim loops about the active
            % u-domain, not the raw knot span.
            ua = obj.domainU(1);
            ub = obj.domainU(2);
        
            loopsOut = cell(size(obj.trimOuterLoops));
            for k = 1:numel(obj.trimOuterLoops)
                loopIn = obj.trimOuterLoops{k};
                if isa(loopIn, 'geom.NURBSCurve')
                    P = loopIn.P;
                    P(:,1) = ua + ub - P(:,1);
                    loopsOut{k} = geom.NURBSCurve(P, loopIn.p, loopIn.U, loopIn.W);
                elseif isnumeric(loopIn)
                    loopTmp = loopIn;
                    loopTmp(:,1) = ua + ub - loopTmp(:,1);
                    loopsOut{k} = loopTmp;
                else
                    error('NURBSSurface:reverseU', ...
                        'Unsupported trim loop type: %s', class(loopIn));
                end
            end
            S2.trimOuterLoops = loopsOut;
        
            loopsOut = cell(size(obj.trimInnerLoops));
            for k = 1:numel(obj.trimInnerLoops)
                loopIn = obj.trimInnerLoops{k};
                if isa(loopIn, 'geom.NURBSCurve')
                    P = loopIn.P;
                    P(:,1) = ua + ub - P(:,1);
                    loopsOut{k} = geom.NURBSCurve(P, loopIn.p, loopIn.U, loopIn.W);
                elseif isnumeric(loopIn)
                    loopTmp = loopIn;
                    loopTmp(:,1) = ua + ub - loopTmp(:,1);
                    loopsOut{k} = loopTmp;
                else
                    error('NURBSSurface:reverseU', ...
                        'Unsupported trim loop type: %s', class(loopIn));
                end
            end
            S2.trimInnerLoops = loopsOut;
        
            if isprop(obj, 'trimTolerance')
                S2.trimTolerance = obj.trimTolerance;
            end
        end
        
        function S2 = reverseV(obj)
            %REVERSEV Reverse the surface parameterization in the v-direction.
            %
            % S2(u,v) = S(u, c+d-v), where [c,d] is the active v knot span.
            cK = obj.V(1);
            dK = obj.V(end);
        
            P2 = flip(obj.P, 2);
            W2 = flip(obj.W, 2);
            V2 = cK + dK - fliplr(obj.V);
        
            S2 = geom.NURBSSurface(P2, obj.p, obj.q, obj.U, V2, W2);
        
            % Preserve trim metadata by reflecting UV trim loops about the active
            % v-domain, not the raw knot span.
            va = obj.domainV(1);
            vb = obj.domainV(2);
        
            loopsOut = cell(size(obj.trimOuterLoops));
            for k = 1:numel(obj.trimOuterLoops)
                loopIn = obj.trimOuterLoops{k};
                if isa(loopIn, 'geom.NURBSCurve')
                    P = loopIn.P;
                    P(:,2) = va + vb - P(:,2);
                    loopsOut{k} = geom.NURBSCurve(P, loopIn.p, loopIn.U, loopIn.W);
                elseif isnumeric(loopIn)
                    loopTmp = loopIn;
                    loopTmp(:,2) = va + vb - loopTmp(:,2);
                    loopsOut{k} = loopTmp;
                else
                    error('NURBSSurface:reverseV', ...
                        'Unsupported trim loop type: %s', class(loopIn));
                end
            end
            S2.trimOuterLoops = loopsOut;
        
            loopsOut = cell(size(obj.trimInnerLoops));
            for k = 1:numel(obj.trimInnerLoops)
                loopIn = obj.trimInnerLoops{k};
                if isa(loopIn, 'geom.NURBSCurve')
                    P = loopIn.P;
                    P(:,2) = va + vb - P(:,2);
                    loopsOut{k} = geom.NURBSCurve(P, loopIn.p, loopIn.U, loopIn.W);
                elseif isnumeric(loopIn)
                    loopTmp = loopIn;
                    loopTmp(:,2) = va + vb - loopTmp(:,2);
                    loopsOut{k} = loopTmp;
                else
                    error('NURBSSurface:reverseV', ...
                        'Unsupported trim loop type: %s', class(loopIn));
                end
            end
            S2.trimInnerLoops = loopsOut;
        
            if isprop(obj, 'trimTolerance')
                S2.trimTolerance = obj.trimTolerance;
            end
        end
        
        function [S2, maxErr] = reduceU(obj, numTimes, tol, nSample)
            %REDUCEU Alias for reduceDegreeU for API symmetry with elevateU.
            if nargin < 2 || isempty(numTimes), numTimes = 1; end
            if nargin < 3 || isempty(tol),      tol      = inf; end
            if nargin < 4 || isempty(nSample),  nSample  = 200; end
            [S2, maxErr] = obj.reduceDegreeU(numTimes, tol, nSample);
        end
        
        function [S2, maxErr] = reduceV(obj, numTimes, tol, nSample)
            %REDUCEV Alias for reduceDegreeV for API symmetry with elevateV.
            if nargin < 2 || isempty(numTimes), numTimes = 1; end
            if nargin < 3 || isempty(tol),      tol      = inf; end
            if nargin < 4 || isempty(nSample),  nSample  = 200; end
            [S2, maxErr] = obj.reduceDegreeV(numTimes, tol, nSample);
        end

        
        function loopsOut = mapTrimLoopsReverseU(loopsIn, ua, ub)
            loopsOut = cell(size(loopsIn));
            for k = 1:numel(loopsIn)
                loopsOut{k} = geom.NURBSSurface.mapSingleTrimLoopReverseU(loopsIn{k}, ua, ub);
            end
        end
        
        function loopsOut = mapTrimLoopsReverseV(loopsIn, va, vb)
            loopsOut = cell(size(loopsIn));
            for k = 1:numel(loopsIn)
                loopsOut{k} = geom.NURBSSurface.mapSingleTrimLoopReverseV(loopsIn{k}, va, vb);
            end
        end
        
        function loopOut = mapSingleTrimLoopReverseU(loopIn, ua, ub)
            if isa(loopIn, 'geom.NURBSCurve')
                P = loopIn.P;
                P(:,1) = ua + ub - P(:,1);
                loopOut = geom.NURBSCurve(P, loopIn.p, loopIn.U, loopIn.W);
            elseif isnumeric(loopIn)
                loopOut = loopIn;
                loopOut(:,1) = ua + ub - loopOut(:,1);
            else
                error('NURBSSurface:reverseU', ...
                    'Unsupported trim loop type: %s', class(loopIn));
            end
        end
        
        function loopOut = mapSingleTrimLoopReverseV(loopIn, va, vb)
            if isa(loopIn, 'geom.NURBSCurve')
                P = loopIn.P;
                P(:,2) = va + vb - P(:,2);
                loopOut = geom.NURBSCurve(P, loopIn.p, loopIn.U, loopIn.W);
            elseif isnumeric(loopIn)
                loopOut = loopIn;
                loopOut(:,2) = va + vb - loopOut(:,2);
            else
                error('NURBSSurface:reverseV', ...
                    'Unsupported trim loop type: %s', class(loopIn));
            end
        end


        function [Su, Sv] = partialDerivatives(obj, u, v)
            SKL = obj.derivatives(u, v, 1);
            Su = SKL{2,1};
            Sv = SKL{1,2};
        end

        function [S, Su, Sv] = firstPartials(obj, u, v)
            SKL = obj.derivatives(u, v, 1);
            S  = SKL{1,1};
            Su = SKL{2,1};
            Sv = SKL{1,2};
        end

        function [Suu, Suv, Svv] = secondPartials(obj, u, v)
            SKL = obj.derivatives(u, v, 2);
            Suu = SKL{3,1};
            Suv = SKL{2,2};
            Svv = SKL{1,3};
        end

        function [Su, Sv, Suu, Suv, Svv] = firstSecond(obj, u, v)
            SKL = obj.derivatives(u, v, 2);
            Su  = SKL{2,1};
            Sv  = SKL{1,2};
            Suu = SKL{3,1};
            Suv = SKL{2,2};
            Svv = SKL{1,3};
        end

        function N = normal(obj, u, v)
            [Su, Sv] = obj.partialDerivatives(u, v);
            N = cross(Su, Sv);
            nrm = norm(N);
            if nrm > eps
                N = N / nrm;
            else
                N = [0 0 0];
            end
        end

        function [K, H, Nhat] = curvatures(obj, u, v)
            [Su, Sv, Suu, Suv, Svv] = obj.firstSecond(u, v);

            N = cross(Su, Sv);
            nrm = norm(N);
            if nrm < eps
                Nhat = [0 0 0];
                K = NaN;
                H = NaN;
                return;
            end
            Nhat = N / nrm;

            E = dot(Su, Su);
            F = dot(Su, Sv);
            G = dot(Sv, Sv);
            L = dot(Suu, Nhat);
            M = dot(Suv, Nhat);
            N2 = dot(Svv, Nhat);

            den = E * G - F * F;
            if abs(den) < eps
                K = NaN;
                H = NaN;
                return;
            end

            K = (L * N2 - M * M) / den;
            H = (E * N2 - 2 * F * M + G * L) / (2 * den);
        end

        function [u_c, v_c, pt_c, d_c] = closestPoint(obj, P, u0, v0)
            P = P(:).';

            if nargin < 3 || isempty(u0) || nargin < 4 || isempty(v0)
                [u0, v0] = obj.coarseSearchPoint(P, 12, 12);
            end

            u = obj.clampU(u0);
            v = obj.clampV(v0);

            for iter = 1:50
                S = obj.evaluate(u, v);
                [Su, Sv, Suu, Suv, Svv] = obj.firstSecond(u, v);

                r = S - P;
                F = [dot(r, Su); dot(r, Sv)];
                J = [dot(Su, Su) + dot(r, Suu), dot(Su, Sv) + dot(r, Suv); ...
                     dot(Su, Sv) + dot(r, Suv), dot(Sv, Sv) + dot(r, Svv)];

                if rcond(J) < eps
                    break;
                end

                delta = -J \ F;
                u_new = obj.clampU(u + delta(1));
                v_new = obj.clampV(v + delta(2));

                if norm([u_new-u, v_new-v]) < 1e-12
                    u = u_new;
                    v = v_new;
                    break;
                end

                u = u_new;
                v = v_new;
            end

            pt_c = obj.evaluate(u, v);
            d_c = norm(pt_c - P);
            u_c = u;
            v_c = v;
        end

        function [u_all, v_all, pt_all, d_all] = closestPointBatch(obj, pts, nu_c, nv_c)
            if nargin < 3 || isempty(nu_c), nu_c = 16; end
            if nargin < 4 || isempty(nv_c), nv_c = 16; end

            N = size(pts, 1);
            u_all = zeros(N, 1);
            v_all = zeros(N, 1);
            pt_all = zeros(N, 3);
            d_all = zeros(N, 1);

            [u_g, v_g, pts_g] = obj.isoGrid(nu_c, nv_c);
            pts_flat = reshape(pts_g, [], 3);

            [UG, VG] = meshgrid(u_g, v_g);
            u_flat = UG.';
            v_flat = VG.';
            u_flat = u_flat(:);
            v_flat = v_flat(:);

            for k = 1:N
                P = pts(k, :);
                d2 = sum((pts_flat - P).^2, 2);
                [~, idx] = min(d2);
                [u_all(k), v_all(k), pt_all(k,:), d_all(k)] = ...
                    obj.closestPoint(P, u_flat(idx), v_flat(idx));
            end
        end

        function [u_inv, v_inv, pt_inv, err] = invertPoint(obj, P, u0, v0)
            if nargin < 3, u0 = []; end
            if nargin < 4, v0 = []; end
            [u_inv, v_inv, pt_inv, err] = obj.closestPoint(P, u0, v0);
        end
      


        function [u_iso, v_iso, pts] = isoGrid(obj, nu, nv)
            u_iso = linspace(obj.domainU(1), obj.domainU(2), nu);
            v_iso = linspace(obj.domainV(1), obj.domainV(2), nv);
            pts = zeros(nu, nv, 3);

            for i = 1:nu
                for j = 1:nv
                    pts(i, j, :) = obj.evaluate(u_iso(i), v_iso(j));
                end
            end
        end

        function mesh = isoMesh(obj, nu, nv, varargin)
            pa = inputParser;
            addParameter(pa, 'SpacingU', 'uniform');
            addParameter(pa, 'SpacingV', 'uniform');
            addParameter(pa, 'RespectTrim', true);
            parse(pa, varargin{:});
            opts = pa.Results;

            u_iso = geom.NURBSSurface.makeSpacing(obj.domainU(1), obj.domainU(2), nu, opts.SpacingU);
            v_iso = geom.NURBSSurface.makeSpacing(obj.domainV(1), obj.domainV(2), nv, opts.SpacingV);

            pts = zeros(nu, nv, 3);
            for i = 1:nu
                for j = 1:nv
                    pts(i, j, :) = obj.evaluate(u_iso(i), v_iso(j));
                end
            end

            mesh.X = squeeze(pts(:,:,1));
            mesh.Y = squeeze(pts(:,:,2));
            mesh.Z = squeeze(pts(:,:,3));
            mesh.nu = nu;
            mesh.nv = nv;

            [mesh.u, mesh.v] = meshgrid(u_iso, v_iso);
            mesh.u = mesh.u.';
            mesh.v = mesh.v.';

            if opts.RespectTrim && obj.isTrimmed()
                mesh.trimMask = obj.trimMask(mesh.u, mesh.v);
                X = mesh.X; Y = mesh.Y; Z = mesh.Z;
                X(~mesh.trimMask) = NaN;
                Y(~mesh.trimMask) = NaN;
                Z(~mesh.trimMask) = NaN;
                mesh.X = X;
                mesh.Y = Y;
                mesh.Z = Z;
            else
                mesh.trimMask = true(size(mesh.u));
            end

            Normals = zeros(nu, nv, 3);
            for i = 1:nu
                for j = 1:nv
                    if mesh.trimMask(i,j)
                        Normals(i, j, :) = obj.normal(u_iso(i), v_iso(j));
                    else
                        Normals(i, j, :) = [NaN NaN NaN];
                    end
                end
            end
            mesh.normals = Normals;

            nquads = (nu - 1) * (nv - 1);
            conn = zeros(nquads, 4);
            qid = 1;
            for i = 1:nu-1
                for j = 1:nv-1
                    n1 = (i-1) * nv + j;
                    n2 = i * nv + j;
                    n3 = i * nv + j + 1;
                    n4 = (i-1) * nv + j + 1;
                    conn(qid, :) = [n1, n2, n3, n4];
                    qid = qid + 1;
                end
            end
            mesh.connectivity = conn;
        end

        
        function S2 = refine(obj, Xu, Xv)
        %REFINE Batch knot refinement using A5.5-style surface refinement.
        % Xu and Xv may contain repeated knots.
        
            if nargin < 2 || isempty(Xu), Xu = []; end
            if nargin < 3 || isempty(Xv), Xv = []; end
        
            S2 = obj.copySurfaceOnly();
        
            if ~isempty(Xu)
                S2 = S2.refineKnotVectSurfaceU(Xu);
            end
        
            if ~isempty(Xv)
                S2 = S2.refineKnotVectSurfaceV(Xv);
            end
        end

        function S2 = refineKnotVectSurfaceU(obj, X)
        %REFINEKNOTVECTSURFACEU Batch refine surface knot vector in U.
        % Implements The NURBS Book Algorithm A5.5 in the U direction,
        % using homogeneous control points.
        
            X = sort(X(:).');
            if isempty(X)
                S2 = obj.copySurfaceOnly();
                return;
            end
        
            tol = 1e-12;
        
            for kk = 1:numel(X)
                u = X(kk);
                if u < obj.domainU(1)-tol || u > obj.domainU(2)+tol
                    error('NURBSSurface:refineKnotVectSurfaceU', ...
                        'Refinement knot %.16g lies outside active U domain.', u);
                end
        
                sOld = sum(abs(obj.U - u) < tol);
                sAdd = sum(abs(X - u) < tol);
        
                if u > obj.domainU(1)+tol && u < obj.domainU(2)-tol
                    if sOld + sAdd > obj.p
                        error('NURBSSurface:refineKnotVectSurfaceU', ...
                            'Refinement would exceed degree multiplicity at u = %.16g.', u);
                    end
                end
            end
        
            for j = 1:obj.m+1
                PwCol = squeeze(obj.Pw_(:,j,:));
        
                [QwCol, Ubar] = geom.NURBSSurface.curveKnotRefineHomogeneous( ...
                    obj.p, obj.U, PwCol, X);
        
                if j == 1
                    Qw = zeros(size(QwCol,1), obj.m+1, 4);
                end
        
                Qw(:,j,:) = QwCol;
            end
        
            S2 = geom.NURBSSurface.fromHomogeneous(Qw, obj.p, obj.q, Ubar, obj.V);
        
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance  = obj.trimTolerance;
        end
        
        
        function S2 = refineKnotVectSurfaceV(obj, X)
        %REFINEKNOTVECTSURFACEV Batch refine surface knot vector in V.
        % Implements The NURBS Book Algorithm A5.5 in the V direction,
        % using homogeneous control points.
        
            X = sort(X(:).');
            if isempty(X)
                S2 = obj.copySurfaceOnly();
                return;
            end
        
            tol = 1e-12;
        
            for kk = 1:numel(X)
                v = X(kk);
                if v < obj.domainV(1)-tol || v > obj.domainV(2)+tol
                    error('NURBSSurface:refineKnotVectSurfaceV', ...
                        'Refinement knot %.16g lies outside active V domain.', v);
                end
        
                sOld = sum(abs(obj.V - v) < tol);
                sAdd = sum(abs(X - v) < tol);
        
                if v > obj.domainV(1)+tol && v < obj.domainV(2)-tol
                    if sOld + sAdd > obj.q
                        error('NURBSSurface:refineKnotVectSurfaceV', ...
                            'Refinement would exceed degree multiplicity at v = %.16g.', v);
                    end
                end
            end
        
            for i = 1:obj.n+1
                PwRow = squeeze(obj.Pw_(i,:,:));
        
                [QwRow, Vbar] = geom.NURBSSurface.curveKnotRefineHomogeneous( ...
                    obj.q, obj.V, PwRow, X);
        
                if i == 1
                    Qw = zeros(obj.n+1, size(QwRow,1), 4);
                end
        
                Qw(i,:,:) = QwRow;
            end
        
            S2 = geom.NURBSSurface.fromHomogeneous(Qw, obj.p, obj.q, obj.U, Vbar);
        
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance  = obj.trimTolerance;
        end
        
        
        
        function S2 = flipNormals(obj)
            P2 = flipud(obj.P);
            W2 = flipud(obj.W);
            U2 = (obj.U(1) + obj.U(end)) - fliplr(obj.U);
            S2 = geom.NURBSSurface(P2, obj.p, obj.q, U2, obj.V, W2);
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance = obj.trimTolerance;
        end


        function [S_lo, S_hi] = splitU(obj, u0)
        %SPLITU Exact surface split in U by full-multiplicity knot insertion.
        
            tol = 1e-12;
            u0 = obj.clampU(u0);
        
            if u0 <= obj.domainU(1) + tol || u0 >= obj.domainU(2) - tol
                error('NURBSSurface:splitU', ...
                    'u0 must lie strictly inside the active U domain.');
            end
        
            s = sum(abs(obj.U - u0) < tol);
            if s > obj.p
                error('NURBSSurface:splitU', ...
                    'Split knot multiplicity exceeds degree.');
            end
        
            reps = obj.p - s;
            if reps > 0
                S_ref = obj.refine(repmat(u0, 1, reps), []);
            else
                S_ref = obj.copySurfaceOnly();
            end
        
            U2 = S_ref.U;
            idx = find(abs(U2 - u0) < tol);
        
            if numel(idx) < obj.p
                error('NURBSSurface:splitU', ...
                    'Failed to create full-multiplicity split knot.');
            end
        
            first = idx(1);
            last  = idx(end);
        
            % splitU
            nLo = last - obj.p;
            
            P_lo = S_ref.P(1:nLo, :, :);
            W_lo = S_ref.W(1:nLo, :);
            U_lo = [U2(1:last), u0];
            
            P_hi = S_ref.P(nLo:end, :, :);
            W_hi = S_ref.W(nLo:end, :);
            U_hi = [u0, U2(first:end)];
        
            % Clamp child knot vectors at the split while keeping absolute domains.
            U_lo(end-obj.p:end) = u0;
            U_hi(1:obj.p+1) = u0;
        
            S_lo = geom.NURBSSurface(P_lo, obj.p, obj.q, U_lo, obj.V, W_lo);
            S_hi = geom.NURBSSurface(P_hi, obj.p, obj.q, U_hi, obj.V, W_hi);
        
            % Conservative trim handling: preserve existing trim metadata.
            % This does not do topological trim splitting.
            S_lo.trimOuterLoops = obj.trimOuterLoops;
            S_lo.trimInnerLoops = obj.trimInnerLoops;
            S_lo.trimTolerance  = obj.trimTolerance;
        
            S_hi.trimOuterLoops = obj.trimOuterLoops;
            S_hi.trimInnerLoops = obj.trimInnerLoops;
            S_hi.trimTolerance  = obj.trimTolerance;
        end

        function [S_lo, S_hi] = splitV(obj, v0)
        %SPLITV Exact surface split in V by full-multiplicity knot insertion.
        
            tol = 1e-12;
            v0 = obj.clampV(v0);
        
            if v0 <= obj.domainV(1) + tol || v0 >= obj.domainV(2) - tol
                error('NURBSSurface:splitV', ...
                    'v0 must lie strictly inside the active V domain.');
            end
        
            s = sum(abs(obj.V - v0) < tol);
            if s > obj.q
                error('NURBSSurface:splitV', ...
                    'Split knot multiplicity exceeds degree.');
            end
        
            reps = obj.q - s;
            if reps > 0
                S_ref = obj.refine([], repmat(v0, 1, reps));
            else
                S_ref = obj.copySurfaceOnly();
            end
        
            V2 = S_ref.V;
            idx = find(abs(V2 - v0) < tol);
        
            if numel(idx) < obj.q
                error('NURBSSurface:splitV', ...
                    'Failed to create full-multiplicity split knot.');
            end
        
            first = idx(1);
            last  = idx(end);
        
            % Number of control points on low side.
            % splitV
            mLo = last - obj.q;
            
            P_lo = S_ref.P(:, 1:mLo, :);
            W_lo = S_ref.W(:, 1:mLo);
            V_lo = [V2(1:last), v0];
            
            P_hi = S_ref.P(:, mLo:end, :);
            W_hi = S_ref.W(:, mLo:end);
            V_hi = [v0, V2(first:end)];
        
            % Clamp child knot vectors at the split while keeping absolute domains.
            V_lo(end-obj.q:end) = v0;
            V_hi(1:obj.q+1) = v0;
        
            S_lo = geom.NURBSSurface(P_lo, obj.p, obj.q, obj.U, V_lo, W_lo);
            S_hi = geom.NURBSSurface(P_hi, obj.p, obj.q, obj.U, V_hi, W_hi);
        
            % Conservative trim handling: preserve existing trim metadata.
            % This does not do topological trim splitting.
            S_lo.trimOuterLoops = obj.trimOuterLoops;
            S_lo.trimInnerLoops = obj.trimInnerLoops;
            S_lo.trimTolerance  = obj.trimTolerance;
        
            S_hi.trimOuterLoops = obj.trimOuterLoops;
            S_hi.trimInnerLoops = obj.trimInnerLoops;
            S_hi.trimTolerance  = obj.trimTolerance;
        end

        function C = isoCurveU(obj, v0)
            v0 = obj.clampV(v0);
        
            vspan = geom.BasisFunctions.FindSpan(obj.m, obj.q, v0, obj.V);
            Nv = geom.BasisFunctions.BasisFuns(vspan, v0, obj.q, obj.V);
        
            Pw = zeros(obj.n + 1, 4);
        
            for i = 1:obj.n+1
                for l = 0:obj.q
                    jj = vspan - obj.q + l;
                    Pw_ij = reshape(obj.Pw_(i, jj, :), 1, 4);
                    Pw(i,:) = Pw(i,:) + Nv(l+1) * Pw_ij;
                end
            end
        
            W2 = Pw(:,4);
            P2 = Pw(:,1:3) ./ W2;
        
            C = geom.NURBSCurve(P2, obj.p, obj.U, W2);
        end

        function C = isoCurveV(obj, u0)
            u0 = obj.clampU(u0);
        
            uspan = geom.BasisFunctions.FindSpan(obj.n, obj.p, u0, obj.U);
            Nu = geom.BasisFunctions.BasisFuns(uspan, u0, obj.p, obj.U);
        
            Pw = zeros(obj.m + 1, 4);
        
            for j = 1:obj.m+1
                for l = 0:obj.p
                    ii = uspan - obj.p + l;
                    Pw_ij = reshape(obj.Pw_(ii, j, :), 1, 4);
                    Pw(j,:) = Pw(j,:) + Nu(l+1) * Pw_ij;
                end
            end
        
            W2 = Pw(:,4);
            P2 = Pw(:,1:3) ./ W2;
        
            C = geom.NURBSCurve(P2, obj.q, obj.V, W2);
        end


        function Ssub = extractUV(obj, urange, vrange)
        %EXTRACTUV Exact rectangular subpatch extraction using absolute UV splits.
        
            ua = urange(1);
            ub = urange(2);
            va = vrange(1);
            vb = vrange(2);
        
            if ~(ua < ub && va < vb)
                error('NURBSSurface:extractUV', ...
                    'Require ua < ub and va < vb.');
            end
        
            tol = 1e-12;
        
            ua = obj.clampU(ua);
            ub = obj.clampU(ub);
            va = obj.clampV(va);
            vb = obj.clampV(vb);
        
            Swork = obj.copySurfaceOnly();
        
            if ua > Swork.domainU(1) + tol
                [~, Swork] = Swork.splitU(ua);
            end
        
            if ub < Swork.domainU(2) - tol
                [Swork, ~] = Swork.splitU(ub);
            end
        
            if va > Swork.domainV(1) + tol
                [~, Swork] = Swork.splitV(va);
            end
        
            if vb < Swork.domainV(2) - tol
                [Swork, ~] = Swork.splitV(vb);
            end
        
            Ssub = Swork;
        end

        function B = boundaryCurves(obj)
            B.u0 = obj.isoCurveV(obj.domainU(1));
            B.u1 = obj.isoCurveV(obj.domainU(2));
            B.v0 = obj.isoCurveU(obj.domainV(1));
            B.v1 = obj.isoCurveU(obj.domainV(2));
        end

        function C = extractInnerSpanCurve(obj)
            u_pick = obj.domainU(1) + 0.05 * (obj.domainU(2) - obj.domainU(1));
            C = obj.isoCurveV(u_pick);
        end

        function C = extractInnerCircCurve(obj)
            v_pick = obj.domainV(1) + 0.05 * (obj.domainV(2) - obj.domainV(1));
            C = obj.isoCurveU(v_pick);
        end

        function St = swapUV(obj)
            P2 = permute(obj.P, [2 1 3]);
            W2 = obj.W.';
            St = geom.NURBSSurface(P2, obj.q, obj.p, obj.V, obj.U, W2);
            St.trimOuterLoops = obj.trimOuterLoops;
            St.trimInnerLoops = obj.trimInnerLoops;
            St.trimTolerance = obj.trimTolerance;
        end

        function patches = decomposeBezier(obj)
            Su = obj;

            Ui = unique(Su.U(Su.p+2:end-Su.p-1));
            for k = 1:numel(Ui)
                s = sum(abs(Su.U - Ui(k)) < 1e-12);
                if s < Su.p
                    Su = Su.refine(repmat(Ui(k), 1, Su.p - s), []);
                end
            end

            Vi = unique(Su.V(Su.q+2:end-Su.q-1));
            for k = 1:numel(Vi)
                s = sum(abs(Su.V - Vi(k)) < 1e-12);
                if s < Su.q
                    Su = Su.refine([], repmat(Vi(k), 1, Su.q - s));
                end
            end

            uBreaks = unique(Su.U);
            vBreaks = unique(Su.V);

            nuSeg = numel(uBreaks) - 1;
            nvSeg = numel(vBreaks) - 1;

            patches = cell(nuSeg, nvSeg);

            for iu = 1:nuSeg
                for jv = 1:nvSeg
                    i0 = (iu-1)*Su.p + 1;
                    i1 = i0 + Su.p;
                    j0 = (jv-1)*Su.q + 1;
                    j1 = j0 + Su.q;

                    Pij = Su.P(i0:i1, j0:j1, :);
                    Wij = Su.W(i0:i1, j0:j1);
                    Uij = [zeros(1, Su.p+1), ones(1, Su.p+1)];
                    Vij = [zeros(1, Su.q+1), ones(1, Su.q+1)];

                    patches{iu,jv} = struct( ...
                        'surface', geom.NURBSSurface(Pij, Su.p, Su.q, Uij, Vij, Wij), ...
                        'u0', uBreaks(iu), 'u1', uBreaks(iu+1), ...
                        'v0', vBreaks(jv), 'v1', vBreaks(jv+1));
                end
            end
        end

        function S2 = elevateU(obj, t)
            if nargin < 2 || isempty(t), t = 1; end
            if t < 0 || t ~= floor(t)
                error('NURBSSurface:elevateU', 't must be a nonnegative integer.');
            end
            if t == 0
                S2 = obj.copySurfaceOnly();
                return;
            end

            for j = 1:size(obj.P,2)
                Cj = geom.NURBSCurve(squeeze(obj.P(:,j,:)), obj.p, obj.U, obj.W(:,j));
                Ce = Cj.elevate(t);
                if j == 1
                    nnew = size(Ce.P,1);
                    P2 = zeros(nnew, size(obj.P,2), 3);
                    W2 = zeros(nnew, size(obj.P,2));
                    U2 = Ce.U;
                end
                P2(:,j,:) = Ce.P;
                W2(:,j) = Ce.W;
            end

            S2 = geom.NURBSSurface(P2, obj.p + t, obj.q, U2, obj.V, W2);
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance = obj.trimTolerance;
        end

        function S2 = elevateV(obj, t)
            if nargin < 2 || isempty(t), t = 1; end
            if t < 0 || t ~= floor(t)
                error('NURBSSurface:elevateV', 't must be a nonnegative integer.');
            end
            if t == 0
                S2 = obj.copySurfaceOnly();
                return;
            end

            for i = 1:size(obj.P,1)
                Ci = geom.NURBSCurve(squeeze(obj.P(i,:,:)), obj.q, obj.V, obj.W(i,:).');
                Ce = Ci.elevate(t);
                if i == 1
                    mnew = size(Ce.P,1);
                    P2 = zeros(size(obj.P,1), mnew, 3);
                    W2 = zeros(size(obj.P,1), mnew);
                    V2 = Ce.U;
                end
                P2(i,:,:) = Ce.P;
                W2(i,:) = Ce.W.';
            end

            S2 = geom.NURBSSurface(P2, obj.p, obj.q + t, obj.U, V2, W2);
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance = obj.trimTolerance;
        end

        function [S2, maxErr] = reduceDegreeU(obj, numTimes, tol, nSample)
            if nargin < 2 || isempty(numTimes), numTimes = 1; end
            if nargin < 3 || isempty(tol), tol = inf; end
            if nargin < 4 || isempty(nSample), nSample = 200; end

            maxErr = 0;
            Scur = obj;

            for step = 1:numTimes
                if Scur.p <= 0
                    error('NURBSSurface:reduceDegreeU', 'Cannot reduce p below zero.');
                end

                for j = 1:size(Scur.P,2)
                    Cj = geom.NURBSCurve(squeeze(Scur.P(:,j,:)), Scur.p, Scur.U, Scur.W(:,j));
                    [Cr, errj] = Cj.reduceDegree(1, tol, nSample);
                    if j == 1
                        nnew = size(Cr.P,1);
                        P2 = zeros(nnew, size(Scur.P,2), 3);
                        W2 = zeros(nnew, size(Scur.P,2));
                        U2 = Cr.U;
                    end
                    P2(:,j,:) = Cr.P;
                    W2(:,j) = Cr.W;
                    maxErr = max(maxErr, errj);
                end

                Scur = geom.NURBSSurface(P2, Scur.p - 1, Scur.q, U2, Scur.V, W2);
            end

            S2 = Scur;
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance = obj.trimTolerance;
        end

        function [S2, maxErr] = reduceDegreeV(obj, numTimes, tol, nSample)
            if nargin < 2 || isempty(numTimes), numTimes = 1; end
            if nargin < 3 || isempty(tol), tol = inf; end
            if nargin < 4 || isempty(nSample), nSample = 200; end

            maxErr = 0;
            Scur = obj;

            for step = 1:numTimes
                if Scur.q <= 0
                    error('NURBSSurface:reduceDegreeV', 'Cannot reduce q below zero.');
                end

                for i = 1:size(Scur.P,1)
                    Ci = geom.NURBSCurve(squeeze(Scur.P(i,:,:)), Scur.q, Scur.V, Scur.W(i,:).');
                    [Cr, erri] = Ci.reduceDegree(1, tol, nSample);
                    if i == 1
                        mnew = size(Cr.P,1);
                        P2 = zeros(size(Scur.P,1), mnew, 3);
                        W2 = zeros(size(Scur.P,1), mnew);
                        V2 = Cr.U;
                    end
                    P2(i,:,:) = Cr.P;
                    W2(i,:) = Cr.W.';
                    maxErr = max(maxErr, erri);
                end

                Scur = geom.NURBSSurface(P2, Scur.p, Scur.q - 1, Scur.U, V2, W2);
            end

            S2 = Scur;
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance = obj.trimTolerance;
        end

        function [S2, removed, maxErr] = removeKnotU(obj, u, numRemove, tol, nSample)
            if nargin < 3 || isempty(numRemove), numRemove = 1; end
            if nargin < 4 || isempty(tol), tol = 1e-10; end
            if nargin < 5 || isempty(nSample), nSample = 200; end

            Scur = obj;
            removed = inf;
            maxErr = 0;

            for j = 1:size(Scur.P,2)
                Cj = geom.NURBSCurve(squeeze(Scur.P(:,j,:)), Scur.p, Scur.U, Scur.W(:,j));
                [Cr, remj, errj] = Cj.removeKnot(u, numRemove, tol, nSample);

                if j == 1
                    removed = remj;
                    nnew = size(Cr.P,1);
                    P2 = zeros(nnew, size(Scur.P,2), 3);
                    W2 = zeros(nnew, size(Scur.P,2));
                    U2 = Cr.U;
                else
                    removed = min(removed, remj);
                end

                P2(:,j,:) = Cr.P;
                W2(:,j) = Cr.W;
                maxErr = max(maxErr, errj);
            end

            S2 = geom.NURBSSurface(P2, Scur.p, Scur.q, U2, Scur.V, W2);
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance = obj.trimTolerance;
        end

        function [S2, removed, maxErr] = removeKnotV(obj, v, numRemove, tol, nSample)
            if nargin < 3 || isempty(numRemove), numRemove = 1; end
            if nargin < 4 || isempty(tol), tol = 1e-10; end
            if nargin < 5 || isempty(nSample), nSample = 200; end

            Scur = obj;
            removed = inf;
            maxErr = 0;

            for i = 1:size(Scur.P,1)
                Ci = geom.NURBSCurve(squeeze(Scur.P(i,:,:)), Scur.q, Scur.V, Scur.W(i,:).');
                [Cr, remi, erri] = Ci.removeKnot(v, numRemove, tol, nSample);

                if i == 1
                    removed = remi;
                    mnew = size(Cr.P,1);
                    P2 = zeros(size(Scur.P,1), mnew, 3);
                    W2 = zeros(size(Scur.P,1), mnew);
                    V2 = Cr.U;
                else
                    removed = min(removed, remi);
                end

                P2(i,:,:) = Cr.P;
                W2(i,:) = Cr.W.';
                maxErr = max(maxErr, erri);
            end

            S2 = geom.NURBSSurface(P2, Scur.p, Scur.q, Scur.U, V2, W2);
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance = obj.trimTolerance;
        end



        function plot(obj, nu, nv, varargin)
            if nargin < 2, nu = 30; end
            if nargin < 3, nv = 30; end

            pa = inputParser;
            addParameter(pa, 'ShowCP', false);
            addParameter(pa, 'ShowIso', false);
            addParameter(pa, 'ShowTrims', true);
            addParameter(pa, 'Alpha', 0.85);
            addParameter(pa, 'EdgeAlpha', 0.2);
            addParameter(pa, 'FaceColor', [0.3 0.6 0.9]);
            parse(pa, varargin{:});
            opts = pa.Results;

            mesh = obj.isoMesh(nu, nv, 'RespectTrim', true);
            surf(mesh.X, mesh.Y, mesh.Z, ...
                'FaceColor', opts.FaceColor, ...
                'FaceAlpha', opts.Alpha, ...
                'EdgeAlpha', opts.EdgeAlpha, ...
                'EdgeColor', 'k');
            hold on;

            if opts.ShowIso
                n_iso = 8;
                u_lines = linspace(obj.domainU(1), obj.domainU(2), n_iso + 2);
                u_lines = u_lines(2:end-1);
                v_lines = linspace(obj.domainV(1), obj.domainV(2), n_iso + 2);
                v_lines = v_lines(2:end-1);

                for i = 1:numel(u_lines)
                    C = obj.isoCurveV(u_lines(i));
                    C.plot(50, 'ShowCP', false, 'Color', [0.1 0.1 0.5], 'LineWidth', 0.7);
                end
                for j = 1:numel(v_lines)
                    C = obj.isoCurveU(v_lines(j));
                    C.plot(50, 'ShowCP', false, 'Color', [0.5 0.1 0.1], 'LineWidth', 0.7);
                end
            end

            if opts.ShowCP
                for i = 1:size(obj.P, 1)
                    plot3(squeeze(obj.P(i,:,1)), squeeze(obj.P(i,:,2)), squeeze(obj.P(i,:,3)), ...
                        '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
                end
                for j = 1:size(obj.P, 2)
                    plot3(squeeze(obj.P(:,j,1)), squeeze(obj.P(:,j,2)), squeeze(obj.P(:,j,3)), ...
                        '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
                end
                Pf = reshape(obj.P, [], 3);
                plot3(Pf(:,1), Pf(:,2), Pf(:,3), 'o', ...
                    'MarkerSize', 4, 'MarkerFaceColor', 'w', 'Color', [0.3 0.3 0.3]);
            end

            if opts.ShowTrims && obj.isTrimmed()
                T = obj.sampleTrimLoopsXYZ(160);
                for k = 1:numel(T.outer)
                    xyz = T.outer{k};
                    plot3(xyz(:,1), xyz(:,2), xyz(:,3), '-', ...
                        'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
                end
                for k = 1:numel(T.inner)
                    xyz = T.inner{k};
                    plot3(xyz(:,1), xyz(:,2), xyz(:,3), '-', ...
                        'Color', [0.1 0.1 0.9], 'LineWidth', 2.0);
                end
            end

            axis equal;
            grid on;
            lighting phong;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
        end

        function plotNormals(obj, nu, nv, scale)
            if nargin < 2, nu = 8; end
            if nargin < 3, nv = 8; end
            if nargin < 4, scale = 0.05; end

            mesh = obj.isoMesh(nu, nv, 'RespectTrim', true);
            quiver3(mesh.X, mesh.Y, mesh.Z, ...
                scale * mesh.normals(:,:,1), ...
                scale * mesh.normals(:,:,2), ...
                scale * mesh.normals(:,:,3), ...
                0, 'b');
        end

        function plotTrimUV(obj, nPerLoop)
            if nargin < 2 || isempty(nPerLoop), nPerLoop = 160; end
            UV = obj.sampleTrimLoops(nPerLoop);
            hold on; grid on; axis equal;
            xlabel('u'); ylabel('v');
            xlim(obj.domainU);
            ylim(obj.domainV);

            rectangle('Position', [obj.domainU(1), obj.domainV(1), ...
                                   obj.domainU(2)-obj.domainU(1), ...
                                   obj.domainV(2)-obj.domainV(1)], ...
                      'EdgeColor', [0.5 0.5 0.5], 'LineStyle', '--');

            for k = 1:numel(UV.outer)
                uv = UV.outer{k};
                plot(uv(:,1), uv(:,2), 'r-', 'LineWidth', 2.0);
            end
            for k = 1:numel(UV.inner)
                uv = UV.inner{k};
                plot(uv(:,1), uv(:,2), 'b-', 'LineWidth', 2.0);
            end
            title('Trim loops in UV space');
        end
    

        function obj2 = copySurfaceOnly(obj)
            obj2 = geom.NURBSSurface.fromHomogeneous( ...
                obj.Pw_, obj.p, obj.q, obj.U, obj.V);
        
            obj2.trimOuterLoops = obj.trimOuterLoops;
            obj2.trimInnerLoops = obj.trimInnerLoops;
            obj2.trimTolerance = obj.trimTolerance;
        end

        function S2 = insertKnotU(obj, u, r)
        %INSERTKNOTU Exact surface knot insertion in u-direction.
        % Implements The NURBS Book Algorithm A5.3 by applying homogeneous
        % curve knot insertion to every v-column.
        
            if nargin < 3 || isempty(r), r = 1; end
            if r < 0 || r ~= floor(r)
                error('NURBSSurface:insertKnotU', 'r must be a nonnegative integer.');
            end
            if r == 0
                S2 = obj.copySurfaceOnly();
                return;
            end
        
            u = obj.clampU(u);
        
            tol = 1e-12;
            s = sum(abs(obj.U - u) < tol);
        
            if u > obj.domainU(1) + tol && u < obj.domainU(2) - tol
                if s + r > obj.p
                    error('NURBSSurface:insertKnotU', ...
                        'Cannot insert knot %g %d times: multiplicity would exceed degree.', u, r);
                end
            end
        
            for j = 1:obj.m+1
                PwCol = squeeze(obj.Pw(:,j,:));
                [QwCol, Uq] = geom.NURBSSurface.curveKnotInsertHomogeneous( ...
                    obj.p, obj.U, PwCol, u, r);
        
                if j == 1
                    Qw = zeros(size(QwCol,1), obj.m+1, 4);
                end
        
                Qw(:,j,:) = QwCol;
            end
        
            S2 = geom.NURBSSurface.fromHomogeneous(Qw, obj.p, obj.q, Uq, obj.V);
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance = obj.trimTolerance;
        end
        
        function S2 = insertKnotV(obj, v, r)
        %INSERTKNOTV Exact surface knot insertion in v-direction.
        % Implements The NURBS Book Algorithm A5.3 in v-direction.
        
            if nargin < 3 || isempty(r), r = 1; end
            if r < 0 || r ~= floor(r)
                error('NURBSSurface:insertKnotV', 'r must be a nonnegative integer.');
            end
            if r == 0
                S2 = obj.copySurfaceOnly();
                return;
            end
        
            v = obj.clampV(v);
        
            tol = 1e-12;
            s = sum(abs(obj.V - v) < tol);
        
            if v > obj.domainV(1) + tol && v < obj.domainV(2) - tol
                if s + r > obj.q
                    error('NURBSSurface:insertKnotV', ...
                        'Cannot insert knot %g %d times: multiplicity would exceed degree.', v, r);
                end
            end
        
            for i = 1:obj.n+1
                PwRow = squeeze(obj.Pw(i,:,:));
                [QwRow, Vq] = geom.NURBSSurface.curveKnotInsertHomogeneous( ...
                    obj.q, obj.V, PwRow, v, r);
        
                if i == 1
                    Qw = zeros(obj.n+1, size(QwRow,1), 4);
                end
        
                Qw(i,:,:) = QwRow;
            end
        
            S2 = geom.NURBSSurface.fromHomogeneous(Qw, obj.p, obj.q, obj.U, Vq);
            S2.trimOuterLoops = obj.trimOuterLoops;
            S2.trimInnerLoops = obj.trimInnerLoops;
            S2.trimTolerance = obj.trimTolerance;
        end
    end


    methods (Access = private)
        function u = clampU(obj, u)
            u = max(obj.domainU(1), min(obj.domainU(2), u));
        end

        function v = clampV(obj, v)
            v = max(obj.domainV(1), min(obj.domainV(2), v));
        end

        function [u0, v0] = coarseSearchPoint(obj, P, nu, nv)
            u_c = linspace(obj.domainU(1), obj.domainU(2), nu);
            v_c = linspace(obj.domainV(1), obj.domainV(2), nv);

            best_d2 = inf;
            u0 = u_c(1);
            v0 = v_c(1);

            for i = 1:nu
                for j = 1:nv
                    pt = obj.evaluate(u_c(i), v_c(j));
                    d2 = sum((pt - P).^2);
                    if d2 < best_d2
                        best_d2 = d2;
                        u0 = u_c(i);
                        v0 = v_c(j);
                    end
                end
            end
        end
    end




    
    methods (Static)

        function S = fromHomogeneous(Pw, p, q, U, V)
            if ndims(Pw) ~= 3 || size(Pw,3) ~= 4
                error('NURBSSurface:fromHomogeneous', ...
                    'Pw must be [n+1 x m+1 x 4].');
            end
        
            W = Pw(:,:,4);
        
            if any(W(:) <= 0)
                error('NURBSSurface:fromHomogeneous', ...
                    'Weights must be strictly positive.');
            end
        
            P = Pw(:,:,1:3) ./ W;
        
            S = geom.NURBSSurface(P, p, q, U, V, W);
        end

        function t = makeSpacing(t0, t1, n, type)
            switch lower(type)
                case 'cosine'
                    beta = linspace(0, pi, n);
                    s = 0.5 * (1 - cos(beta));
                case 'cosine_half'
                    beta = linspace(0, pi/2, n);
                    s = sin(beta);
                case 'cosine_wrap'
                    n1 = ceil(n/2);
                    n2 = n - n1 + 1;
                    b1 = linspace(0, pi, n1);
                    s1 = 0.5 * (1 - cos(b1)) * 0.5;
                    b2 = linspace(0, pi, n2);
                    s2 = 0.5 + 0.5 * (1 - cos(b2)) * 0.5;
                    s = [s1, s2(2:end)];
                otherwise
                    s = linspace(0, 1, n);
            end
            t = t0 + s * (t1 - t0);
        end

        function U = normalizeKnotVector(U)
            a = U(1);
            b = U(end);
            if abs(b-a) < eps
                U = zeros(size(U));
            else
                U = (U - a) / (b - a);
            end
        end

        function loops = normalizeTrimLoopSet(loopsIn)
            if isempty(loopsIn)
                loops = {};
                return;
            end
            if ~iscell(loopsIn)
                loopsIn = {loopsIn};
            end
            loops = cell(size(loopsIn));
            for k = 1:numel(loopsIn)
                loops{k} = geom.NURBSSurface.loopToPolylineUV(loopsIn{k}, 200);
            end
        end

        function uv = loopToPolylineUV(loop, nSamp)
            if nargin < 2 || isempty(nSamp), nSamp = 200; end

            if isa(loop, 'geom.NURBSCurve')
                u = linspace(loop.domain(1), loop.domain(2), nSamp);
                xyz = loop.evaluate(u);
                uv = xyz(:,1:2);
            elseif isnumeric(loop) && size(loop,2) == 2
                uv = loop;
            else
                error('NURBSSurface:loopToPolylineUV', ...
                    'Trim loop must be Nx2 numeric or geom.NURBSCurve in UV plane.');
            end

            if size(uv,1) < 3
                error('NURBSSurface:loopToPolylineUV', ...
                    'Trim loop needs at least 3 points.');
            end

            if norm(uv(1,:) - uv(end,:)) > 1e-12
                uv = [uv; uv(1,:)];
            end
        end

        function tf = pointInPolygon2D(poly, pt, tol)
            if nargin < 3 || isempty(tol), tol = 1e-10; end
            x = pt(1);
            y = pt(2);
            xv = poly(:,1);
            yv = poly(:,2);

            % Boundary check
            for i = 1:size(poly,1)-1
                a = poly(i,:);
                b = poly(i+1,:);
                if geom.NURBSSurface.pointOnSegment2D(pt, a, b, tol)
                    tf = true;
                    return;
                end
            end

            tf = inpolygon(x, y, xv, yv);
        end

        function tf = pointOnSegment2D(p, a, b, tol)
            ap = p - a;
            ab = b - a;
            lab2 = dot(ab, ab);
            if lab2 < eps
                tf = norm(ap) <= tol;
                return;
            end
            t = dot(ap, ab) / lab2;
            if t < -tol || t > 1+tol
                tf = false;
                return;
            end
            proj = a + min(max(t,0),1) * ab;
            tf = norm(p - proj) <= tol;
        end

        function uv = resampleTrimLoop(loop, nPerLoop)
            if isa(loop, 'geom.NURBSCurve')
                us = linspace(loop.domain(1), loop.domain(2), nPerLoop);
                xyz = loop.evaluate(us);
                uv = xyz(:,1:2);
                if norm(uv(1,:) - uv(end,:)) > 1e-12
                    uv = [uv; uv(1,:)];
                end
            else
                uv0 = geom.NURBSSurface.loopToPolylineUV(loop, nPerLoop);
                d = [0; cumsum(vecnorm(diff(uv0,1,1),2,2))];
                if d(end) < eps
                    uv = uv0;
                    return;
                end
                s = linspace(0, d(end), nPerLoop).';
                uv = [interp1(d, uv0(:,1), s, 'linear'), ...
                      interp1(d, uv0(:,2), s, 'linear')];
                if norm(uv(1,:) - uv(end,:)) > 1e-12
                    uv = [uv; uv(1,:)];
                end
            end
        end


        function S = globalInterpNet(Q, p, q, methodU, methodV)
            if nargin < 4 || isempty(methodU), methodU = 'centripetal'; end
            if nargin < 5 || isempty(methodV), methodV = 'centripetal'; end

            if size(Q,3) ~= 3
                error('NURBSSurface:globalInterpNet', 'Q must be nu x nv x 3.');
            end

            nv = size(Q,2);

            for j = 1:nv
                Cj = geom.NURBSCurve.globalInterp(squeeze(Q(:,j,:)), p, methodU);
                if j == 1
                    Pu = zeros(size(Cj.P,1), nv, 3);
                    Wu = zeros(size(Cj.P,1), nv);
                    U = Cj.U;
                end
                Pu(:,j,:) = Cj.P;
                Wu(:,j) = Cj.W;
            end

            for i = 1:size(Pu,1)
                Ci = geom.NURBSCurve.globalInterp(squeeze(Pu(i,:,:)), q, methodV);
                if i == 1
                    P = zeros(size(Pu,1), size(Ci.P,1), 3);
                    W = zeros(size(Pu,1), size(Ci.P,1));
                    V = Ci.U;
                end
                P(i,:,:) = Ci.P;
                W(i,:) = Ci.W.';
            end

            S = geom.NURBSSurface(P, p, q, U, V, W);
        end

        function S = globalLeastSquaresFitNet(Q, p, q, nCtrlU, nCtrlV, methodU, methodV)
            if nargin < 6 || isempty(methodU), methodU = 'centripetal'; end
            if nargin < 7 || isempty(methodV), methodV = 'centripetal'; end

            if size(Q,3) ~= 3
                error('NURBSSurface:globalLeastSquaresFitNet', 'Q must be nu x nv x 3.');
            end

            nv = size(Q,2);

            for j = 1:nv
                Cj = geom.NURBSCurve.globalLeastSquaresFit(squeeze(Q(:,j,:)), p, nCtrlU, methodU);
                if j == 1
                    Pu = zeros(size(Cj.P,1), nv, 3);
                    Wu = zeros(size(Cj.P,1), nv);
                    U = Cj.U;
                end
                Pu(:,j,:) = Cj.P;
                Wu(:,j) = Cj.W;
            end

            for i = 1:size(Pu,1)
                Ci = geom.NURBSCurve.globalLeastSquaresFit(squeeze(Pu(i,:,:)), q, nCtrlV, methodV);
                if i == 1
                    P = zeros(size(Pu,1), size(Ci.P,1), 3);
                    W = zeros(size(Pu,1), size(Ci.P,1));
                    V = Ci.U;
                end
                P(i,:,:) = Ci.P;
                W(i,:) = Ci.W.';
            end

            S = geom.NURBSSurface(P, p, q, U, V, W);
        end

        function S = ruled(C1, C2)
            if C1.p ~= C2.p || numel(C1.U) ~= numel(C2.U) || any(abs(C1.U - C2.U) > 1e-12)
                error('NURBSSurface:ruled', 'Curves must have matching degree and knot vector.');
            end
            if size(C1.P,1) ~= size(C2.P,1)
                error('NURBSSurface:ruled', 'Curves must have same number of control points.');
            end

            n = size(C1.P,1);
            P = zeros(n, 2, 3);
            W = zeros(n, 2);

            P(:,1,:) = C1.P;
            P(:,2,:) = C2.P;
            W(:,1) = C1.W;
            W(:,2) = C2.W;

            U = C1.U;
            V = [0 0 1 1];

            S = geom.NURBSSurface(P, C1.p, 1, U, V, W);
        end

        function S = bilinearCoons(P00, P10, P01, P11)
            P = zeros(2,2,3);
            P(1,1,:) = reshape(P00,1,1,3);
            P(2,1,:) = reshape(P10,1,1,3);
            P(1,2,:) = reshape(P01,1,1,3);
            P(2,2,:) = reshape(P11,1,1,3);

            W = ones(2,2);
            U = [0 0 1 1];
            V = [0 0 1 1];

            S = geom.NURBSSurface(P, 1, 1, U, V, W);
        end

        function [curvesOut, pmax, Uunion] = makeCompatibleCurves(curves)
            if ~iscell(curves), curves = {curves}; end
            if isempty(curves)
                error('NURBSSurface:makeCompatibleCurves', 'Need at least one curve.');
            end

            pmax = max(cellfun(@(c) c.p, curves));
            curvesOut = curves;

            for k = 1:numel(curvesOut)
                if curvesOut{k}.p < pmax
                    curvesOut{k} = curvesOut{k}.elevate(pmax - curvesOut{k}.p);
                end
            end

            Uunion = curvesOut{1}.U;
            for k = 2:numel(curvesOut)
                Uunion = union(Uunion, curvesOut{k}.U);
            end
            Uunion = sort(Uunion(:).');

            for val = Uunion
                multTarget = max(cellfun(@(c) sum(abs(c.U - val) < 1e-12), curvesOut));
                for k = 1:numel(curvesOut)
                    multK = sum(abs(curvesOut{k}.U - val) < 1e-12);
                    if multK < multTarget
                        curvesOut{k} = curvesOut{k}.insertKnot(val, multTarget - multK);
                    end
                end
            end

            Uunion = curvesOut{1}.U;
        end

        function [profilesC, guidesC, uPar, vPar] = makeCompatibleNetwork(profileCurves, guideCurves)
        % Practical compatibility helper for rectangular curve networks.
        %
        % Returns:
        %   profilesC : compatible profile curves
        %   guidesC   : compatible guide curves
        %   uPar      : normalized guide stations
        %   vPar      : normalized profile stations
            if ~iscell(profileCurves), profileCurves = {profileCurves}; end
            if ~iscell(guideCurves), guideCurves = {guideCurves}; end
            if isempty(profileCurves) || isempty(guideCurves)
                error('NURBSSurface:makeCompatibleNetwork', ...
                    'Need profile and guide curves.');
            end

            [profilesC, ~, ~] = geom.NURBSSurface.makeCompatibleCurves(profileCurves);
            [guidesC, ~, ~] = geom.NURBSSurface.makeCompatibleCurves(guideCurves);

            uPar = linspace(0,1,numel(guideCurves)).';
            vPar = linspace(0,1,numel(profileCurves)).';
        end

        function v = sectionParameters(curves, method)
            if nargin < 2 || isempty(method), method = 'centripetal'; end
            if ~iscell(curves), curves = {curves}; end

            nsec = numel(curves);
            Q = zeros(nsec, 3);
            for k = 1:nsec
                uk = 0.5 * (curves{k}.domain(1) + curves{k}.domain(2));
                Q(k,:) = curves{k}.evaluate(uk);
            end
            v = geom.NURBSCurve.parameterizeData(Q, method, min(3, nsec-1));
        end

        function S = loft(curves, q, method, sectionParams)
            if nargin < 2 || isempty(q), q = 3; end
            if nargin < 3 || isempty(method), method = 'centripetal'; end
            if ~iscell(curves), curves = {curves}; end
            if numel(curves) < 2
                error('NURBSSurface:loft', 'Need at least two section curves.');
            end

            [curvesC, pU, U] = geom.NURBSSurface.makeCompatibleCurves(curves);
            nsec = numel(curvesC);
            nCtrlU = size(curvesC{1}.P,1);

            if q > nsec-1
                error('NURBSSurface:loft', ...
                    'Loft degree q=%d exceeds number of sections-1 = %d.', q, nsec-1);
            end

            if nargin < 4 || isempty(sectionParams)
                vpar = geom.NURBSSurface.sectionParameters(curvesC, method);
            else
                vpar = sectionParams(:);
                if numel(vpar) ~= nsec
                    error('NURBSSurface:loft', 'sectionParams must match number of curves.');
                end
            end

            V = geom.NURBSCurve.averagingKnotVector(vpar, q);

            A = zeros(nsec, nsec);
            nV = nsec - 1;
            for r = 1:nsec
                span = geom.BasisFunctions.FindSpan(nV, q, vpar(r), V);
                N = geom.BasisFunctions.BasisFuns(span, vpar(r), q, V);
                cols = (span-q):span;
                A(r, cols) = N;
            end

            P = zeros(nCtrlU, nsec, 3);
            W = zeros(nCtrlU, nsec);

            for i = 1:nCtrlU
                H = zeros(nsec, 4);
                for s = 1:nsec
                    w = curvesC{s}.W(i);
                    H(s,:) = [w * curvesC{s}.P(i,:), w];
                end

                Qh = A \ H;
                if any(Qh(:,4) <= 0)
                    error('NURBSSurface:loft', ...
                        'Lofting produced nonpositive weights.');
                end

                W(i,:) = Qh(:,4).';
                P(i,:,:) = Qh(:,1:3) ./ Qh(:,4);
            end

            S = geom.NURBSSurface(P, pU, q, U, V, W);
        end


        function S = coons(Cu0, Cu1, Cv0, Cv1, p, q, nu, nv)
        %COONS Exact tensor-product Coons patch for polynomial B-spline boundaries.
        %
        % Boundary convention:
        %   Cu0(u) = S(u,0)
        %   Cu1(u) = S(u,1)
        %   Cv0(v) = S(0,v)
        %   Cv1(v) = S(1,v)
        %
        % This replaces the previous sampled/interpolated Coons implementation.
        %
        % Notes:
        %   - nu, nv are accepted for backward compatibility but unused.
        %   - This exact control-net form assumes non-rational / unit-weight
        %     boundary curves after compatibility. General rational Coons blending
        %     is not just P = Su + Sv - Sb in Cartesian control-net space.
        
            %#ok<INUSD> nu nv
        
            if nargin < 5 || isempty(p), p = max(Cu0.p, Cu1.p); end
            if nargin < 6 || isempty(q), q = max(Cv0.p, Cv1.p); end
        
            tol = 1e-10;
        
            % Corner consistency.
            P00a = Cu0.evaluate(Cu0.domain(1));
            P10a = Cu0.evaluate(Cu0.domain(2));
            P01a = Cu1.evaluate(Cu1.domain(1));
            P11a = Cu1.evaluate(Cu1.domain(2));
        
            P00b = Cv0.evaluate(Cv0.domain(1));
            P01b = Cv0.evaluate(Cv0.domain(2));
            P10b = Cv1.evaluate(Cv1.domain(1));
            P11b = Cv1.evaluate(Cv1.domain(2));
        
            if max([norm(P00a-P00b), norm(P10a-P10b), ...
                    norm(P01a-P01b), norm(P11a-P11b)]) > 1e-7
                error('NURBSSurface:coons', ...
                    'Boundary curve corners are inconsistent.');
            end
        
            % Make opposite boundary pairs compatible.
            [Cu, pU, U] = geom.NURBSSurface.makeCompatibleCurves({Cu0, Cu1});
            [Cv, qV, V] = geom.NURBSSurface.makeCompatibleCurves({Cv0, Cv1});
        
            Cu0c = Cu{1};
            Cu1c = Cu{2};
            Cv0c = Cv{1};
            Cv1c = Cv{2};
        
            % Optional requested degree elevation.
            if p > pU
                Cu0c = Cu0c.elevate(p - pU);
                Cu1c = Cu1c.elevate(p - pU);
                [Cu, pU, U] = geom.NURBSSurface.makeCompatibleCurves({Cu0c, Cu1c});
                Cu0c = Cu{1};
                Cu1c = Cu{2};
            else
                p = pU;
            end
        
            if q > qV
                Cv0c = Cv0c.elevate(q - qV);
                Cv1c = Cv1c.elevate(q - qV);
                [Cv, qV, V] = geom.NURBSSurface.makeCompatibleCurves({Cv0c, Cv1c});
                Cv0c = Cv{1};
                Cv1c = Cv{2};
            else
                q = qV;
            end
        
            % This simple exact Coons construction is polynomial.
            if any(abs(Cu0c.W - 1) > tol) || any(abs(Cu1c.W - 1) > tol) || ...
               any(abs(Cv0c.W - 1) > tol) || any(abs(Cv1c.W - 1) > tol)
                error('NURBSSurface:coons', ...
                    ['Exact control-net Coons currently supports unit-weight ', ...
                     'boundary curves only. Rational boundaries need a separate ', ...
                     'product-space rational Coons implementation.']);
            end
        
            nuCtrl = size(Cu0c.P,1);
            nvCtrl = size(Cv0c.P,1);
        
            % Corner points.
            P00 = P00a;
            P10 = P10a;
            P01 = P01a;
            P11 = P11a;
        
            % Tensor-product Coons control net:
            %
            % S(u,v) = (1-v) Cu0(u) + v Cu1(u)
            %        + (1-u) Cv0(v) + u Cv1(v)
            %        - bilinear_corners(u,v)
            %
            Pnet = zeros(nuCtrl, nvCtrl, 3);
        
            for i = 1:nuCtrl
                uhat = grevilleLocal(i, p, U);
        
                for j = 1:nvCtrl
                    vhat = grevilleLocal(j, q, V);
        
                    Pu0 = reshape(Cu0c.P(i,:), 1, 3);
                    Pu1 = reshape(Cu1c.P(i,:), 1, 3);
                    Pv0 = reshape(Cv0c.P(j,:), 1, 3);
                    Pv1 = reshape(Cv1c.P(j,:), 1, 3);
        
                    Su = (1 - vhat) * Pu0 + vhat * Pu1;
                    Sv = (1 - uhat) * Pv0 + uhat * Pv1;
        
                    Sb = (1 - uhat) * (1 - vhat) * P00 + ...
                          uhat      * (1 - vhat) * P10 + ...
                         (1 - uhat) * vhat       * P01 + ...
                          uhat      * vhat       * P11;
        
                    Pnet(i,j,:) = reshape(Su + Sv - Sb, 1, 1, 3);
                end
            end
        
            W = ones(nuCtrl, nvCtrl);
            S = geom.NURBSSurface(Pnet, p, q, U, V, W);
        
            function g = grevilleLocal(ii, deg, K)
                % ii is MATLAB 1-based control-point index.
                % For clamped Bezier this gives the Bernstein control-position
                % parameter, and for general B-splines gives the Greville abscissa.
                if deg == 0
                    g = 0.0;
                else
                    g = sum(K(ii+1:ii+deg)) / deg;
                end
            end
        end

        function S = gordon(profileCurves, guideCurves, p, q, nu, nv)
        %GORDON Exact Gordon surface, Algorithm A10.3-style.
        %
        % profileCurves: curves running in u-direction, placed at v stations
        % guideCurves:   curves running in v-direction, placed at u stations
        %
        % S = L1 + L2 - T
        %
        % nu,nv are accepted for backward compatibility but unused.
        
            %#ok<INUSD> nu nv
        
            if nargin < 3 || isempty(p), p = 3; end
            if nargin < 4 || isempty(q), q = 3; end
        
            if ~iscell(profileCurves), profileCurves = {profileCurves}; end
            if ~iscell(guideCurves),   guideCurves   = {guideCurves};   end
        
            nProf = numel(profileCurves);
            nGuid = numel(guideCurves);
        
            if nProf < 2 || nGuid < 2
                error('NURBSSurface:gordon', ...
                    'Need at least two profile curves and two guide curves.');
            end
        
            if p > nGuid - 1
                error('NURBSSurface:gordon', ...
                    'Requested p=%d exceeds number of guide stations-1=%d.', p, nGuid-1);
            end
            if q > nProf - 1
                error('NURBSSurface:gordon', ...
                    'Requested q=%d exceeds number of profile stations-1=%d.', q, nProf-1);
            end
        
            tol = 1e-9;
        
            % Uniform Gordon network parameters.
            uPar = linspace(0, 1, nGuid).';
            vPar = linspace(0, 1, nProf).';
        
            % L1(u,v): loft profile curves in v direction.
            L1 = geom.NURBSSurface.loft(profileCurves, q, 'uniform', vPar);
        
            % L2(u,v): loft guide curves in u direction, then swap.
            L2 = geom.NURBSSurface.loft(guideCurves, p, 'uniform', uPar).swapUV();
        
            % Intersection net T(u,v).
            X = zeros(nGuid, nProf, 3);
        
            for i = 1:nGuid
                for j = 1:nProf
                    Pi = profileCurves{j}.evaluate( ...
                        map01ToDomain(uPar(i), profileCurves{j}.domain));
        
                    Pj = guideCurves{i}.evaluate( ...
                        map01ToDomain(vPar(j), guideCurves{i}.domain));
        
                    if norm(Pi - Pj) > tol
                        error('NURBSSurface:gordon', ...
                            'Curve network mismatch at guide %d, profile %d: err = %.3e.', ...
                            i, j, norm(Pi - Pj));
                    end
        
                    X(i,j,:) = reshape(0.5*(Pi + Pj), 1, 1, 3);
                end
            end
        
            T = interpSurfaceNetWithParams(X, p, q, uPar, vPar);
        
            % Final common degree.
            pFinal = max([L1.p, L2.p, T.p]);
            qFinal = max([L1.q, L2.q, T.q]);
        
            L1 = elevateTo(L1, pFinal, qFinal);
            L2 = elevateTo(L2, pFinal, qFinal);
            T  = elevateTo(T,  pFinal, qFinal);
        
            % Merge knot vectors, preserving maximum multiplicities.
            U = mergeKnotsMaxMult(L1.U, L2.U, T.U, 1e-12);
            V = mergeKnotsMaxMult(L1.V, L2.V, T.V, 1e-12);
        
            L1 = refineToKnots(L1, U, V, 1e-12);
            L2 = refineToKnots(L2, U, V, 1e-12);
            T  = refineToKnots(T,  U, V, 1e-12);
        
            % This exact control-net addition assumes polynomial/unit weights.
            if max(abs(L1.W(:)-1)) > tol || max(abs(L2.W(:)-1)) > tol || max(abs(T.W(:)-1)) > tol
                error('NURBSSurface:gordon', ...
                    'Exact Gordon control-net blend currently supports unit-weight surfaces only.');
            end
        
            P = L1.P + L2.P - T.P;
            W = ones(size(P,1), size(P,2));
        
            S = geom.NURBSSurface(P, pFinal, qFinal, U, V, W);
        
            function t = map01ToDomain(a, dom)
                t = dom(1) + a*(dom(2)-dom(1));
            end
        
            function Sout = elevateTo(Sin, pp, qq)
                Sout = Sin;
                if Sout.p < pp
                    Sout = Sout.elevateU(pp - Sout.p);
                end
                if Sout.q < qq
                    Sout = Sout.elevateV(qq - Sout.q);
                end
            end
        
            function Sfit = interpSurfaceNetWithParams(Q, pp, qq, up, vp)
                nuQ = size(Q,1);
                nvQ = size(Q,2);
        
                Ufit = geom.NURBSCurve.averagingKnotVector(up, pp);
                Vfit = geom.NURBSCurve.averagingKnotVector(vp, qq);
        
                Au = zeros(nuQ, nuQ);
                for a = 1:nuQ
                    span = geom.BasisFunctions.FindSpan(nuQ-1, pp, up(a), Ufit);
                    N = geom.BasisFunctions.BasisFuns(span, up(a), pp, Ufit);
                    cols = (span-pp):span;
                    Au(a, cols) = N;
                end
        
                Av = zeros(nvQ, nvQ);
                for b = 1:nvQ
                    span = geom.BasisFunctions.FindSpan(nvQ-1, qq, vp(b), Vfit);
                    N = geom.BasisFunctions.BasisFuns(span, vp(b), qq, Vfit);
                    cols = (span-qq):span;
                    Av(b, cols) = N;
                end
        
                Pfit = zeros(nuQ, nvQ, 3);
        
                for d = 1:3
                    Qd = Q(:,:,d);
                    Pfit(:,:,d) = Au \ Qd / Av.';
                end
        
                Sfit = geom.NURBSSurface(Pfit, pp, qq, Ufit, Vfit, ones(nuQ,nvQ));
            end
        
            function K = mergeKnotsMaxMult(varargin)
                tolK = varargin{end};
                knots = varargin(1:end-1);
        
                vals = [];
                for kk = 1:numel(knots)
                    vals = [vals, knots{kk}(:).']; %#ok<AGROW>
                end
                vals = sort(vals);
        
                uniqueVals = [];
                while ~isempty(vals)
                    v0 = vals(1);
                    uniqueVals(end+1) = v0; %#ok<AGROW>
                    vals(abs(vals - v0) < tolK) = [];
                end
        
                K = [];
                for a = 1:numel(uniqueVals)
                    v0 = uniqueVals(a);
                    maxMult = 0;
                    for kk = 1:numel(knots)
                        maxMult = max(maxMult, sum(abs(knots{kk} - v0) < tolK));
                    end
                    K = [K, repmat(v0, 1, maxMult)]; %#ok<AGROW>
                end
            end
        
            function Sout = refineToKnots(Sin, Ut, Vt, tolK)
                Sout = Sin;
        
                Xu = [];
                uInt = unique(Ut);
                for a = 1:numel(uInt)
                    val = uInt(a);
                    if val <= Sout.domainU(1)+tolK || val >= Sout.domainU(2)-tolK
                        continue;
                    end
                    mt = sum(abs(Ut - val) < tolK);
                    ms = sum(abs(Sout.U - val) < tolK);
                    if mt > ms
                        Xu = [Xu, repmat(val, 1, mt-ms)]; %#ok<AGROW>
                    end
                end
        
                Xv = [];
                vInt = unique(Vt);
                for a = 1:numel(vInt)
                    val = vInt(a);
                    if val <= Sout.domainV(1)+tolK || val >= Sout.domainV(2)-tolK
                        continue;
                    end
                    mt = sum(abs(Vt - val) < tolK);
                    ms = sum(abs(Sout.V - val) < tolK);
                    if mt > ms
                        Xv = [Xv, repmat(val, 1, mt-ms)]; %#ok<AGROW>
                    end
                end
        
                Sout = Sout.refine(Xu, Xv);
            end
        end


        function S = revolve(C, axisPoint, axisDir, thetaTotal, q, nSections)
            %REVOLVE Exact rational NURBS surface of revolution.
            %
            % Drop-in replacement for the old approximate loft-based revolve.
            %
            % Signature keeps q and nSections for backward compatibility, but they
            % are intentionally unused. The revolution direction is exactly quadratic
            % rational, per Piegl & Tiller Algorithm A8.1.
            %
            % Surface directions:
            %   u = original curve direction
            %   v = revolution direction
        
            %#ok<INUSD> q nSections
        
            if nargin < 4 || isempty(thetaTotal)
                thetaTotal = 2*pi;
            end
        
            axisPoint = axisPoint(:).';
            axisDir   = axisDir(:).';
        
            if norm(axisDir) < eps
                error('NURBSSurface:revolve', 'axisDir must be nonzero.');
            end
        
            T = axisDir / norm(axisDir);
        
            if abs(thetaTotal) < eps
                error('NURBSSurface:revolve', 'thetaTotal must be nonzero.');
            end
        
            if abs(thetaTotal) > 2*pi + 1e-12
                error('NURBSSurface:revolve', ...
                    'Exact revolve supports |thetaTotal| <= 2*pi.');
            end
        
            sgn = sign(thetaTotal);
            theta = abs(thetaTotal);
        
            % Algorithm A8.1 arc count.
            if theta <= pi/2 + 1e-12
                narcs = 1;
            elseif theta <= pi + 1e-12
                narcs = 2;
            elseif theta <= 3*pi/2 + 1e-12
                narcs = 3;
            else
                narcs = 4;
            end
        
            dtheta = sgn * theta / narcs;
            wm = cos(dtheta / 2);
        
            % Quadratic clamped knot vector in revolution direction.
            switch narcs
                case 1
                    V = [0 0 0 1 1 1];
                case 2
                    V = [0 0 0 0.5 0.5 1 1 1];
                case 3
                    V = [0 0 0 1/3 1/3 2/3 2/3 1 1 1];
                case 4
                    V = [0 0 0 0.25 0.25 0.5 0.5 0.75 0.75 1 1 1];
            end
        
            qrev = 2;
            nCurve = size(C.P, 1);
            nRev   = 2 * narcs + 1;
        
            P = zeros(nCurve, nRev, 3);
            W = zeros(nCurve, nRev);
        
            Pc = C.P;
            Wc = C.W(:);
        
            for j = 1:nCurve
                Pj = Pc(j,:);
        
                % Center of rotation for this control point.
                O = axisPoint + dot(Pj - axisPoint, T) * T;
        
                X = Pj - O;
                r = norm(X);
        
                % Point lies on axis: revolution degenerates to the same point.
                if r < 1e-14
                    for k = 1:nRev
                        P(j,k,:) = reshape(Pj, 1, 1, 3);
                        W(j,k) = Wc(j);
                    end
                    continue;
                end
        
                X = X / r;
                Y = cross(T, X);
                Y = Y / norm(Y);
        
                angle0 = 0.0;
        
                P0 = Pj;
                P(j,1,:) = reshape(P0, 1, 1, 3);
                W(j,1) = Wc(j);
        
                idx = 1;
        
                for iarc = 1:narcs
                    angle1 = angle0 + dtheta;
                    angleM = angle0 + 0.5*dtheta;
        
                    P2 = O + r*cos(angle1)*X + r*sin(angle1)*Y;
        
                    % Exact rational conic middle control point.
                    % Equivalent to intersecting endpoint tangent lines.
                    P1 = O + (r / wm) * (cos(angleM)*X + sin(angleM)*Y);
        
                    P(j,idx+1,:) = reshape(P1, 1, 1, 3);
                    W(j,idx+1) = Wc(j) * wm;
        
                    P(j,idx+2,:) = reshape(P2, 1, 1, 3);
                    W(j,idx+2) = Wc(j);
        
                    idx = idx + 2;
                    angle0 = angle1;
                end
            end
        
            S = geom.NURBSSurface(P, C.p, qrev, C.U, V, W);
        end

        function S = sweep(profileCurve, spineCurve, q, nStations, upVec)
            if nargin < 3 || isempty(q), q = 3; end
            if nargin < 4 || isempty(nStations), nStations = 9; end
            if nargin < 5 || isempty(upVec), upVec = [0 0 1]; end

            upVec = upVec(:).';
            upVec = upVec / norm(upVec);

            ts = spineCurve.arcLengthParam(nStations);
            sections = cell(1, nStations);

            for k = 1:nStations
                tk = ts(k);
                P0 = spineCurve.evaluate(tk);
                T0 = spineCurve.tangent(tk);
                T0 = T0 / norm(T0);

                B0 = cross(T0, upVec);
                if norm(B0) < 1e-8
                    alt = [1 0 0];
                    if abs(dot(T0, alt)) > 0.9
                        alt = [0 1 0];
                    end
                    B0 = cross(T0, alt);
                end
                B0 = B0 / norm(B0);
                N0 = cross(B0, T0);
                N0 = N0 / norm(N0);

                R = [T0(:), N0(:), B0(:)];
                H = eye(4);
                H(1:3,1:3) = R;
                H(1:3,4) = P0(:);

                sections{k} = profileCurve.transform(H);
            end

            S = geom.NURBSSurface.loft(sections, q, 'chord');
        end
    end

    methods (Static, Access = private)

        function [Qw, Ubar] = curveKnotRefineHomogeneous(p, U, Pw, X)
        %CURVEKNOTREFINEHOMOGENEOUS Batch curve knot refinement.
        % Implements Piegl & Tiller Algorithm A5.4, adapted to MATLAB indexing.
        %
        % Inputs:
        %   p  : degree
        %   U  : original knot vector
        %   Pw : [n+1 x dim] homogeneous control points
        %   X  : refinement knot vector, sorted or unsorted
        %
        % Outputs:
        %   Qw   : refined homogeneous control points
        %   Ubar : refined knot vector
        
            X = sort(X(:).');
            U = U(:).';
        
            if isempty(X)
                Qw = Pw;
                Ubar = U;
                return;
            end
        
            n = size(Pw,1) - 1;        % zero-based last control index
            m = n + p + 1;             % zero-based last knot index
            r = numel(X) - 1;          % zero-based last refinement-knot index
            dim = size(Pw,2);
        
            a = geom.BasisFunctions.FindSpan(n, p, X(1), U) - 1;
            b = geom.BasisFunctions.FindSpan(n, p, X(end), U) - 1;
            b = b + 1;
        
            Qw   = zeros(n + r + 2, dim);
            Ubar = zeros(1, m + r + 2);
        
            % Save unaltered control points.
            for j = 0:(a-p)
                Qw(j+1,:) = Pw(j+1,:);
            end
        
            for j = (b-1):n
                Qw(j+r+2,:) = Pw(j+1,:);
            end
        
            % Save unaltered knots.
            for j = 0:a
                Ubar(j+1) = U(j+1);
            end
        
            for j = (b+p):m
                Ubar(j+r+2) = U(j+1);
            end
        
            i = b + p - 1;
            k = b + p + r;
        
            for j = r:-1:0
                xj = X(j+1);
        
                while xj <= U(i+1) && i > a
                    Qw(k-p,:) = Pw(i-p,:);
                    Ubar(k+1) = U(i+1);
                    k = k - 1;
                    i = i - 1;
                end
        
                Qw(k-p,:) = Qw(k-p+1,:);
        
                for l = 1:p
                    ind = k - p + l;
        
                    denom = Ubar(k+l+1) - U(i-p+l+1);
        
                    if abs(denom) < eps
                        alpha = 0.0;
                    else
                        alpha = (Ubar(k+l+1) - xj) / denom;
                    end
        
                    if abs(alpha) < eps
                        Qw(ind,:) = Qw(ind+1,:);
                    else
                        Qw(ind,:) = alpha * Qw(ind,:) + ...
                                    (1.0 - alpha) * Qw(ind+1,:);
                    end
                end
        
                Ubar(k+1) = xj;
                k = k - 1;
            end
        end

        function T = axisRotationTransform(point, axisDir, theta)
            axisDir = axisDir(:) / norm(axisDir);
            x = axisDir(1); y = axisDir(2); z = axisDir(3);
            c = cos(theta);
            s = sin(theta);
            C = 1 - c;

            R = [x*x*C + c,   x*y*C - z*s, x*z*C + y*s; ...
                 y*x*C + z*s, y*y*C + c,   y*z*C - x*s; ...
                 z*x*C - y*s, z*y*C + x*s, z*z*C + c];

            p = point(:);
            t = p - R*p;

            T = eye(4);
            T(1:3,1:3) = R;
            T(1:3,4) = t;
        end

        function N = allBasisFuns(i, u, p, U)
            N = zeros(p+1,p+1);
            N(1,1) = 1.0;
            left = zeros(1,p+1);
            right = zeros(1,p+1);
        
            for j = 1:p
            left(j+1)  = u - U(i+1-j);
            right(j+1) = U(i+j) - u;
                saved = 0.0;
        
                for r = 0:j-1
                    temp = N(j,r+1) / (right(r+2) + left(j-r+1));
                    N(j+1,r+1) = saved + right(r+2) * temp;
                    saved = left(j-r+1) * temp;
                end
        
                N(j+1,j+1) = saved;
            end
        end


        function PK = curveDerivCpts1D(p, U, P, d, r1, r2)
            % P is local block P(r1:r2), 1-based input indices.
            r = r2 - r1;
            du = min(d,p);
        
            PK = cell(du+1,1);
            PK{1} = P;
        
            for k = 1:du
                prev = PK{k};
                cur = zeros(r-k+1, size(P,2));
        
                for i = 0:r-k
                    global0 = (r1-1) + i; % book-style zero-based control index
                    denom = U(global0+p+2) - U(global0+k+1);
                    cur(i+1,:) = (p-k+1) * (prev(i+2,:) - prev(i+1,:)) / denom;
                end
        
                PK{k+1} = cur;
            end
        end

        function SKL = rationalSurfaceDerivativesFromHomogeneous(Aders,wders,d)
            SKL = cell(d+1,d+1);
            for k = 0:d
                for l = 0:d
                    SKL{k+1,l+1} = zeros(1,3);
                end
            end
        
            w00 = wders(1,1);
            if abs(w00) < eps
                error('NURBSSurface:ZeroWeight', ...
                    'Surface weight function vanished.');
            end
        
            for k = 0:d
                for l = 0:(d-k)
                    vkl = Aders{k+1,l+1};
        
                    for j = 1:l
                        vkl = vkl - nchoosek(l,j) * wders(1,j+1) * SKL{k+1,l-j+1};
                    end
        
                    for i = 1:k
                        vkl = vkl - nchoosek(k,i) * wders(i+1,1) * SKL{k-i+1,l+1};
                        for j = 1:l
                            vkl = vkl - nchoosek(k,i) * nchoosek(l,j) * ...
                                wders(i+1,j+1) * SKL{k-i+1,l-j+1};
                        end
                    end
        
                    SKL{k+1,l+1} = vkl / w00;
                end
            end
        end

        function [Qw, Uq] = curveKnotInsertHomogeneous(p, U, Pw, u, r)
        %CURVEKNOTINSERTHOMOGENEOUS Exact homogeneous curve knot insertion.
        % Implements Piegl & Tiller Algorithm A5.1, CurveKnotIns.
        %
        % Pw is [n+1 x dim], usually dim=4.
        
            if nargin < 5 || isempty(r), r = 1; end
            r = floor(r);
        
            if r < 0
                error('NURBSSurface:curveKnotInsertHomogeneous', ...
                    'r must be a nonnegative integer.');
            end
        
            U = U(:).';
        
            if r == 0
                Qw = Pw;
                Uq = U;
                return;
            end
        
            n  = size(Pw,1) - 1;       % book zero-based n
            mp = n + p + 1;            % book m = n+p+1
        
            % Your FindSpan returns MATLAB 1-based span.
            k1 = geom.BasisFunctions.FindSpan(n, p, u, U);
            k  = k1 - 1;               % convert to book zero-based span
        
            s = sum(abs(U - u) < 1e-12);
        
            if s + r > p
                error('NURBSSurface:curveKnotInsertHomogeneous', ...
                    'Cannot insert knot: multiplicity s+r exceeds degree p.');
            end
        
            nq = n + r;
        
            Uq = zeros(1, mp + r + 1);
            Qw = zeros(nq + 1, size(Pw,2));
        
            % ---- Load new knot vector ----
            for i = 0:k
                Uq(i+1) = U(i+1);
            end
        
            for i = 1:r
                Uq(k+i+1) = u;
            end
        
            for i = k+1:mp
                Uq(i+r+1) = U(i+1);
            end
        
            % ---- Save unaltered control points ----
            for i = 0:k-p
                Qw(i+1,:) = Pw(i+1,:);
            end
        
            for i = k-s:n
                Qw(i+r+1,:) = Pw(i+1,:);
            end
        
            % ---- Local affected control points ----
            Rw = zeros(p-s+1, size(Pw,2));
            for i = 0:p-s
                Rw(i+1,:) = Pw(k-p+i+1,:);
            end
        
            % ---- Insert the knot r times ----
            L = 0;
            for j = 1:r
                L = k - p + j;
        
                for i = 0:p-j-s
                    alpha = (u - U(L+i+1)) / (U(i+k+2) - U(L+i+1));
                    Rw(i+1,:) = alpha * Rw(i+2,:) + (1.0 - alpha) * Rw(i+1,:);
                end
        
                Qw(L+1,:) = Rw(1,:);
                Qw(k+r-j-s+1,:) = Rw(p-j-s+1,:);
            end
        
            % ---- Load remaining control points ----
            for i = L+1:k-s-1
                Qw(i+1,:) = Rw(i-L+1,:);
            end
        end





    end





end