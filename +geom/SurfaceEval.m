classdef SurfaceEval
% SURFACEEVAL  Analytic NURBS surface derivative utilities.
%
% Drop-in utility for projects that already have geom.NURBSSurface and
% geom.BasisFunctions. This avoids finite differencing for Suu, Suv, Svv,
% curvature, and Newton projection/intersection work.
%
% Primary entry points:
%   SKL  = geom.SurfaceEval.derivatives(S, u, v, d)
%   [Su,Sv,Suu,Suv,Svv] = geom.SurfaceEval.firstSecond(S, u, v)
%   [K,H,N] = geom.SurfaceEval.curvatures(S, u, v)
%
% SKL is returned as a cell array where SKL{k+1,l+1} is the (k,l)
% derivative d^(k+l)S / du^k dv^l.
%
% References:
%   Piegl & Tiller, The NURBS Book, 2nd ed., Algorithms A3.6 / A4.4.

    methods (Static)
        function SKL = derivatives(S, u, v, d)
        % DERIVATIVES  Analytic rational surface derivatives up to order d.
        %
        %   SKL = geom.SurfaceEval.derivatives(S, u, v, d)
        %
        % Inputs
        %   S  - geom.NURBSSurface
        %   u  - scalar parameter in u direction
        %   v  - scalar parameter in v direction
        %   d  - maximum total derivative order in each direction
        %
        % Output
        %   SKL - cell(d+1,d+1); SKL{k+1,l+1} is [1x3]
        %         derivative d^(k+l)S/du^k dv^l

            if nargin < 4 || isempty(d)
                d = 1;
            end

            u = geom.SurfaceEval.clampToDomain(u, S.domainU);
            v = geom.SurfaceEval.clampToDomain(v, S.domainV);

            du = min(d, S.p);
            dv = min(d, S.q);

            uspan = geom.BasisFunctions.FindSpan(S.n, S.p, u, S.U);
            vspan = geom.BasisFunctions.FindSpan(S.m, S.q, v, S.V);
            Nu = geom.BasisFunctions.DersBasisFuns(uspan, u, S.p, du, S.U);
            Nv = geom.BasisFunctions.DersBasisFuns(vspan, v, S.q, dv, S.V);

            % Homogeneous derivatives: Aders{k+1,l+1} and wders(k+1,l+1)
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
                    for i = 0:S.p
                        ii = uspan - S.p + i;
                        for j = 0:S.q
                            jj = vspan - S.q + j;
                            B = Nu(k+1, i+1) * Nv(l+1, j+1);
                            wij = S.W(ii, jj);
                            Pij = squeeze(S.P(ii, jj, :))';
                            Akl = Akl + B * wij * Pij;
                            wkl = wkl + B * wij;
                        end
                    end
                    Aders{k+1,l+1} = Akl;
                    wders(k+1,l+1) = wkl;
                end
            end

            SKL = cell(d+1, d+1);
            for k = 0:d
                for l = 0:d
                    SKL{k+1,l+1} = zeros(1,3);
                end
            end

            w00 = wders(1,1);
            if abs(w00) < eps
                error('geom.SurfaceEval:ZeroWeight', ...
                    'Surface weight function vanished at (u,v).');
            end

            for k = 0:d
                for l = 0:d
                    vkl = Aders{k+1,l+1};

                    % subtract v-direction rational terms
                    for j = 1:l
                        vkl = vkl - nchoosek(l,j) * wders(1,j+1) * SKL{k+1,l-j+1};
                    end

                    % subtract u-direction and mixed terms
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

        function [Su, Sv, Suu, Suv, Svv] = firstSecond(S, u, v)
        % FIRSTSECOND  Convenience wrapper for first and second derivatives.
            SKL = geom.SurfaceEval.derivatives(S, u, v, 2);
            Su  = SKL{2,1};
            Sv  = SKL{1,2};
            Suu = SKL{3,1};
            Suv = SKL{2,2};
            Svv = SKL{1,3};
        end

        function [K, H, Nhat] = curvatures(S, u, v)
        % CURVATURES  Gaussian and mean curvature from analytic derivatives.
            SKL = geom.SurfaceEval.derivatives(S, u, v, 2);
            Su  = SKL{2,1};
            Sv  = SKL{1,2};
            Suu = SKL{3,1};
            Suv = SKL{2,2};
            Svv = SKL{1,3};

            N = cross(Su, Sv);
            nrm = norm(N);
            if nrm < eps
                Nhat = [0 0 0];
                K = NaN; H = NaN;
                return;
            end
            Nhat = N / nrm;

            E = dot(Su, Su);
            F = dot(Su, Sv);
            G = dot(Sv, Sv);
            L = dot(Suu, Nhat);
            M = dot(Suv, Nhat);
            N2 = dot(Svv, Nhat);
            den = E*G - F*F;
            if abs(den) < eps
                K = NaN; H = NaN;
                return;
            end
            K = (L*N2 - M*M) / den;
            H = (E*N2 - 2*F*M + G*L) / (2*den);
        end

        function S = point(Surf, u, v)
        % POINT  Evaluate point using derivative utility storage path.
            SKL = geom.SurfaceEval.derivatives(Surf, u, v, 0);
            S = SKL{1,1};
        end
    end

    methods (Static, Access = private)
        function u = clampToDomain(u, dom)
            u = max(dom(1), min(dom(2), u));
        end
    end
end
