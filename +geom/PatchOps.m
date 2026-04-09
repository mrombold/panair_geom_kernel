classdef PatchOps
% PATCHOPS  Exact patch split/extract helpers for tensor-product NURBS surfaces.
%
% Includes:
%   [Slo, Shi] = splitU(S, u0)
%   [Slo, Shi] = splitV(S, v0)
%   Ssub       = extractUV(S, urange, vrange)
%   B          = boundaryCurves(S)
%   C          = extractInnerSpanCurve(S)
%   C          = extractInnerCircCurve(S)

    methods (Static)

        function [Slo, Shi] = splitU(S, u0)
            tol = 1e-12;
            if u0 <= S.U(1) + tol || u0 >= S.U(end) - tol
                error('geom.PatchOps:splitU', ...
                      'u0 must lie strictly inside the U domain.');
            end

            p = S.p;
            U = S.U;

            s = sum(abs(U - u0) < tol);
            if s > p
                error('geom.PatchOps:splitU', ...
                      'Split knot multiplicity exceeds degree.');
            end

            reps = p - s;
            if reps > 0
                Sref = S.refine(repmat(u0, 1, reps), []);
            else
                Sref = S;
            end

            U2 = Sref.U;
            first = find(abs(U2 - u0) < tol, 1, 'first');
            last  = find(abs(U2 - u0) < tol, 1, 'last');
            if isempty(first) || isempty(last)
                error('geom.PatchOps:splitU', ...
                      'Failed to create split knot in refined surface.');
            end

            cp_split = last - p;

            P2 = Sref.P;
            W2 = Sref.W;

            Plo = P2(1:cp_split, :, :);
            Wlo = W2(1:cp_split, :);
            Ulo = [U2(1:last), u0];

            Phi = P2(cp_split:end, :, :);
            Whi = W2(cp_split:end, :);
            Uhi = [u0, U2(first:end)];

            Ulo = (Ulo - Ulo(1)) / (Ulo(end) - Ulo(1));
            Uhi = (Uhi - Uhi(1)) / (Uhi(end) - Uhi(1));

            Slo = geom.NURBSSurface(Plo, S.p, S.q, Ulo, S.V, Wlo);
            Shi = geom.NURBSSurface(Phi, S.p, S.q, Uhi, S.V, Whi);
        end

        function [Slo, Shi] = splitV(S, v0)
            tol = 1e-12;
            if v0 <= S.V(1) + tol || v0 >= S.V(end) - tol
                error('geom.PatchOps:splitV', ...
                      'v0 must lie strictly inside the V domain.');
            end

            St = geom.PatchOps.swapUV(S);
            [A, B] = geom.PatchOps.splitU(St, v0);
            Slo = geom.PatchOps.swapUV(A);
            Shi = geom.PatchOps.swapUV(B);
        end

        function Ssub = extractUV(S, urange, vrange)
            ua = urange(1); ub = urange(2);
            va = vrange(1); vb = vrange(2);

            if ~(ua < ub && va < vb)
                error('geom.PatchOps:extractUV', ...
                      'Require ua < ub and va < vb.');
            end

            tol = 1e-12;

            if ua > S.U(1) + tol
                [~, S] = geom.PatchOps.splitU(S, ua);
            end

            if ub < 1 - tol
                ub_local = (ub - ua) / (1 - ua);
                [S, ~] = geom.PatchOps.splitU(S, ub_local);
            end

            if va > S.V(1) + tol
                [~, S] = geom.PatchOps.splitV(S, va);
            end

            if vb < 1 - tol
                vb_local = (vb - va) / (1 - va);
                [S, ~] = geom.PatchOps.splitV(S, vb_local);
            end

            Ssub = S;
        end

        function B = boundaryCurves(S)
            B.u0 = S.isoCurveV(S.U(1));
            B.u1 = S.isoCurveV(S.U(end));
            B.v0 = S.isoCurveU(S.V(1));
            B.v1 = S.isoCurveU(S.V(end));
        end

        function C = extractInnerSpanCurve(S)
        % EXTRACTINNERSPANCURVE  Return a representative spanwise/inboard-side curve.
        %
        % Intended for wing-like patches. We pick a curve slightly inside the
        % lower-U boundary to avoid edge degeneracy after repeated splitting.

            u0 = S.domainU(1);
            u1 = S.domainU(2);
            u_pick = u0 + 0.05 * (u1 - u0);
            u_pick = max(u0, min(u1, u_pick));
            C = S.isoCurveV(u_pick);
        end

        function C = extractInnerCircCurve(S)
        % EXTRACTINNERCIRCCURVE  Return a representative circumferential-side curve.
        %
        % Intended for fuselage-like mid patches. We pick a curve slightly
        % inside the lower-V boundary to avoid edge degeneracy after repeated splitting.

            v0 = S.domainV(1);
            v1 = S.domainV(2);
            v_pick = v0 + 0.05 * (v1 - v0);
            v_pick = max(v0, min(v1, v_pick));
            C = S.isoCurveU(v_pick);
        end

        function St = swapUV(S)
            P = permute(S.P, [2 1 3]);
            W = S.W.';
            St = geom.NURBSSurface(P, S.q, S.p, S.V, S.U, W);
        end
    end
end