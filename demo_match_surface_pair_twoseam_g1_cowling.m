function demo_match_surface_pair_twoseam_g1_cowling()
clear; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir), addpath(script_dir); end
if ~exist('+geom', 'dir') && isempty(which('geom.NURBSCurve'))
    error(['Cannot find +geom package.\n' ...
           'Run this from the repo root or add the geometry kernel to path.']);
end
fprintf('geom package found at: %s\n', script_dir);

[Upr,~] = geom.Loft.limingConic([4.6328 0 2.9060], [18 0 6.31], ...
    [4.6328+(6.31-2.9060)*tan(24*pi/180) 0 6.31], [9 0 5.35]);
[MaxBp,~] = geom.Loft.limingConic([5.7352 2.941 0], [18 6.31 0], ...
    [5.7352+(6.31-2.941)*tan(30*pi/180) 6.31 0], [9 5.14 0]);
[MaxBn,~] = geom.Loft.limingConic([5.7352 -2.941 0], [18 -6.31 0], ...
    [5.7352+(6.31-2.941)*tan(30*pi/180) -6.31 0], [9 -5.14 0]);
[Lwr,~] = geom.Loft.limingConic([6.1476 0 -2.9739], [18 0 -6.31], ...
    [6.1476+(6.31-2.9739)*tan(36*pi/180) 0 -6.31], [9 0 -4.93]);
[UprShldrP,~] = geom.Loft.limingConic([4.8509 2.9237*sin(pi/4) 2.9237*cos(pi/4)], ...
    [18 6.31*sin(pi/4) 6.31*cos(pi/4)], ...
    [4.6328+(6.31-2.9237)*tan(27*pi/180) 6.31*sin(pi/4) 6.31*cos(pi/4)], ...
    [9 5.35*sin(pi/4) 5.35*cos(pi/4)]);
[UprShldrN,~] = geom.Loft.limingConic([4.8509 2.9237*sin(-pi/4) 2.9237*cos(-pi/4)], ...
    [18 6.31*sin(-pi/4) 6.31*cos(-pi/4)], ...
    [4.6328+(6.31-2.9237)*tan(27*pi/180) 6.31*sin(-pi/4) 6.31*cos(-pi/4)], ...
    [9 5.35*sin(-pi/4) 5.35*cos(-pi/4)]);
[LwrShldrP,~] = geom.Loft.limingConic([5.9139 2.9739*sin(pi/4) -2.9739*cos(pi/4)], ...
    [18 6.31*sin(pi/4) -6.31*cos(pi/4)], ...
    [4.6328+(6.31-2.9739)*tan(27*pi/180) 6.31*sin(pi/4) -6.31*cos(pi/4)], ...
    [9 5.035*sin(pi/4) -5.035*cos(pi/4)]);
[LwrShldrN,~] = geom.Loft.limingConic([5.9139 2.9739*sin(-pi/4) -2.9739*cos(-pi/4)], ...
    [18 6.31*sin(-pi/4) -6.31*cos(-pi/4)], ...
    [4.6328+(6.31-2.9739)*tan(27*pi/180) 6.31*sin(-pi/4) -6.31*cos(-pi/4)], ...
    [9 5.035*sin(-pi/4) -5.035*cos(-pi/4)]);

S1 = geom.NURBSSurface.loft({Upr, UprShldrP, MaxBp, LwrShldrP, Lwr}, 2);
S2 = geom.NURBSSurface.loft({Upr, UprShldrN, MaxBn, LwrShldrN, Lwr}, 2);

[S1g1, S2g1, info] = geom.matchSurfacePairTwoSeamG1Analytic(S1, S2);

fprintf('Two-seam solve: %d scalar constraints, %d scalar unknowns, residual %.3e\n', ...
    info.numScalarConstraints, info.numScalarUnknowns, info.residualNorm);

figure('Name', 'Two-seam simultaneous G1 analytic match', 'Color', 'w');

subplot(1,2,1); hold on;
title('Before simultaneous two-seam G1 match');
S1.plot(70,35); S2.plot(70,35);
localPlotAllSeamArrows(S1, S2, info.edgePairs, info.reverseAlongEdge);
xlabel('X'); ylabel('Y'); zlabel('Z'); axis equal; view(3); grid on;

subplot(1,2,2); hold on;
title('After simultaneous two-seam G1 match');
S1g1.plot(70,35); S2g1.plot(70,35);
localPlotAllSeamArrows(S1g1, S2g1, info.edgePairs, info.reverseAlongEdge);
xlabel('X'); ylabel('Y'); zlabel('Z'); axis equal; view(3); grid on;

sgtitle('Cowling half-patches: simultaneous two-seam analytic G1 solve');
end

function localPlotAllSeamArrows(S1, S2, edgePairs, reverseFlags)
for k = 1:2
    localPlotSeamArrows(S1, edgePairs{k,1}, [0.85 0.10 0.10], false);
    localPlotSeamArrows(S2, edgePairs{k,2}, [0.10 0.10 0.85], reverseFlags(k));
end
end

function localPlotSeamArrows(S, edgeName, colorVec, reverseAlong)
if nargin < 4, reverseAlong = false; end
nSamp = 7;
switch lower(edgeName)
    case 'u0'
        v = linspace(S.domainV(1), S.domainV(2), nSamp); if reverseAlong, v = fliplr(v); end
        u = S.domainU(1);
        P0 = zeros(nSamp,3); Dn = zeros(nSamp,3);
        for k = 1:nSamp
            D = S.derivatives(u, v(k), 1);
            P0(k,:) = D{1,1};
            Dn(k,:) = D{2,1};
        end
    case 'u1'
        v = linspace(S.domainV(1), S.domainV(2), nSamp); if reverseAlong, v = fliplr(v); end
        u = S.domainU(2);
        P0 = zeros(nSamp,3); Dn = zeros(nSamp,3);
        for k = 1:nSamp
            D = S.derivatives(u, v(k), 1);
            P0(k,:) = D{1,1};
            Dn(k,:) = -D{2,1};
        end
    case 'v0'
        u = linspace(S.domainU(1), S.domainU(2), nSamp); if reverseAlong, u = fliplr(u); end
        v = S.domainV(1);
        P0 = zeros(nSamp,3); Dn = zeros(nSamp,3);
        for k = 1:nSamp
            D = S.derivatives(u(k), v, 1);
            P0(k,:) = D{1,1};
            Dn(k,:) = D{1,2};
        end
    case 'v1'
        u = linspace(S.domainU(1), S.domainU(2), nSamp); if reverseAlong, u = fliplr(u); end
        v = S.domainV(2);
        P0 = zeros(nSamp,3); Dn = zeros(nSamp,3);
        for k = 1:nSamp
            D = S.derivatives(u(k), v, 1);
            P0(k,:) = D{1,1};
            Dn(k,:) = -D{1,2};
        end
    otherwise
        error('Unknown edge');
end
quiver3(P0(:,1), P0(:,2), P0(:,3), Dn(:,1), Dn(:,2), Dn(:,3), 0.2, ...
    'Color', colorVec, 'LineWidth', 1.4, 'MaxHeadSize', 0.6);
end
