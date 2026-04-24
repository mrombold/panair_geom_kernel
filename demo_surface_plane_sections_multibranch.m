function demo_surface_plane_sections_multibranch()
% DEMO_SURFACE_PLANE_SECTIONS_MULTIBRANCH
% Demonstrates multi-branch plane/surface section extraction.
%
% Requires:
%   +geom/intersectCurvePlaneBisection.m
%   +geom/surfacePlaneSection_multibranch.m
%
% Optional:
%   If you also have the original surfacePlaneSection.m, you can compare
%   single-branch behavior against multi-branch behavior.

clear; close all; clc;

%% ------------------------------------------------------------------------
% Build a simple lofted surface that can give multiple z-plane intersections
%% ------------------------------------------------------------------------
x = linspace(0,1,21).';

sec1 = geom.NURBSCurve.globalInterp([0.0+0*x, 0.70*sin(pi*x),  0.10*cos(pi*x)], 3, 'centripetal');
sec2 = geom.NURBSCurve.globalInterp([1.0+0*x, 0.90*sin(pi*x),  0.12+0.08*cos(pi*x)], 3, 'centripetal');
sec3 = geom.NURBSCurve.globalInterp([2.0+0*x, 1.00*sin(pi*x),  0.18+0.06*cos(pi*x)], 3, 'centripetal');
sec4 = geom.NURBSCurve.globalInterp([3.0+0*x, 0.60*sin(pi*x),  0.08*cos(pi*x)], 3, 'centripetal');

S = geom.NURBSSurface.loft({sec1, sec2, sec3, sec4}, 3, 'centripetal');

%% ------------------------------------------------------------------------
% Example planes
%% ------------------------------------------------------------------------
xStation = 1.5;
yButtock = 0.35;
zWater   = 0.10;

%% ------------------------------------------------------------------------
% Extract sections
%% ------------------------------------------------------------------------
Cstation = geom.surfacePlaneSection_multibranch(S, 'x', xStation, ...
    'Direction', 'v', ...
    'MultiBranchAction', 'warn_primary');

Cbuttock = geom.surfacePlaneSection_multibranch(S, 'y', yButtock, ...
    'Direction', 'u', ...
    'MultiBranchAction', 'warn_primary');

% This one is the interesting case: ask for all branches
[Cwater, infoWater] = geom.surfacePlaneSection_multibranch(S, 'z', zWater, ...
    'Direction', 'u', ...
    'MultiBranchAction', 'all');

Ccanted = geom.surfacePlaneSection_multibranch(S, [1.5 0 0], [1 0.25 0.15], ...
    'Offset', 0.0, ...
    'Direction', 'auto', ...
    'MultiBranchAction', 'warn_primary');

%% ------------------------------------------------------------------------
% Diagnostics
%% ------------------------------------------------------------------------
fprintf('=== Multi-branch section demo ===\n');
fprintf('Station X=%g -> one primary curve\n', xStation);
fprintf('Buttock Y=%g -> one primary curve\n', yButtock);

if iscell(Cwater)
    fprintf('Waterline Z=%g -> %d branches returned\n', zWater, numel(Cwater));
else
    fprintf('Waterline Z=%g -> single branch returned\n', zWater);
end

fprintf('infoWater.numBranches = %d\n', infoWater.numBranches);
fprintf('Chosen direction      = %s\n', infoWater.direction);
fprintf('Strategy              = %s\n', infoWater.strategy);

%% ------------------------------------------------------------------------
% Plot
%% ------------------------------------------------------------------------
figure('Color','w','Name','Surface plane sections (multi-branch)');
hold on; grid on; axis equal; view(3);

S.plot(50, 26, 'ShowCP', false, 'Alpha', 0.55, 'FaceColor', [0.70 0.80 0.98]);

% Station
Cstation.plot(200, 'Color', [0.90 0.10 0.10], 'LineWidth', 2.0);

% Buttock
Cbuttock.plot(200, 'Color', [0.10 0.60 0.10], 'LineWidth', 2.0);

% Waterline branches
if iscell(Cwater)
    waterColors = [ ...
        0.10 0.10 0.90
        0.00 0.70 0.90
        0.55 0.20 0.90
        0.10 0.10 0.50];
    for k = 1:numel(Cwater)
        ck = waterColors(1+mod(k-1,size(waterColors,1)), :);
        Cwater{k}.plot(220, 'Color', ck, 'LineWidth', 2.2);
    end
else
    Cwater.plot(220, 'Color', [0.10 0.10 0.90], 'LineWidth', 2.2);
end

% Canted plane
Ccanted.plot(220, 'Color', [0.80 0.40 0.00], 'LineWidth', 2.0);

xlabel('X'); ylabel('Y'); zlabel('Z');
title('NURBS surface section curves from planes (multi-branch extraction)');

legendEntries = {'Surface', ...
                 sprintf('Station X=%g', xStation), ...
                 sprintf('Buttock Y=%g', yButtock)};

if iscell(Cwater)
    for k = 1:numel(Cwater)
        legendEntries{end+1} = sprintf('Waterline Z=%g (branch %d)', zWater, k); %#ok<AGROW>
    end
else
    legendEntries{end+1} = sprintf('Waterline Z=%g', zWater);
end

legendEntries{end+1} = 'Canted plane';
legend(legendEntries, 'Location', 'best');

%% ------------------------------------------------------------------------
% Optional branch point plots for debugging
%% ------------------------------------------------------------------------
if isfield(infoWater, 'branchPoints') && ~isempty(infoWater.branchPoints)
    figure('Color','w','Name','Waterline branch sample points');
    hold on; grid on; axis equal; view(3);
    S.plot(50, 26, 'ShowCP', false, 'Alpha', 0.35, 'FaceColor', [0.85 0.90 1.00]);

    for k = 1:numel(infoWater.branchPoints)
        Pk = infoWater.branchPoints{k};
        plot3(Pk(:,1), Pk(:,2), Pk(:,3), 'o-', 'LineWidth', 1.5, 'MarkerSize', 5);
    end

    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('Sample points used to build waterline branches at Z=%g', zWater));
end

end
