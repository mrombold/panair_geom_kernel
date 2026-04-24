function demo_surface_plane_sections()
clear; close all; clc;

x = linspace(0,1,21).';
sec1 = geom.NURBSCurve.globalInterp([0*x,        0.7*sin(pi*x), 0.10*cos(pi*x)], 3, 'centripetal');
sec2 = geom.NURBSCurve.globalInterp([1.0+0*x,    0.9*sin(pi*x), 0.12+0.08*cos(pi*x)], 3, 'centripetal');
sec3 = geom.NURBSCurve.globalInterp([2.0+0*x,    1.0*sin(pi*x), 0.18+0.06*cos(pi*x)], 3, 'centripetal');
sec4 = geom.NURBSCurve.globalInterp([3.0+0*x,    0.6*sin(pi*x), 0.08*cos(pi*x)], 3, 'centripetal');

S = geom.NURBSSurface.loft({sec1,sec2,sec3,sec4}, 3, 'centripetal');

Csta = geom.surfacePlaneSection(S, 'x', 1.5, 'Direction', 'v');
Cbut = geom.surfacePlaneSection(S, 'y', 0.35, 'Direction', 'u');
Cwat = geom.surfacePlaneSection(S, 'z', 0.10, 'Direction', 'u');
Ccanted = geom.surfacePlaneSection(S, [1.5 0 0], [1 0.25 0.15], ...
    'Offset', 0.0, 'Direction', 'auto');

figure('Color','w','Name','Surface plane sections');
hold on; grid on; axis equal; view(3);

S.plot(40, 24, 'ShowCP', false, 'Alpha', 0.55, 'FaceColor', [0.65 0.75 0.95]);
Csta.plot(200, 'Color', [0.90 0.10 0.10], 'LineWidth', 2.0);
Cbut.plot(200, 'Color', [0.10 0.60 0.10], 'LineWidth', 2.0);
Cwat.plot(200, 'Color', [0.15 0.15 0.85], 'LineWidth', 2.0);
Ccanted.plot(200, 'Color', [0.80 0.40 0.00], 'LineWidth', 2.0);

legend({'Surface','Station X=1.5','Buttock Y=0.35','Waterline Z=0.10','Canted plane'}, ...
    'Location','best');
title('NURBS surface section curves from planes');
xlabel('X'); ylabel('Y'); zlabel('Z');
end
