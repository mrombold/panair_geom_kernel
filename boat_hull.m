%% demo_loft_liming_fuselage_sections.m
clear; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir)
    addpath(script_dir);
end
if ~exist('+geom', 'dir') && isempty(which('geom.NURBSCurve'))
    error(['Cannot find +geom package.\n' ...
           'Run this from the repo root or add the geometry kernel to path.']);
end
fprintf('geom package found at: %s\n', script_dir);


    

[lwr,~]=geom.Loft.limingConic([0 0 0], [40 0 20],[35 0 0], [35 0 9])
[upr,~]=geom.Loft.limingConic([0 0 10], [40 0 20],[30 0 10], [30 0 13])
[shldr,~]=geom.Loft.limingConic([0 0 9], [40 0 20],[30 0 9], [30 0 12])
[tan,~]=geom.Loft.limingConic([0 11 0], [40 0 0], [20 11 0], [20 7 0])
[width,~]=geom.Loft.limingConic([0 15 0], [40 0 0], [20 15 0], [20 12 0])
[MaxB,~]=geom.Loft.combinePlanarGuidesTo3D(width,upr)
[tancurve,~]=geom.Loft.combinePlanarGuidesTo3D(tan,upr)
[shldr3d,~]=geom.Loft.combinePlanarGuidesTo3D(tan,shldr)

lwr.plot()
%upr.plot()
%width.plot()
%tan.plot()
MaxB.plot()
%tancurve.plot()
shldr3d.plot()



% for sta=0:1:39
% 
%     [~,l]=geom.Loft.sampleCurveAtStation(lwr,sta)
%     [~,mb]=geom.Loft.sampleCurveAtStation(MaxB,sta)
%     [~,sh]=geom.Loft.sampleCurveAtStation(shldr3d,sta)
%     [~,t]=geom.Loft.sampleCurveAtStation(tancurve,sta)
% 
%     [frame0,~]=geom.Loft.limingConic(l,mb,t,sh)
%     frame0.plot()
% end

Sconic = geom.ConicSurface( ...
    'UpperGuide', MaxB, ...
    'LowerGuide', lwr, ...
    'TangencyGuide', tancurve, ...
    'ShoulderGuide', shldr3d, ...
    'SweepOrigin', [0 0 0], ...
    'SweepVector', [1 0 0], ...
    'StationRange', [0 39], ...
    'EnableSectionCache', true);


M = Sconic.isoMesh(101, 41);

hSurf = surf(M.X, M.Y, M.Z, ...
    'FaceAlpha', 0.75, ...
    'EdgeColor', [0.25 0.25 0.25], ...
    'FaceColor', [0.60 0.78 0.96]);


ax = gca;
ax.Clipping = 'off';