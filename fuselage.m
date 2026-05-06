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


    

[fwd,~]=geom.Loft.limingConic([0 0 25], [58 0 5],[15 0 5], [22 0 12])
[mid,~]=geom.Loft.limingConic([58 0 5], [110 0 20],[86 0 5], [84 0 10])
[aft,~]=geom.Loft.limingConic([110 0 20], [155 0 25],[118  0 25], [120 0 23])

fwd.plot()
mid.plot()
aft.plot()
axis equal




