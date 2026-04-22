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

% -------------------------------------------------------------------------
% User-editable Liming conic construction points
% -------------------------------------------------------------------------


[Upr,~]=geom.Loft.limingConic([4.6328 0 2.9060], [18 0 6.31],[4.6328+(6.31-2.9060)*tan(24*pi/180) 0 6.31], [9 0 5.35])
[MaxBp,~]=geom.Loft.limingConic([5.7352 2.941 0], [18 6.31 0],[5.7352+(6.31-2.941)*tan(30*pi/180) 6.31 0], [9 5.14 0])
[MaxBn,~]=geom.Loft.limingConic([5.7352 -2.941 0], [18 -6.31 0],[5.7352+(6.31-2.941)*tan(30*pi/180) -6.31 0], [9 -5.14 0])
[Lwr,~]=geom.Loft.limingConic([6.1476 0 -2.9739], [18 0 -6.31],[6.1476+(6.31-2.9739)*tan(36*pi/180) 0 -6.31], [9 0 -4.93])
[UprShldrP,~]=geom.Loft.limingConic([4.8509 2.9237*sin(pi/4) 2.9237*cos(pi/4)], ...
    [18 6.31*sin(pi/4) 6.31*cos(pi/4)], ...
    [4.6328+(6.31-2.9237)*tan(27*pi/180) 6.31*sin(pi/4) 6.31*cos(pi/4)], ...
    [9 5.35*sin(pi/4) 5.35*cos(pi/4)])
[UprShldrN,~]=geom.Loft.limingConic([4.8509 2.9237*sin(-pi/4) 2.9237*cos(-pi/4)], ...
    [18 6.31*sin(-pi/4) 6.31*cos(-pi/4)], ...
    [4.6328+(6.31-2.9237)*tan(27*pi/180) 6.31*sin(-pi/4) 6.31*cos(-pi/4)], ...
    [9 5.35*sin(-pi/4) 5.35*cos(-pi/4)])

[LwrShldrP,~]=geom.Loft.limingConic([5.9139 2.9739*sin(pi/4) -2.9739*cos(pi/4)], ...
    [18 6.31*sin(pi/4) -6.31*cos(pi/4)], ...
    [4.6328+(6.31-2.9739)*tan(27*pi/180) 6.31*sin(pi/4) -6.31*cos(pi/4)], ...
    [9 5.035*sin(pi/4) -5.035*cos(pi/4)])

[LwrShldrN,~]=geom.Loft.limingConic([5.9139 2.9739*sin(-pi/4) -2.9739*cos(-pi/4)], ...
    [18 6.31*sin(-pi/4) -6.31*cos(-pi/4)], ...
    [4.6328+(6.31-2.9739)*tan(27*pi/180) 6.31*sin(-pi/4) -6.31*cos(-pi/4)], ...
    [9 5.035*sin(-pi/4) -5.035*cos(-pi/4)])


%[MaxBTop,~]=geom.Loft.limingConic([0 5 0], [100 10 0],[25 10 0], [25 8 0])
%[UprShldrXZ,~]=geom.Loft.limingConic([0 0 38], [100 0 52],[25 0 52], [25 0 47])
%[UprShldrYZ,~]=geom.Loft.limingConic([0 4 0], [100 8 0],[25 8 0], [25 7 0])







%[MaxB,~]=geom.Loft.combinePlanarGuidesTo3D(MaxBTop,MaxBSide)
%[UprShldr,~]=geom.Loft.combinePlanarGuidesTo3D(UprShldrTop,UprShldrSide)
%sta=0
%%[~,temp1]=geom.Loft.sampleCurveAtStation(Upr,sta)
%[~,temp2]=geom.Loft.sampleCurveAtStation(MaxBSide,sta)
%[~,temp3]=geom.Loft.sampleCurveAtStation(MaxBTop,sta)
%[~,temp4]=geom.Loft.sampleCurveAtStation(UprShldrSide,sta)
%[~,temp5]=geom.Loft.sampleCurveAtStation(UprShldrTop,sta)

%staP0=[temp1]
%staP1=[temp2(1) temp3(2) temp2(3)]
%staT=[temp1(1) temp3(2) temp1(3)]
%staS=[temp4(1) temp5(2) temp4(3)]

%[frame0,~]=geom.Loft.limingConic(staP0,staP1,staT,staS)



S1=geom.NURBSSurface.loft({Upr, UprShldrP, MaxBp, LwrShldrP, Lwr},2)
S2=geom.NURBSSurface.loft({Upr, UprShldrN, MaxBn, LwrShldrN, Lwr},2)

% -------------------------------------------------------------------------
% Plot 3D construction
% -------------------------------------------------------------------------

figure('Name','Fuselage curves')
hold on
Upr.plot()
MaxBp.plot()
MaxBn.plot()
Lwr.plot()
UprShldrP.plot()
UprShldrN.plot()
LwrShldrP.plot()
LwrShldrN.plot()
S1.plot()
S2.plot()
ax = gca;
ax.Clipping = 'off';

fprintf('\nDone.\n');