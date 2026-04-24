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
[Upr,~]=geom.Loft.limingConic([0 0 40], [100 0 55],[25 0 55], [25 0 50])
[MaxBSide,~]=geom.Loft.limingConic([0 0 35], [100 0 40],[25 0 40], [25 0 36.5])
[MaxBTop,~]=geom.Loft.limingConic([0 5 0], [100 10 0],[25 10 0], [25 8 0])
[UprShldrSide,~]=geom.Loft.limingConic([0 0 38], [100 0 52],[25 0 52], [25 0 47])
[UprShldrTop,~]=geom.Loft.limingConic([0 4 0], [100 8 0],[25 8 0], [25 7 0])

[MaxB,~]=geom.Loft.combinePlanarGuidesTo3D(MaxBTop,MaxBSide)
[UprShldr,~]=geom.Loft.combinePlanarGuidesTo3D(UprShldrTop,UprShldrSide)
[TangencyGuide,~]=geom.Loft.combinePlanarGuidesTo3D(MaxBTop,Upr)


sta=0
[~,temp1]=geom.Loft.sampleCurveAtStation(Upr,sta)
[~,temp2]=geom.Loft.sampleCurveAtStation(MaxBSide,sta)
[~,temp3]=geom.Loft.sampleCurveAtStation(MaxBTop,sta)
[~,temp4]=geom.Loft.sampleCurveAtStation(UprShldrSide,sta)
[~,temp5]=geom.Loft.sampleCurveAtStation(UprShldrTop,sta)

staP0=[temp1]
staP1=[temp2(1) temp3(2) temp2(3)]
staT=[temp1(1) temp3(2) temp1(3)]
staS=[temp4(1) temp5(2) temp4(3)]

[frame0,~]=geom.Loft.limingConic(staP0,staP1,staT,staS)

sta=20
[~,temp1]=geom.Loft.sampleCurveAtStation(Upr,sta)
[~,temp2]=geom.Loft.sampleCurveAtStation(MaxBSide,sta)
[~,temp3]=geom.Loft.sampleCurveAtStation(MaxBTop,sta)
[~,temp4]=geom.Loft.sampleCurveAtStation(UprShldrSide,sta)
[~,temp5]=geom.Loft.sampleCurveAtStation(UprShldrTop,sta)

staP0=[temp1]
staP1=[temp2(1) temp3(2) temp2(3)]
staT=[temp1(1) temp3(2) temp1(3)]
staS=[temp4(1) temp5(2) temp4(3)]

[frame20,~]=geom.Loft.limingConic(staP0,staP1,staT,staS)

sta=40
[~,temp1]=geom.Loft.sampleCurveAtStation(Upr,sta)
[~,temp2]=geom.Loft.sampleCurveAtStation(MaxBSide,sta)
[~,temp3]=geom.Loft.sampleCurveAtStation(MaxBTop,sta)
[~,temp4]=geom.Loft.sampleCurveAtStation(UprShldrSide,sta)
[~,temp5]=geom.Loft.sampleCurveAtStation(UprShldrTop,sta)

staP0=[temp1]
staP1=[temp2(1) temp3(2) temp2(3)]
staT=[temp1(1) temp3(2) temp1(3)]
staS=[temp4(1) temp5(2) temp4(3)]

[frame40,~]=geom.Loft.limingConic(staP0,staP1,staT,staS)

sta=60
[~,temp1]=geom.Loft.sampleCurveAtStation(Upr,sta)
[~,temp2]=geom.Loft.sampleCurveAtStation(MaxBSide,sta)
[~,temp3]=geom.Loft.sampleCurveAtStation(MaxBTop,sta)
[~,temp4]=geom.Loft.sampleCurveAtStation(UprShldrSide,sta)
[~,temp5]=geom.Loft.sampleCurveAtStation(UprShldrTop,sta)

staP0=[temp1]
staP1=[temp2(1) temp3(2) temp2(3)]
staT=[temp1(1) temp3(2) temp1(3)]
staS=[temp4(1) temp5(2) temp4(3)]

[frame60,~]=geom.Loft.limingConic(staP0,staP1,staT,staS)

sta=80
[~,temp1]=geom.Loft.sampleCurveAtStation(Upr,sta)
[~,temp2]=geom.Loft.sampleCurveAtStation(MaxBSide,sta)
[~,temp3]=geom.Loft.sampleCurveAtStation(MaxBTop,sta)
[~,temp4]=geom.Loft.sampleCurveAtStation(UprShldrSide,sta)
[~,temp5]=geom.Loft.sampleCurveAtStation(UprShldrTop,sta)

staP0=[temp1]
staP1=[temp2(1) temp3(2) temp2(3)]
staT=[temp1(1) temp3(2) temp1(3)]
staS=[temp4(1) temp5(2) temp4(3)]

[frame80,~]=geom.Loft.limingConic(staP0,staP1,staT,staS)

sta=100
[~,temp1]=geom.Loft.sampleCurveAtStation(Upr,sta)
[~,temp2]=geom.Loft.sampleCurveAtStation(MaxBSide,sta)
[~,temp3]=geom.Loft.sampleCurveAtStation(MaxBTop,sta)
[~,temp4]=geom.Loft.sampleCurveAtStation(UprShldrSide,sta)
[~,temp5]=geom.Loft.sampleCurveAtStation(UprShldrTop,sta)

staP0=[temp1]
staP1=[temp2(1) temp3(2) temp2(3)]
staT=[temp1(1) temp3(2) temp1(3)]
staS=[temp4(1) temp5(2) temp4(3)]

[frame100,~]=geom.Loft.limingConic(staP0,staP1,staT,staS)


S=geom.NURBSSurface.loft({frame0 frame20 frame40 frame60 frame80 frame100},3)
%Sgordon=geom.NURBSSurface.gordon({frame0 frame20 frame40 frame60 frame80 frame100},{Upr, UprShldr, MaxB},2,3)
Sconic = geom.ConicSurface( ...
    'UpperGuide', Upr, ...
    'LowerGuide', MaxB, ...
    'TangencyGuide', TangencyGuide, ...
    'ShoulderGuide', UprShldr, ...
    'SweepOrigin', [0 0 0], ...
    'SweepVector', [1 0 0], ...
    'StationRange', [0 100], ...
    'EnableSectionCache', true);


M = Sconic.isoMesh(101, 41);

%[Lwr,~]=geom.Loft.limingConic([0 0 30], [100 0 8],[25 0 8], [25 0 15])



% -------------------------------------------------------------------------
% Plot 3D construction
% -------------------------------------------------------------------------

figure('Name','Fuselage curves')
frame0.plot()
hold on
Upr.plot()
MaxB.plot()
UprShldr.plot()
frame20.plot()
frame40.plot()
frame60.plot()
frame80.plot()
frame100.plot()
S.plot()
Sconic_nurbs.plot()

hSurf = surf(M.X, M.Y, M.Z, ...
    'FaceAlpha', 0.75, ...
    'EdgeColor', [0.25 0.25 0.25], ...
    'FaceColor', [0.60 0.78 0.96]);

%Sgordon.plot()
ax = gca;
ax.Clipping = 'off';

fprintf('\nDone.\n');