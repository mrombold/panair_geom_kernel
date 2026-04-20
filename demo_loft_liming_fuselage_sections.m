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
P0 = [0.0, 0.0, 1.0];
P1 = [1.2, 0.0, 0.0];
T  = [0.0, 0.0, 0.0];
S  = [0.42, 0.0, 0.38];

fprintf('\n--- Build single Liming conic ---\n');

[C, meta] = geom.Loft.limingConic(P0, P1, T, S);

[Upr,~]=geom.Loft.limingConic([0 0 40], [100 0 55],[25 0 55], [25 0 50])
[MaxBSide,~]=geom.Loft.limingConic([0 0 35], [100 0 40],[25 0 40], [25 0 36.5])
[MaxBTop,~]=geom.Loft.limingConic([0 5 0], [100 10 0],[25 10 0], [25 8 0])
[UprShldrSide,~]=geom.Loft.limingConic([0 0 38], [100 0 52],[25 0 52], [25 0 47])
[UprShldrTop,~]=geom.Loft.limingConic([0 4 0], [100 8 0],[25 8 0], [25 7 0])

sta0=0
[~,temp1]=geom.Loft.sampleCurveAtStation(Upr,sta0)
[~,temp2]=geom.Loft.sampleCurveAtStation(MaxBSide,sta0)
[~,temp3]=geom.Loft.sampleCurveAtStation(MaxBTop,sta0)
[~,temp4]=geom.Loft.sampleCurveAtStation(UprShldrSide,sta0)
[~,temp5]=geom.Loft.sampleCurveAtStation(UprShldrTop,sta0)

staP0=[temp1]
staP1=[temp2(1) temp3(2) temp2(3)]
staT=[temp1(1) temp3(2) temp1(3)]
staS=[temp4(1) temp5(2) temp4(3)]

[frame0,~]=geom.Loft.limingConic(staP0,staP1,staT,staS)

sta20=20
[~,temp1]=geom.Loft.sampleCurveAtStation(Upr,sta0)
[~,temp2]=geom.Loft.sampleCurveAtStation(MaxBSide,sta0)
[~,temp3]=geom.Loft.sampleCurveAtStation(MaxBTop,sta0)
[~,temp4]=geom.Loft.sampleCurveAtStation(UprShldrSide,sta0)
[~,temp5]=geom.Loft.sampleCurveAtStation(UprShldrTop,sta0)

staP0=[temp1]
staP1=[temp2(1) temp3(2) temp2(3)]
staT=[temp1(1) temp3(2) temp1(3)]
staS=[temp4(1) temp5(2) temp4(3)]

[frame20,~]=geom.Loft.limingConic(staP0,staP1,staT,staS)

sta40=40
[~,temp1]=geom.Loft.sampleCurveAtStation(Upr,sta0)
[~,temp2]=geom.Loft.sampleCurveAtStation(MaxBSide,sta0)
[~,temp3]=geom.Loft.sampleCurveAtStation(MaxBTop,sta0)
[~,temp4]=geom.Loft.sampleCurveAtStation(UprShldrSide,sta0)
[~,temp5]=geom.Loft.sampleCurveAtStation(UprShldrTop,sta0)

staP0=[temp1]
staP1=[temp2(1) temp3(2) temp2(3)]
staT=[temp1(1) temp3(2) temp1(3)]
staS=[temp4(1) temp5(2) temp4(3)]

[frame40,~]=geom.Loft.limingConic(staP0,staP1,staT,staS)

sta60=60
[~,temp1]=geom.Loft.sampleCurveAtStation(Upr,sta0)
[~,temp2]=geom.Loft.sampleCurveAtStation(MaxBSide,sta0)
[~,temp3]=geom.Loft.sampleCurveAtStation(MaxBTop,sta0)
[~,temp4]=geom.Loft.sampleCurveAtStation(UprShldrSide,sta0)
[~,temp5]=geom.Loft.sampleCurveAtStation(UprShldrTop,sta0)

staP0=[temp1]
staP1=[temp2(1) temp3(2) temp2(3)]
staT=[temp1(1) temp3(2) temp1(3)]
staS=[temp4(1) temp5(2) temp4(3)]

[frame60,~]=geom.Loft.limingConic(staP0,staP1,staT,staS)

sta80=80
[~,temp1]=geom.Loft.sampleCurveAtStation(Upr,sta0)
[~,temp2]=geom.Loft.sampleCurveAtStation(MaxBSide,sta0)
[~,temp3]=geom.Loft.sampleCurveAtStation(MaxBTop,sta0)
[~,temp4]=geom.Loft.sampleCurveAtStation(UprShldrSide,sta0)
[~,temp5]=geom.Loft.sampleCurveAtStation(UprShldrTop,sta0)

staP0=[temp1]
staP1=[temp2(1) temp3(2) temp2(3)]
staT=[temp1(1) temp3(2) temp1(3)]
staS=[temp4(1) temp5(2) temp4(3)]

[frame80,~]=geom.Loft.limingConic(staP0,staP1,staT,staS)

sta100=100
[~,temp1]=geom.Loft.sampleCurveAtStation(Upr,sta0)
[~,temp2]=geom.Loft.sampleCurveAtStation(MaxBSide,sta0)
[~,temp3]=geom.Loft.sampleCurveAtStation(MaxBTop,sta0)
[~,temp4]=geom.Loft.sampleCurveAtStation(UprShldrSide,sta0)
[~,temp5]=geom.Loft.sampleCurveAtStation(UprShldrTop,sta0)

staP0=[temp1]
staP1=[temp2(1) temp3(2) temp2(3)]
staT=[temp1(1) temp3(2) temp1(3)]
staS=[temp4(1) temp5(2) temp4(3)]

[frame100,~]=geom.Loft.limingConic(staP0,staP1,staT,staS)


[Lwr,~]=geom.Loft.limingConic([0 0 30], [100 0 8],[25 0 8], [25 0 15])

% -------------------------------------------------------------------------
% Plot 3D construction
% -------------------------------------------------------------------------

figure('Name','Fuselage curves')
frame0.plot()
hold on
%Upr.plot()
%Lwr.plot()
%MaxBSide.plot()
%MaxBTop.plot()
%UprShldrSide.plot()
%UprShldrTop.plot()
frame20.plot()
frame40.plot()
frame60.plot()
frame80.plot()
frame100.plot()

fprintf('\nDone.\n');