function demo_surface_tools_fill_g1()
%DEMO_SURFACE_TOOLS_FILL_G1 Exercise CATIA-like SurfaceTools fills.
%
% Put this file in examples/demo_surface_tools_fill_g1.m
% Put SurfaceTools.m in +geom/SurfaceTools.m
%
% Run:
%   demo_surface_tools_fill_g1

    clc;
    close all;

    thisFile = mfilename('fullpath');
    examplesDir = fileparts(thisFile);
    root = fileparts(examplesDir);
    addpath(root);

    fprintf('geom package found at: %s\n', root);
    fprintf('\n--- SurfaceTools two-edge fill demo ---\n');

    U = [0 0 0 0 1 1 1 1];
    V = [0 0 0 0 1 1 1 1];
    p = 3; q = 3;

    P1 = zeros(4,4,3);
    P2 = zeros(4,4,3);

    for i = 1:4
        u = (i-1)/3;
        for j = 1:4
            v = (j-1)/3;

            x = 3*u;
            y = 1.2*v;
            z = 0.18*sin(pi*u) + 0.10*v + 0.05*u*v;
            P1(i,j,:) = [x, y, z];

            x2 = 3*u;
            y2 = 1.8 + 1.2*v;
            z2 = 0.18*sin(pi*u) + 0.28 + 0.12*v - 0.03*u*v;
            P2(i,j,:) = [x2, y2, z2];
        end
    end

    S1 = geom.NURBSSurface(P1, p, q, U, V);
    S2 = geom.NURBSSurface(P2, p, q, U, V);

    e1 = geom.SurfaceTools.edge(S1, 'v1');
    e2 = geom.SurfaceTools.edge(S2, 'v0');

    % Basic fills
    Fg0 = geom.SurfaceTools.twoEdgeFill(e1, e2, ...
        'Continuity1','G0', ...
        'Continuity2','G0', ...
        'DegreeV',3);

    Fg1 = geom.SurfaceTools.twoEdgeFill(e1, e2, ...
        'Continuity1','G1', ...
        'Continuity2','G1', ...
        'DegreeV',3, ...
        'Scale1',0.65, ...
        'Scale2',0.65);

    Ffair = geom.SurfaceTools.twoEdgeFillFair(e1, e2, ...
        'Continuity1','G1', ...
        'Continuity2','G1', ...
        'DegreeV',5, ...
        'Scale1',0.65, ...
        'Scale2',0.65, ...
        'Fairness',0.35);

    % Variable tangent law fill
    Flaw = geom.SurfaceTools.twoEdgeFillLaw(e1, e2, ...
        'Continuity1','G1', ...
        'Continuity2','G1', ...
        'DegreeV',7, ...
        'ScaleLaw1',@(t) 0.45 + 0.25*sin(pi*t), ...
        'ScaleLaw2',@(t) 0.55 + 0.20*sin(pi*t), ...
        'Fairness',0.25);

    % Relaxed version, preserving boundary and tangent rows
    FlawRelax = geom.SurfaceTools.relaxInterior(Flaw, ...
        'Iterations',20, ...
        'Lambda',0.30, ...
        'PreserveTangencyRows',true);

    % Reports
    fprintf('\n--- Continuity reports ---\n');

    printReport('G0    -> S1', Fg0, 'v0', S1, 'v1');
    printReport('G0    -> S2', Fg0, 'v1', S2, 'v0');

    printReport('G1    -> S1', Fg1, 'v0', S1, 'v1');
    printReport('G1    -> S2', Fg1, 'v1', S2, 'v0');

    printReport('Fair  -> S1', Ffair, 'v0', S1, 'v1');
    printReport('Fair  -> S2', Ffair, 'v1', S2, 'v0');

    printReport('Law   -> S1', Flaw, 'v0', S1, 'v1');
    printReport('Law   -> S2', Flaw, 'v1', S2, 'v0');

    printReport('Relax -> S1', FlawRelax, 'v0', S1, 'v1');
    printReport('Relax -> S2', FlawRelax, 'v1', S2, 'v0');

    % Main visualization
    figure('Name','SurfaceTools G1 Fill');
    hold on; grid on; axis equal;
    title('Two-edge G1 fill');
    xlabel('X'); ylabel('Y'); zlabel('Z');

    S1.plot(25,25,'ShowCP',false,'FaceColor',[0.25 0.55 0.95],'Alpha',0.65);
    S2.plot(25,25,'ShowCP',false,'FaceColor',[0.25 0.95 0.55],'Alpha',0.65);
    Fg1.plot(25,15,'ShowCP',true,'FaceColor',[0.95 0.65 0.25],'Alpha',0.85);
    view(35,22);

    % Compare primary fill modes
    figure('Name','SurfaceTools Fill Comparison');

    plotCase(1,5,1,'G0 fill',S1,S2,Fg0,20,12);
    plotCase(1,5,2,'G1 fill',S1,S2,Fg1,20,12);
    plotCase(1,5,3,'Fair G1 fill',S1,S2,Ffair,20,16);
    plotCase(1,5,4,'Tangent law fill',S1,S2,Flaw,20,18);
    plotCase(1,5,5,'Law + relaxed fill',S1,S2,FlawRelax,20,18);

    % Control-net-only comparison
    figure('Name','Fill Control Nets');
    hold on; grid on; axis equal;
    title('Fill control nets');
    xlabel('X'); ylabel('Y'); zlabel('Z');

    Fg1.plot(5,5,'ShowCP',true,'FaceColor',[0.95 0.65 0.25],'Alpha',0.25);
    Flaw.plot(5,5,'ShowCP',true,'FaceColor',[0.65 0.65 0.95],'Alpha',0.25);
    FlawRelax.plot(5,5,'ShowCP',true,'FaceColor',[0.95 0.35 0.35],'Alpha',0.25);
    view(35,22);

    fprintf('\nDemo complete.\n');

end


function printReport(label, Sa, edgeA, Sb, edgeB)
    r = geom.SurfaceTools.edgeReport(Sa, edgeA, Sb, edgeB, 51);
    fprintf('%-12s gap max %.3e, rms %.3e, normal angle max %.3f deg, rms %.3f deg\n', ...
        label, r.maxGap, r.rmsGap, r.maxNormalAngleDeg, r.rmsNormalAngleDeg);
end


function plotCase(nr,nc,idx,name,S1,S2,F,ns,nf)
    subplot(nr,nc,idx);
    hold on; grid on; axis equal;
    title(name);

    S1.plot(ns,ns,'ShowCP',false,'FaceColor',[0.25 0.55 0.95],'Alpha',0.50);
    S2.plot(ns,ns,'ShowCP',false,'FaceColor',[0.25 0.95 0.55],'Alpha',0.50);
    F.plot(ns,nf,'ShowCP',true,'FaceColor',[0.95 0.65 0.25],'Alpha',0.85);

    xlabel('X'); ylabel('Y'); zlabel('Z');
    view(35,22);
end