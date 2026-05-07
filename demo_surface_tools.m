function demo_surface_tools()
%DEMO_SURFACE_TOOLS Exercise CATIA-like geom.SurfaceTools utilities.
%
% Put this file in examples/demo_surface_tools.m
% Put SurfaceTools.m in +geom/SurfaceTools.m
%
% Run:
%   demo_surface_tools

    clc;
    close all;

    thisFile = mfilename('fullpath');
    examplesDir = fileparts(thisFile);
    root = fileparts(examplesDir);
    addpath(root);

    fprintf('geom package found at: %s\n', root);
    fprintf('\n=== SurfaceTools Demo ===\n');

    [S1, S2] = buildParentSurfaces();

    e1 = geom.SurfaceTools.edge(S1, 'v1');
    e2 = geom.SurfaceTools.edge(S2, 'v0');

    fprintf('\n--- 1. Two-edge fill / blend tools ---\n');

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

    Flaw = geom.SurfaceTools.twoEdgeFillLaw(e1, e2, ...
        'Continuity1','G1', ...
        'Continuity2','G1', ...
        'DegreeV',7, ...
        'ScaleLaw1',@(t) 0.45 + 0.25*sin(pi*t), ...
        'ScaleLaw2',@(t) 0.55 + 0.20*sin(pi*t), ...
        'Fairness',0.25);

    FlawRelax = geom.SurfaceTools.relaxInterior(Flaw, ...
        'Iterations',20, ...
        'Lambda',0.30, ...
        'PreserveTangencyRows',true);

    fprintf('\n--- 2. Continuity reports ---\n');

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

    if ismethod('geom.SurfaceTools','edgeG1Report')
        fprintf('\n--- 3. Detailed G1 reports ---\n');
        printG1Report('G1    -> S1', Fg1, 'v0', S1, 'v1');
        printG1Report('G1    -> S2', Fg1, 'v1', S2, 'v0');
        printG1Report('Relax -> S1', FlawRelax, 'v0', S1, 'v1');
        printG1Report('Relax -> S2', FlawRelax, 'v1', S2, 'v0');
    else
        fprintf('\n--- 3. Detailed G1 reports skipped: edgeG1Report not found ---\n');
    end

    fprintf('\n--- 4. Match-surface / G1 seam conditioning ---\n');

    Fmatched = [];
    if ismethod('geom.SurfaceTools','matchG1EdgeRobust')
        [Fmatched0, rpt0] = geom.SurfaceTools.matchG1EdgeRobust( ...
            Fg0, 'v0', S1, 'v1', ...
            'MagnitudeMode','target', ...
            'Scale',1.0);

        [Fmatched, rpt1] = geom.SurfaceTools.matchG1EdgeRobust( ...
            Fmatched0, 'v1', S2, 'v0', ...
            'MagnitudeMode','target', ...
            'Scale',1.0);

        fprintf('matchG1EdgeRobust v0->S1: gap %.3e, normal angle %.3f deg\n', ...
            rpt0.maxGap, rpt0.maxNormalAngleDeg);
        fprintf('matchG1EdgeRobust v1->S2: gap %.3e, normal angle %.3f deg\n', ...
            rpt1.maxGap, rpt1.maxNormalAngleDeg);

        printReport('Matched -> S1', Fmatched, 'v0', S1, 'v1');
        printReport('Matched -> S2', Fmatched, 'v1', S2, 'v0');
    else
        fprintf('matchG1EdgeRobust not found. Add it to SurfaceTools to exercise this section.\n');
    end

    fprintf('\n--- 5. Visualization ---\n');

    figure('Name','SurfaceTools Main G1 Fill');
    hold on; grid on; axis equal;
    title('Two-edge G1 fill');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    plotParents(S1,S2,25);
    Fg1.plot(25,15,'ShowCP',true,'FaceColor',[0.95 0.65 0.25],'Alpha',0.85);
    view(35,22);

    figure('Name','G0 Fill');
    plotSingleCase('G0 fill',S1,S2,Fg0,20,12);
    
    figure('Name','G1 Fill');
    plotSingleCase('G1 fill',S1,S2,Fg1,20,12);
    
    figure('Name','Fair G1 Fill');
    plotSingleCase('Fair G1 fill',S1,S2,Ffair,20,16);
    
    figure('Name','Tangent Law Fill');
    plotSingleCase('Tangent law fill',S1,S2,Flaw,20,18);
    
    figure('Name','Law + Relaxed Fill');
    plotSingleCase('Law + relaxed fill',S1,S2,FlawRelax,20,18);
    
    if ~isempty(Fmatched)
        figure('Name','Matched G1 from G0');
        plotSingleCase('Matched G1 from G0',S1,S2,Fmatched,20,12);
    end

    figure('Name','SurfaceTools Control Nets');
    hold on; grid on; axis equal;
    title('Fill control nets');
    xlabel('X'); ylabel('Y'); zlabel('Z');

    Fg1.plot(5,5,'ShowCP',true,'FaceColor',[0.95 0.65 0.25],'Alpha',0.20);
    Ffair.plot(5,5,'ShowCP',true,'FaceColor',[0.25 0.85 0.85],'Alpha',0.20);
    Flaw.plot(5,5,'ShowCP',true,'FaceColor',[0.65 0.65 0.95],'Alpha',0.20);
    FlawRelax.plot(5,5,'ShowCP',true,'FaceColor',[0.95 0.35 0.35],'Alpha',0.20);

    if ~isempty(Fmatched)
        Fmatched.plot(5,5,'ShowCP',true,'FaceColor',[0.35 0.95 0.35],'Alpha',0.20);
    end

    view(35,22);


    fprintf('\n--- 6. Four-edge fill demo ---\n');
    
    % Build a rectangular hole from four surrounding parent surfaces.
    [Sb, St, Sl, Sr] = buildFourEdgeParents();
    
    eBottom = geom.SurfaceTools.edge(Sb, 'v1');
    eTop    = geom.SurfaceTools.edge(St, 'v0');
    eLeft   = geom.SurfaceTools.edge(Sl, 'u1');
    eRight  = geom.SurfaceTools.edge(Sr, 'u0');
    
    F4g0 = geom.SurfaceTools.fourEdgeFill( ...
        eBottom, eTop, eLeft, eRight, ...
        'ContinuityV0','G0', ...
        'ContinuityV1','G0', ...
        'ContinuityU0','G0', ...
        'ContinuityU1','G0');
    
    F4g1 = geom.SurfaceTools.fourEdgeFill( ...
        eBottom, eTop, eLeft, eRight, ...
        'ContinuityV0','G1', ...
        'ContinuityV1','G1', ...
        'ContinuityU0','G1', ...
        'ContinuityU1','G1', ...
        'ScaleV0',0.60, ...
        'ScaleV1',0.60, ...
        'ScaleU0',0.60, ...
        'ScaleU1',0.60);
    
    fprintf('\nFour-edge G0 reports:\n');
    printReport('F4G0 -> bottom', F4g0, 'v0', Sb, 'v1');
    printReport('F4G0 -> top',    F4g0, 'v1', St, 'v0');
    printReport('F4G0 -> left',   F4g0, 'u0', Sl, 'u1');
    printReport('F4G0 -> right',  F4g0, 'u1', Sr, 'u0');
    
    fprintf('\nFour-edge G1 reports:\n');
    printReport('F4G1 -> bottom', F4g1, 'v0', Sb, 'v1');
    printReport('F4G1 -> top',    F4g1, 'v1', St, 'v0');
    printReport('F4G1 -> left',   F4g1, 'u0', Sl, 'u1');
    printReport('F4G1 -> right',  F4g1, 'u1', Sr, 'u0');
    
    figure('Name','Four-edge G0 Fill');
    plotFourEdgeCase('Four-edge G0 fill', Sb, St, Sl, Sr, F4g0);
    
    figure('Name','Four-edge G1 Fill');
    plotFourEdgeCase('Four-edge G1 fill', Sb, St, Sl, Sr, F4g1);

    fprintf('\nFour-edge G1 interior-only reports, endpoints trimmed:\n');
    printInteriorG1Report('F4G1 int bottom', F4g1, 'v0', Sb, 'v1');
    printInteriorG1Report('F4G1 int top',    F4g1, 'v1', St, 'v0');
    printInteriorG1Report('F4G1 int left',   F4g1, 'u0', Sl, 'u1');
    printInteriorG1Report('F4G1 int right',  F4g1, 'u1', Sr, 'u0');



    fprintf('\nDemo complete.\n');

end


function [S1, S2] = buildParentSurfaces()
    U = [0 0 0 0 1 1 1 1];
    V = [0 0 0 0 1 1 1 1];
    p = 3;
    q = 3;

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
end


function printReport(label, Sa, edgeA, Sb, edgeB)
    r = geom.SurfaceTools.edgeReport(Sa, edgeA, Sb, edgeB, 51);
    fprintf('%-14s gap max %.3e, rms %.3e, normal angle max %.3f deg, rms %.3f deg\n', ...
        label, r.maxGap, r.rmsGap, r.maxNormalAngleDeg, r.rmsNormalAngleDeg);
end


function printG1Report(label, Sa, edgeA, Sb, edgeB)
    r = geom.SurfaceTools.edgeG1Report(Sa, edgeA, Sb, edgeB, 101);
    fprintf(['%-14s gap %.3e, normal max %.3f deg, ', ...
             'cross-tangent parallel max %.3f deg\n'], ...
        label, r.maxGap, r.maxNormalAngleDeg, ...
        r.maxCrossTangentParallelAngleDeg);
end


function plotParents(S1,S2,n)
    S1.plot(n,n,'ShowCP',false,'FaceColor',[0.25 0.55 0.95],'Alpha',0.50);
    S2.plot(n,n,'ShowCP',false,'FaceColor',[0.25 0.95 0.55],'Alpha',0.50);
end


function plotCase(nr,nc,idx,name,S1,S2,F,ns,nf)
    subplot(nr,nc,idx);
    hold on; grid on; axis equal;
    title(name);

    plotParents(S1,S2,ns);
    F.plot(ns,nf,'ShowCP',true,'FaceColor',[0.95 0.65 0.25],'Alpha',0.85);

    xlabel('X'); ylabel('Y'); zlabel('Z');
    view(35,22);
end



function plotSingleCase(name,S1,S2,F,ns,nf)

    hold on;
    grid on;
    axis equal;

    title(name);

    plotParents(S1,S2,ns);

    F.plot(ns,nf, ...
        'ShowCP',true, ...
        'FaceColor',[0.95 0.65 0.25], ...
        'Alpha',0.85);

    xlabel('X');
    ylabel('Y');
    zlabel('Z');

    view(35,22);

end


function [Sb, St, Sl, Sr] = buildFourEdgeParents()
%BUILDFOUREDGE PARENTS Build four surfaces around a rectangular hole.
%
% Hole region is approximately:
%   x in [1,2]
%   y in [1,2]
%
% Bottom parent: y <= 1
% Top parent:    y >= 2
% Left parent:   x <= 1
% Right parent:  x >= 2

    U = [0 0 0 0 1 1 1 1];
    V = [0 0 0 0 1 1 1 1];
    p = 3;
    q = 3;

    zfun = @(x,y) 0.18*sin(pi*x/3).*cos(pi*y/3) + ...
                  0.05*x.*y/3 + ...
                  0.06*sin(0.7*x + 0.4*y);

    % Bottom patch: x 1->2, y 0->1, fill edge is v1.
    Pb = zeros(4,4,3);
    for i = 1:4
        u = (i-1)/3;
        x = 1 + u;
        for j = 1:4
            v = (j-1)/3;
            y = v;
            Pb(i,j,:) = [x, y, zfun(x,y)];
        end
    end

    % Top patch: x 1->2, y 2->3, fill edge is v0.
    Pt = zeros(4,4,3);
    for i = 1:4
        u = (i-1)/3;
        x = 1 + u;
        for j = 1:4
            v = (j-1)/3;
            y = 2 + v;
            Pt(i,j,:) = [x, y, zfun(x,y)];
        end
    end

    % Left patch: x 0->1, y 1->2, fill edge is u1.
    Pl = zeros(4,4,3);
    for i = 1:4
        u = (i-1)/3;
        x = u;
        for j = 1:4
            v = (j-1)/3;
            y = 1 + v;
            Pl(i,j,:) = [x, y, zfun(x,y)];
        end
    end

    % Right patch: x 2->3, y 1->2, fill edge is u0.
    Pr = zeros(4,4,3);
    for i = 1:4
        u = (i-1)/3;
        x = 2 + u;
        for j = 1:4
            v = (j-1)/3;
            y = 1 + v;
            Pr(i,j,:) = [x, y, zfun(x,y)];
        end
    end

    Sb = geom.NURBSSurface(Pb, p, q, U, V);
    St = geom.NURBSSurface(Pt, p, q, U, V);
    Sl = geom.NURBSSurface(Pl, p, q, U, V);
    Sr = geom.NURBSSurface(Pr, p, q, U, V);
end


function plotFourEdgeCase(name, Sb, St, Sl, Sr, F)
    hold on;
    grid on;
    axis equal;
    title(name);
    xlabel('X'); ylabel('Y'); zlabel('Z');

    Sb.plot(18,18,'ShowCP',false,'FaceColor',[0.25 0.55 0.95],'Alpha',0.45);
    St.plot(18,18,'ShowCP',false,'FaceColor',[0.25 0.95 0.55],'Alpha',0.45);
    Sl.plot(18,18,'ShowCP',false,'FaceColor',[0.75 0.55 0.95],'Alpha',0.45);
    Sr.plot(18,18,'ShowCP',false,'FaceColor',[0.55 0.95 0.95],'Alpha',0.45);

    F.plot(24,24,'ShowCP',true,'FaceColor',[0.95 0.65 0.25],'Alpha',0.90);

    view(35,28);
end

function printInteriorG1Report(label, Sa, edgeA, Sb, edgeB)
    r = geom.SurfaceTools.edgeG1ReportInterior(Sa, edgeA, Sb, edgeB, ...
        'NSample',101, ...
        'TrimFraction',0.05);

    fprintf(['%-16s gap %.3e, normal max %.3f deg, rms %.3f deg, ', ...
             'cross-tangent max %.3f deg\n'], ...
        label, r.maxGap, r.maxNormalAngleDeg, r.rmsNormalAngleDeg, ...
        r.maxCrossTangentParallelAngleDeg);
end