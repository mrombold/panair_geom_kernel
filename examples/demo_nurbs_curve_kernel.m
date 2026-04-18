%% DEMO_NURBS_CURVE_KERNEL.m
% Demonstration / validation script for the NURBS curve kernel.
%
% Exercises:
%   1. Construction, evaluation, derivatives, tangent, curvature
%   2. Arc length and arc-length-based parameter sampling
%   3. Knot refinement
%   4. Reverse / transform / translate
%   5. Closest-point projection and inversion
%   6. Diagnostics: validation, closure, continuity, Greville points
%   7. Exact split and Bezier decomposition
%   8. Global interpolation
%   9. Global least-squares fit
%  10. Knot removal
%  11. Degree elevation / degree reduction
%  12. Simple constructors
%  13. Parameterization method comparison
%  14. Weighted least-squares fit
%  15. Constrained least-squares fit

clear; close all; clc;

%% ---- PATH SETUP -------------------------------------------------------
script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir)
    addpath(script_dir);
end

if ~exist('+geom', 'dir') && isempty(which('geom.NURBSCurve'))
    error(['Cannot find +geom package.\n' ...
           'Expected it at: %s\n' ...
           'Add the geom kernel folder to your MATLAB path:\n' ...
           '  addpath(''%s'')'], ...
           fullfile(script_dir, '+geom'), script_dir);
end
fprintf('geom package found at: %s\n', script_dir);

%% ======================================================================
%  1. CONSTRUCTION, EVALUATION, DERIVATIVES
%% ======================================================================
fprintf('\n--- 1. Construction / evaluation / derivatives ---\n');

P = [0.0  0.0   0.0;
     0.8  1.2   0.2;
     1.8  1.6  -0.1;
     2.8  0.6   0.4;
     4.0  0.0   0.0];
p = 3;

C = geom.NURBSCurve(P, p);

u_test = linspace(C.domain(1), C.domain(2), 6);
pts = C.evaluate(u_test);

fprintf('  Degree: %d\n', C.p);
fprintf('  # control points: %d\n', size(C.P,1));
fprintf('  Evaluate at u = %s\n', mat2str(u_test, 3));
disp(pts);

try
    D1 = C.derivative(0.50, 1);
    D2 = C.derivative(0.50, 2);
    T  = C.tangent(0.50);
    K  = C.curvature(linspace(C.domain(1), C.domain(2), 100));

    fprintf('  D1(0.5) = [%.6f %.6f %.6f]\n', D1);
    fprintf('  D2(0.5) = [%.6f %.6f %.6f]\n', D2);
    fprintf('  T (0.5) = [%.6f %.6f %.6f]\n', T);
    fprintf('  max curvature on sample = %.6f\n', max(K));
catch ME
    fprintf('  Derivative/curvature section failed: %s\n', ME.message);
end

fig1 = figure('Name','1 - Basic NURBS curve');
C.plot(200, 'ShowCP', true, 'ShowKnots', true);
title('Basic cubic NURBS curve');
view(3);

%% ======================================================================
%  2. ARC LENGTH
%% ======================================================================
fprintf('\n--- 2. Arc length ---\n');
try
    L = C.arcLength();
    fprintf('  total arc length = %.8f\n', L);

    u_al = C.arcLengthParam(15);
    pts_al = C.evaluate(u_al);

    fig2 = figure('Name','2 - Arc-length sampling');
    C.plot(300, 'ShowCP', true, 'ShowKnots', false);
    hold on;
    plot3(pts_al(:,1), pts_al(:,2), pts_al(:,3), 'ro', ...
          'MarkerFaceColor','r', 'MarkerSize', 5);
    title('Arc-length-based parameter sampling');
    view(3);
catch ME
    fprintf('  Arc-length section failed: %s\n', ME.message);
end

%% ======================================================================
%  3. KNOT REFINEMENT
%% ======================================================================
fprintf('\n--- 3. Knot refinement ---\n');
try
    X = [0.25 0.50 0.75];
    Cref = C.refine(X);

    fprintf('  original #CP = %d, refined #CP = %d\n', size(C.P,1), size(Cref.P,1));

    fig3 = figure('Name','3 - Knot refinement');
    subplot(1,2,1);
    C.plot(200, 'ShowCP', true, 'ShowKnots', true);
    title('Original');
    view(3);

    subplot(1,2,2);
    Cref.plot(200, 'ShowCP', true, 'ShowKnots', true);
    title('After knot insertion');
    view(3);
catch ME
    fprintf('  Knot refinement section failed: %s\n', ME.message);
end

%% ======================================================================
%  4. REVERSE / TRANSFORM / TRANSLATE
%% ======================================================================
fprintf('\n--- 4. Reverse / transform / translate ---\n');
try
    Crev = C.reverse();
    Ct   = C.translate([0 0 0.8]);

    th = deg2rad(25);
    Rz = [cos(th) -sin(th) 0 0;
          sin(th)  cos(th) 0 0;
          0        0       1 0;
          0        0       0 1];
    Cr = C.transform(Rz);

    fig4 = figure('Name','4 - Reverse / transform / translate');
    hold on; grid on; axis equal;
    C.plot(200, 'ShowCP', true,  'Color', [0.0 0.35 0.85], 'LineWidth', 1.8);
    Crev.plot(200, 'ShowCP', false, 'Color', [0.85 0.2 0.2], 'LineWidth', 1.5);
    Ct.plot(200,   'ShowCP', false, 'Color', [0.1 0.65 0.1], 'LineWidth', 1.5);
    Cr.plot(200,   'ShowCP', false, 'Color', [0.65 0.2 0.65], 'LineWidth', 1.5);
    legend({'Original','Reversed','Translated','Rotated'}, 'Location','best');
    title('Curve transforms');
    view(3);
catch ME
    fprintf('  Reverse/transform section failed: %s\n', ME.message);
end

%% ======================================================================
%  5. CLOSEST-POINT PROJECTION / INVERSION
%% ======================================================================
fprintf('\n--- 5. Closest point / inversion ---\n');
try
    Q = [2.0, 0.8, 0.7];
    [u_cp, Pcp, dcp] = C.closestPoint(Q);

    fprintf('  closest-point u = %.8f\n', u_cp);
    fprintf('  closest-point distance = %.8f\n', dcp);

    fig5 = figure('Name','5 - Closest point');
    C.plot(250, 'ShowCP', true, 'ShowKnots', false);
    hold on;
    plot3(Q(1), Q(2), Q(3), 'ks', 'MarkerFaceColor','y', 'MarkerSize', 8);
    plot3(Pcp(1), Pcp(2), Pcp(3), 'ro', 'MarkerFaceColor','r', 'MarkerSize', 8);
    plot3([Q(1), Pcp(1)], [Q(2), Pcp(2)], [Q(3), Pcp(3)], 'k--', 'LineWidth', 1.2);
    legend({'Curve','Control polygon','Query point','Closest point','Distance'}, 'Location','best');
    title('Closest point projection');
    view(3);
catch ME
    fprintf('  Closest-point section failed: %s\n', ME.message);
end

%% ======================================================================
%  6. DIAGNOSTICS
%% ======================================================================
fprintf('\n--- 6. Diagnostics ---\n');
try
    if ismethod(C, 'validate')
        fprintf('  validate() = %d\n', C.validate());
    end

    if ismethod(C, 'isClamped')
        fprintf('  isClamped() = %d\n', C.isClamped());
    end

    if ismethod(C, 'isClosed')
        fprintf('  isClosed()  = %d\n', C.isClosed(1e-9));
    end

    if ismethod(C, 'grevilleAbscissae')
        g = C.grevilleAbscissae();
        fprintf('  Greville abscissae:\n');
        disp(g(:).');
    end

    if ismethod(C, 'continuityAt')
        Uu = unique(C.U);
        Uu = Uu(Uu > C.domain(1) & Uu < C.domain(2));
        if isempty(Uu)
            fprintf('  No interior knots for continuity checks.\n');
        else
            for k = 1:numel(Uu)
                fprintf('  continuity at u=%.6f : C^%d\n', Uu(k), C.continuityAt(Uu(k)));
            end
        end
    end
catch ME
    fprintf('  Diagnostics section failed: %s\n', ME.message);
end

%% ======================================================================
%  7. EXACT SPLIT / BEZIER DECOMPOSITION
%% ======================================================================
fprintf('\n--- 7. Split / Bezier decomposition ---\n');
try
    if ismethod(C, 'split')
        [CL, CR] = C.split(0.50);

        fig6 = figure('Name','6 - Exact split');
        hold on; grid on; axis equal;
        CL.plot(200, 'ShowCP', true, 'Color', [0.1 0.4 0.9], 'LineWidth', 1.8);
        CR.plot(200, 'ShowCP', true, 'Color', [0.9 0.2 0.2], 'LineWidth', 1.8);
        title('Exact split at u = 0.5');
        legend({'Left','Right'}, 'Location','best');
        view(3);
    end

    segs = {};
    if ismethod(C, 'decomposeToBezier')
        segs = C.decomposeToBezier();
    elseif ismethod(C, 'decomposeBezier')
        parts = C.decomposeBezier();
        segs = cell(size(parts));
        for i = 1:numel(parts)
            segs{i} = parts{i}.curve;
        end
    end

    if ~isempty(segs)
        fprintf('  # Bezier segments = %d\n', numel(segs));

        fig7 = figure('Name','7 - Bezier decomposition');
        hold on; grid on; axis equal;
        cmap = lines(max(numel(segs), 1));
        for i = 1:numel(segs)
            segs{i}.plot(150, 'ShowCP', true, 'Color', cmap(i,:), 'LineWidth', 1.8);
        end
        title('Bezier decomposition');
        view(3);
    end
catch ME
    fprintf('  Split/Bezier section failed: %s\n', ME.message);
end

%% ======================================================================
%  8. GLOBAL INTERPOLATION
%% ======================================================================
fprintf('\n--- 8. Global interpolation ---\n');
try
    if ismethod('geom.NURBSCurve', 'globalInterp')
        th = linspace(0, 2*pi, 15).';
        th(end) = [];
        QI = [cos(th), 0.6*sin(th), 0.15*sin(2*th)];

        Cint = geom.NURBSCurve.globalInterp(QI, 3, 'centripetal');

        fig8 = figure('Name','8 - Global interpolation');
        Cint.plot(400, 'ShowCP', true, 'ShowKnots', false);
        hold on;
        plot3(QI(:,1), QI(:,2), QI(:,3), 'ko', 'MarkerFaceColor','y', 'MarkerSize', 6);
        title('Global interpolation through data');
        legend({'Interpolant','Control polygon','Data points'}, 'Location','best');
        axis equal;
        view(3);
    end
catch ME
    fprintf('  Global interpolation section failed: %s\n', ME.message);
end

%% ======================================================================
%  9. GLOBAL LEAST-SQUARES FIT
%% ======================================================================
fprintf('\n--- 9. Global least-squares fit ---\n');
try
    if ismethod('geom.NURBSCurve', 'globalLeastSquaresFit')
        x = linspace(0, 1, 80).';
        y = 0.25*sin(2*pi*x) + 0.03*cos(10*pi*x);
        z = 0.10*cos(2*pi*x);
        QF = [x, y, z];

        Cfit = geom.NURBSCurve.globalLeastSquaresFit(QF, 3, 12, 'chord');

        ufit = geom.NURBSCurve.parameterizeData(QF, 'chord', 3);
        Pfit = Cfit.evaluate(ufit);
        rms_fit = sqrt(mean(sum((Pfit - QF).^2, 2)));

        fprintf('  LSQ fit RMS (sampled at data parameters) = %.8e\n', rms_fit);

        fig9 = figure('Name','9 - Global LSQ fit');
        Cfit.plot(400, 'ShowCP', true, 'ShowKnots', false);
        hold on;
        plot3(QF(:,1), QF(:,2), QF(:,3), 'k.', 'MarkerSize', 10);
        title(sprintf('Global least-squares fit (RMS ~ %.3e)', rms_fit));
        legend({'Fit','Control polygon','Data'}, 'Location','best');
        axis equal;
        view(3);
    end
catch ME
    fprintf('  Global LSQ fit section failed: %s\n', ME.message);
end

%% ======================================================================
% 10. KNOT REMOVAL
%% ======================================================================
fprintf('\n--- 10. Knot removal ---\n');
try
    if ismethod(C, 'removeKnot')
        Ck = C.refine([0.35 0.35 0.65 0.65]);
        fprintf('  Before removal: #CP = %d, #knots = %d\n', size(Ck.P,1), numel(Ck.U));

        [Crm1, removed1, err1] = Ck.removeKnot(0.35, 1, 1e-8);
        fprintf('  removeKnot(0.35): removed = %d, maxErr = %.3e\n', removed1, err1);

        [Crm2, removed2, err2] = Crm1.removeKnot(0.65, 1, 1e-8);
        fprintf('  removeKnot(0.65): removed = %d, maxErr = %.3e\n', removed2, err2);

        fig10 = figure('Name','10 - Knot removal');
        hold on; grid on; axis equal;
        Ck.plot(300,   'ShowCP', true, 'Color', [0.1 0.35 0.9], 'LineWidth', 1.8);
        Crm2.plot(300, 'ShowCP', true, 'Color', [0.9 0.15 0.15], 'LineWidth', 1.5);
        title('Knot removal comparison');
        legend({'Before removal','After removal'}, 'Location','best');
        view(3);
    end
catch ME
    fprintf('  Knot removal section failed: %s\n', ME.message);
end

%% ======================================================================
% 11. DEGREE ELEVATION / DEGREE REDUCTION
%% ======================================================================
fprintf('\n--- 11. Degree elevation / degree reduction ---\n');
try
    Celev = [];
    Cred  = [];

    if ismethod(C, 'elevate')
        Celev = C.elevate(1);
        fprintf('  elevate(1): degree %d -> %d\n', C.p, Celev.p);
    end

    if ~isempty(Celev) && ismethod(Celev, 'reduceDegree')
        try
            [Cred, errRed] = Celev.reduceDegree(1, 1e-8);
            fprintf('  reduceDegree maxErr = %.3e\n', errRed);
        catch MEred
            fprintf('  reduceDegree present but failed in this run: %s\n', MEred.message);
        end
    end

    if ~isempty(Celev)
        fig11 = figure('Name','11 - Degree elevation / reduction');
        hold on; grid on; axis equal;
        C.plot(300,     'ShowCP', false, 'Color', [0.0 0.2 0.8], 'LineWidth', 2.2);
        Celev.plot(300, 'ShowCP', true,  'Color', [0.1 0.7 0.1], 'LineWidth', 1.6);
        if ~isempty(Cred)
            Cred.plot(300, 'ShowCP', true, 'Color', [0.85 0.2 0.2], 'LineWidth', 1.4);
            legend({'Original','Elevated','Reduced'}, 'Location','best');
        else
            legend({'Original','Elevated'}, 'Location','best');
        end
        title('Degree editing');
        view(3);
    end
catch ME
    fprintf('  Degree editing section failed: %s\n', ME.message);
end

%% ======================================================================
% 12. SIMPLE CONSTRUCTORS
%% ======================================================================
fprintf('\n--- 12. Simple constructors ---\n');
try
    fig12 = figure('Name','12 - Simple constructors');
    hold on; grid on; axis equal;

    legend_entries = {};

    if ismethod('geom.NURBSCurve', 'line')
        Cline = geom.NURBSCurve.line([0 0 0], [1.5 0.2 0.0]);
        Cline.plot(80, 'ShowCP', true, 'Color', [0.1 0.4 0.9], 'LineWidth', 1.8);
        legend_entries{end+1} = 'Line'; %#ok<SAGROW>
    end

    if ismethod('geom.NURBSCurve', 'quadraticConic')
        Ccon = geom.NURBSCurve.quadraticConic([0 0 0], [0.8 1.0 0], [1.8 0 0], 0.55);
        Ccon.plot(150, 'ShowCP', true, 'Color', [0.9 0.2 0.2], 'LineWidth', 1.8);
        legend_entries{end+1} = 'Conic'; %#ok<SAGROW>
    end

    title('Simple constructors: line and quadratic conic');
    if ~isempty(legend_entries)
        legend(legend_entries, 'Location','best');
    end
    view(3);
catch ME
    fprintf('  Simple constructors section failed: %s\n', ME.message);
end

%% ======================================================================
% 13. PARAMETERIZATION METHODS
%% ======================================================================
fprintf('\n--- 13. Parameterization methods ---\n');
try
    x = linspace(0, 1, 16).';
    Qp = [x, 0.30*sin(2*pi*x) + 0.06*cos(7*pi*x), 0.10*cos(2*pi*x)];

    methods_list = {'uniform','chord','centripetal','arc_length'};
    curves_param = cell(size(methods_list));

    fig13 = figure('Name','13 - Parameterization methods');
    hold on; grid on; axis equal;
    plot3(Qp(:,1), Qp(:,2), Qp(:,3), 'ko', 'MarkerFaceColor','y', 'MarkerSize', 5);

    cmap = lines(numel(methods_list));
    for k = 1:numel(methods_list)
        curves_param{k} = geom.NURBSCurve.globalInterp(Qp, 3, methods_list{k});
        curves_param{k}.plot(300, 'ShowCP', false, 'Color', cmap(k,:), 'LineWidth', 1.8);
        fprintf('  %s parameterization: built interpolant\n', methods_list{k});
    end

    title('Interpolation under different parameterizations');
    legend({'Data','Uniform','Chord','Centripetal','Arc-length'}, 'Location','best');
    view(3);
catch ME
    fprintf('  Parameterization section failed: %s\n', ME.message);
end

%% ======================================================================
% 14. WEIGHTED LEAST-SQUARES FIT
%% ======================================================================
fprintf('\n--- 14. Weighted least-squares fit ---\n');
try
    x = linspace(0, 1, 100).';
    Qw = [x, 0.24*sin(2*pi*x) + 0.05*cos(9*pi*x), 0.08*cos(2*pi*x)];

    w = ones(size(Qw,1),1);
    w(1) = 25;
    w(end) = 25;
    w(round(end/2)) = 10;

    Cuw = geom.NURBSCurve.globalLeastSquaresFit(Qw, 3, 12, 'chord');
    Cwt = geom.NURBSCurve.globalWeightedLeastSquaresFit(Qw, 3, 12, w, 'chord');

    udata = geom.NURBSCurve.parameterizeData(Qw, 'chord', 3);
    rms_u = sqrt(mean(sum((Cuw.evaluate(udata) - Qw).^2, 2)));
    rms_w = sqrt(mean(sum((Cwt.evaluate(udata) - Qw).^2, 2)));

    fprintf('  unweighted LSQ RMS = %.8e\n', rms_u);
    fprintf('  weighted   LSQ RMS = %.8e\n', rms_w);

    fig14 = figure('Name','14 - Weighted least-squares fit');
    hold on; grid on; axis equal;
    plot3(Qw(:,1), Qw(:,2), Qw(:,3), 'k.', 'MarkerSize', 9);
    Cuw.plot(400, 'ShowCP', true, 'ShowKnots', false, 'Color', [0.1 0.3 0.9], 'LineWidth', 1.8);
    Cwt.plot(400, 'ShowCP', true, 'ShowKnots', false, 'Color', [0.9 0.2 0.2], 'LineWidth', 1.8);
    title('Weighted vs unweighted LSQ fit');
    legend({'Data','Unweighted LSQ','Weighted LSQ'}, 'Location','best');
    view(3);
catch ME
    fprintf('  Weighted LSQ section failed: %s\n', ME.message);
end

%% ======================================================================
% 15. CONSTRAINED LEAST-SQUARES FIT
%% ======================================================================
fprintf('\n--- 15. Constrained least-squares fit ---\n');
try
    x = linspace(0, 1, 90).';
    Qc = [x, 0.22*sin(2*pi*x) + 0.02*cos(11*pi*x), 0.07*cos(2*pi*x)];

    Cfree = geom.NURBSCurve.globalLeastSquaresFit(Qc, 3, 12, 'centripetal');
    Cfix  = geom.NURBSCurve.globalLeastSquaresFitFixedEnds(Qc, 3, 12, 'centripetal');

    fprintf('  free fit first CP  = [%.6f %.6f %.6f]\n', Cfree.P(1,:));
    fprintf('  fixed fit first CP = [%.6f %.6f %.6f]\n', Cfix.P(1,:));
    fprintf('  data first point   = [%.6f %.6f %.6f]\n', Qc(1,:));

    fprintf('  free fit last CP   = [%.6f %.6f %.6f]\n', Cfree.P(end,:));
    fprintf('  fixed fit last CP  = [%.6f %.6f %.6f]\n', Cfix.P(end,:));
    fprintf('  data last point    = [%.6f %.6f %.6f]\n', Qc(end,:));

    fig15 = figure('Name','15 - Constrained least-squares fit');
    hold on; grid on; axis equal;
    plot3(Qc(:,1), Qc(:,2), Qc(:,3), 'k.', 'MarkerSize', 9);
    Cfree.plot(400, 'ShowCP', true, 'ShowKnots', false, 'Color', [0.1 0.4 0.9], 'LineWidth', 1.8);
    Cfix.plot(400,  'ShowCP', true, 'ShowKnots', false, 'Color', [0.9 0.2 0.2], 'LineWidth', 1.8);
    plot3(Qc(1,1), Qc(1,2), Qc(1,3), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 7);
    plot3(Qc(end,1), Qc(end,2), Qc(end,3), 'ks', 'MarkerFaceColor', 'g', 'MarkerSize', 7);
    title('Constrained LSQ with fixed endpoints');
    legend({'Data','Free LSQ','Fixed-end LSQ','Start point','End point'}, 'Location','best');
    view(3);
catch ME
    fprintf('  Constrained LSQ section failed: %s\n', ME.message);
end

%% ======================================================================
%  SUMMARY
%% ======================================================================
fprintf('\n=== Curve Kernel Demo Complete ===\n');
fprintf('Figures:\n');
fprintf('  (1) Basic curve\n');
fprintf('  (2) Arc length sampling\n');
fprintf('  (3) Knot refinement\n');
fprintf('  (4) Transformations\n');
fprintf('  (5) Closest point\n');
fprintf('  (6) Exact split\n');
fprintf('  (7) Bezier decomposition\n');
fprintf('  (8) Global interpolation\n');
fprintf('  (9) LSQ fit\n');
fprintf(' (10) Knot removal\n');
fprintf(' (11) Degree editing\n');
fprintf(' (12) Constructors\n');
fprintf(' (13) Parameterization methods\n');
fprintf(' (14) Weighted LSQ fit\n');
fprintf(' (15) Constrained LSQ fit\n');