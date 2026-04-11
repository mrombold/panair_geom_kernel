%% DEMO_NURBS_SURFACE_KERNEL.m
% Demonstration / validation script for the NURBS surface kernel.
%
% Exercises:
%   1. Construction / validation
%   2. Evaluation / derivatives / normals / curvatures
%   3. Iso-grid / mesh / plotting
%   4. Closest-point projection
%   5. Iso-curves
%   6. Refinement
%   7. SplitU / SplitV
%   8. Subpatch extraction
%   9. Bezier patch decomposition
%  10. Degree elevation / reduction
%  11. Knot removal
%  12. Global interpolation from rectangular net
%  13. Global least-squares fit from rectangular net
%  14. Ruled surface constructor
%  15. Bilinear patch constructor

clear; close all; clc;

%% ---- PATH SETUP -------------------------------------------------------
script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir)
    addpath(script_dir);
end

if ~exist('+geom', 'dir') && isempty(which('geom.NURBSSurface'))
    error(['Cannot find +geom package.\n' ...
           'Expected it at: %s\n' ...
           'Add the geom kernel folder to your MATLAB path:\n' ...
           '  addpath(''%s'')'], ...
           fullfile(script_dir, '+geom'), script_dir);
end
fprintf('geom package found at: %s\n', script_dir);

%% ======================================================================
%  1. BASE SURFACE CONSTRUCTION
%% ======================================================================
fprintf('\n--- 1. Construction / validation ---\n');

% Build a simple 4x4 cubic-cubic test surface
P = zeros(4,4,3);

xu = [0.0 1.0 2.2 3.2];
yv = [0.0 0.8 1.8 3.0];

for i = 1:4
    for j = 1:4
        x = xu(i);
        y = yv(j);
        z = 0.28 * sin(0.9*x) * cos(0.8*y) + 0.08 * x - 0.05 * y;
        P(i,j,:) = [x, y, z];
    end
end

S = geom.NURBSSurface(P, 3, 3);

fprintf('  p = %d, q = %d\n', S.p, S.q);
fprintf('  size(P) = [%d x %d x %d]\n', size(S.P,1), size(S.P,2), size(S.P,3));
fprintf('  validate() = %d\n', S.validate());
fprintf('  domainU = [%.3f %.3f]\n', S.domainU(1), S.domainU(2));
fprintf('  domainV = [%.3f %.3f]\n', S.domainV(1), S.domainV(2));

%% ======================================================================
%  2. EVALUATION / DERIVATIVES / NORMAL / CURVATURE
%% ======================================================================
fprintf('\n--- 2. Evaluation / derivatives / normal / curvature ---\n');

u0 = 0.37;
v0 = 0.61;

try
    Spt = S.evaluate(u0, v0);
    [Su, Sv] = S.partialDerivatives(u0, v0);
    [Suu, Suv, Svv] = S.secondPartials(u0, v0);
    N = S.normal(u0, v0);
    [K, H, ~] = S.curvatures(u0, v0);

    fprintf('  S(%.3f, %.3f)   = [%.6f %.6f %.6f]\n', u0, v0, Spt);
    fprintf('  Su             = [%.6f %.6f %.6f]\n', Su);
    fprintf('  Sv             = [%.6f %.6f %.6f]\n', Sv);
    fprintf('  Suu            = [%.6f %.6f %.6f]\n', Suu);
    fprintf('  Suv            = [%.6f %.6f %.6f]\n', Suv);
    fprintf('  Svv            = [%.6f %.6f %.6f]\n', Svv);
    fprintf('  N              = [%.6f %.6f %.6f]\n', N);
    fprintf('  Gaussian K     = %.6e\n', K);
    fprintf('  Mean H         = %.6e\n', H);
catch ME
    fprintf('  Derivative section failed: %s\n', ME.message);
end

%% ======================================================================
%  3. PLOT / ISO MESH / NORMALS
%% ======================================================================
fprintf('\n--- 3. Plot / iso-mesh / normals ---\n');

try
    mesh = S.isoMesh(20, 20);
    fprintf('  mesh size = %d x %d\n', mesh.nu, mesh.nv);
    fprintf('  # quads   = %d\n', size(mesh.connectivity,1));

    fig1 = figure('Name','1 - Base NURBS surface');
    S.plot(24, 24, 'ShowCP', true, 'ShowIso', true);
    title('Base NURBS surface');
    view(3);

    fig2 = figure('Name','2 - Surface normals');
    S.plot(18, 18, 'ShowCP', false, 'ShowIso', false);
    hold on;
    S.plotNormals(8, 8, 0.18);
    title('Surface normals');
    view(3);
catch ME
    fprintf('  Plot / mesh section failed: %s\n', ME.message);
end

%% ======================================================================
%  4. CLOSEST POINT
%% ======================================================================
fprintf('\n--- 4. Closest-point projection ---\n');

try
    Q = [1.55 1.25 0.55];
    [u_cp, v_cp, Pcp, dcp] = S.closestPoint(Q);

    fprintf('  closest-point (u,v) = [%.8f %.8f]\n', u_cp, v_cp);
    fprintf('  closest-point dist  = %.8e\n', dcp);

    fig3 = figure('Name','3 - Closest point');
    S.plot(24, 24, 'ShowCP', true, 'ShowIso', true);
    hold on;
    plot3(Q(1), Q(2), Q(3), 'ks', 'MarkerFaceColor','y', 'MarkerSize', 8);
    plot3(Pcp(1), Pcp(2), Pcp(3), 'ro', 'MarkerFaceColor','r', 'MarkerSize', 8);
    plot3([Q(1), Pcp(1)], [Q(2), Pcp(2)], [Q(3), Pcp(3)], 'k--', 'LineWidth', 1.2);
    title('Closest-point projection');
    view(3);
catch ME
    fprintf('  Closest-point section failed: %s\n', ME.message);
end

%% ======================================================================
%  5. ISO-CURVES
%% ======================================================================
fprintf('\n--- 5. Iso-curves ---\n');

try
    Cu = S.isoCurveU(0.35);
    Cv = S.isoCurveV(0.65);

    fig4 = figure('Name','4 - Iso-curves');
    S.plot(18, 18, 'ShowCP', false, 'ShowIso', false);
    hold on;
    Cu.plot(120, 'ShowCP', true, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
    Cv.plot(120, 'ShowCP', true, 'Color', [0.1 0.1 0.9], 'LineWidth', 2.0);
    title('Iso-curves U(v=const) and V(u=const)');
    view(3);
catch ME
    fprintf('  Iso-curve section failed: %s\n', ME.message);
end

%% ======================================================================
%  6. REFINEMENT
%% ======================================================================
fprintf('\n--- 6. Refinement ---\n');

try
    Sref = S.refine([0.25 0.50 0.75], [0.33 0.66]);

    fprintf('  original net size = %d x %d\n', size(S.P,1), size(S.P,2));
    fprintf('  refined  net size = %d x %d\n', size(Sref.P,1), size(Sref.P,2));

    fig5 = figure('Name','5 - Surface refinement');
    subplot(1,2,1);
    S.plot(20, 20, 'ShowCP', true, 'ShowIso', false);
    title('Original');
    view(3);

    subplot(1,2,2);
    Sref.plot(20, 20, 'ShowCP', true, 'ShowIso', false);
    title('Refined');
    view(3);
catch ME
    fprintf('  Refinement section failed: %s\n', ME.message);
end

%% ======================================================================
%  7. SPLIT U / SPLIT V
%% ======================================================================
fprintf('\n--- 7. SplitU / SplitV ---\n');

try
    [Sulo, Suhi] = S.splitU(0.5);
    [Svlo, Svhi] = S.splitV(0.5);

    fprintf('  splitU: low net  = %d x %d, high net = %d x %d\n', ...
        size(Sulo.P,1), size(Sulo.P,2), size(Suhi.P,1), size(Suhi.P,2));
    fprintf('  splitV: low net  = %d x %d, high net = %d x %d\n', ...
        size(Svlo.P,1), size(Svlo.P,2), size(Svhi.P,1), size(Svhi.P,2));

    fig6 = figure('Name','6 - Split U');
    hold on;
    Sulo.plot(18, 18, 'ShowCP', true, 'ShowIso', false, 'FaceColor', [0.9 0.4 0.4]);
    Suhi.plot(18, 18, 'ShowCP', true, 'ShowIso', false, 'FaceColor', [0.4 0.4 0.9]);
    title('Split in U');
    view(3);

    fig7 = figure('Name','7 - Split V');
    hold on;
    Svlo.plot(18, 18, 'ShowCP', true, 'ShowIso', false, 'FaceColor', [0.4 0.8 0.4]);
    Svhi.plot(18, 18, 'ShowCP', true, 'ShowIso', false, 'FaceColor', [0.8 0.7 0.3]);
    title('Split in V');
    view(3);
catch ME
    fprintf('  Split section failed: %s\n', ME.message);
end

%% ======================================================================
%  8. SUBPATCH EXTRACTION
%% ======================================================================
fprintf('\n--- 8. Subpatch extraction ---\n');

try
    Sex = S.extractUV([0.2 0.8], [0.15 0.75]);
    fprintf('  extracted net size = %d x %d\n', size(Sex.P,1), size(Sex.P,2));

    fig8 = figure('Name','8 - Extracted subpatch');
    hold on;
    S.plot(18, 18, 'ShowCP', false, 'ShowIso', false, 'Alpha', 0.25, 'FaceColor', [0.7 0.7 0.7]);
    Sex.plot(18, 18, 'ShowCP', true, 'ShowIso', true, 'Alpha', 0.9, 'FaceColor', [0.2 0.7 0.9]);
    title('Extracted UV subpatch');
    view(3);
catch ME
    fprintf('  Subpatch extraction failed: %s\n', ME.message);
end

%% ======================================================================
%  9. BEZIER PATCH DECOMPOSITION
%% ======================================================================
fprintf('\n--- 9. Bezier patch decomposition ---\n');

try
    patches = Sref.decomposeBezier();
    [nuPatch, nvPatch] = size(patches);
    fprintf('  # Bezier patches = %d x %d = %d\n', nuPatch, nvPatch, numel(patches));

    fig9 = figure('Name','9 - Bezier decomposition');
    hold on;
    cmap = lines(max(numel(patches), 1));
    idx = 1;
    for i = 1:nuPatch
        for j = 1:nvPatch
            patches{i,j}.surface.plot(8, 8, 'ShowCP', true, 'ShowIso', false, ...
                'Alpha', 0.85, 'FaceColor', cmap(idx,:));
            idx = idx + 1;
        end
    end
    title('Bezier patch decomposition');
    view(3);
catch ME
    fprintf('  Bezier decomposition failed: %s\n', ME.message);
end

%% ======================================================================
% 10. DEGREE ELEVATION / REDUCTION
%% ======================================================================
fprintf('\n--- 10. Degree elevation / reduction ---\n');

try
    SeU = S.elevateU(1);
    SeV = S.elevateV(1);

    fprintf('  elevateU: p %d -> %d\n', S.p, SeU.p);
    fprintf('  elevateV: q %d -> %d\n', S.q, SeV.q);

    [SrU, errU] = SeU.reduceDegreeU(1, 1e-8);
    [SrV, errV] = SeV.reduceDegreeV(1, 1e-8);

    fprintf('  reduceDegreeU maxErr = %.3e\n', errU);
    fprintf('  reduceDegreeV maxErr = %.3e\n', errV);

    fig10 = figure('Name','10 - Degree editing');
    subplot(1,2,1);
    hold on;
    S.plot(14, 14, 'ShowCP', false, 'ShowIso', false, 'FaceColor', [0.2 0.4 0.9]);
    SeU.plot(14, 14, 'ShowCP', true, 'ShowIso', false, 'FaceColor', [0.9 0.4 0.4], 'Alpha', 0.5);
    title('Elevate U');
    view(3);

    subplot(1,2,2);
    hold on;
    S.plot(14, 14, 'ShowCP', false, 'ShowIso', false, 'FaceColor', [0.2 0.4 0.9]);
    SeV.plot(14, 14, 'ShowCP', true, 'ShowIso', false, 'FaceColor', [0.4 0.8 0.4], 'Alpha', 0.5);
    title('Elevate V');
    view(3);

    fig11 = figure('Name','11 - Degree reduction');
    subplot(1,2,1);
    hold on;
    S.plot(14, 14, 'ShowCP', false, 'ShowIso', false, 'FaceColor', [0.2 0.4 0.9]);
    SrU.plot(14, 14, 'ShowCP', true, 'ShowIso', false, 'FaceColor', [0.9 0.5 0.2], 'Alpha', 0.5);
    title('Reduce U back');
    view(3);

    subplot(1,2,2);
    hold on;
    S.plot(14, 14, 'ShowCP', false, 'ShowIso', false, 'FaceColor', [0.2 0.4 0.9]);
    SrV.plot(14, 14, 'ShowCP', true, 'ShowIso', false, 'FaceColor', [0.7 0.3 0.9], 'Alpha', 0.5);
    title('Reduce V back');
    view(3);
catch ME
    fprintf('  Degree editing failed: %s\n', ME.message);
end

%% ======================================================================
% 11. KNOT REMOVAL
%% ======================================================================
fprintf('\n--- 11. Knot removal ---\n');

try
    Srm0 = S.refine([0.35 0.35], [0.65 0.65]);
    fprintf('  before removal net size = %d x %d\n', size(Srm0.P,1), size(Srm0.P,2));

    [SrmU, remU, errRemU] = Srm0.removeKnotU(0.35, 1, 1e-8);
    [SrmV, remV, errRemV] = SrmU.removeKnotV(0.65, 1, 1e-8);

    fprintf('  removeKnotU: removed = %d, maxErr = %.3e\n', remU, errRemU);
    fprintf('  removeKnotV: removed = %d, maxErr = %.3e\n', remV, errRemV);

    fig12 = figure('Name','12 - Knot removal');
    subplot(1,2,1);
    Srm0.plot(18, 18, 'ShowCP', true, 'ShowIso', false);
    title('Before removal');
    view(3);

    subplot(1,2,2);
    SrmV.plot(18, 18, 'ShowCP', true, 'ShowIso', false);
    title('After removal');
    view(3);
catch ME
    fprintf('  Knot removal failed: %s\n', ME.message);
end

%% ======================================================================
% 12. GLOBAL INTERPOLATION OF RECTANGULAR NET
%% ======================================================================
fprintf('\n--- 12. Global interpolation of rectangular net ---\n');

try
    uu = linspace(0, 1, 7);
    vv = linspace(0, 1, 6);
    QI = zeros(numel(uu), numel(vv), 3);
    for i = 1:numel(uu)
        for j = 1:numel(vv)
            u = uu(i);
            v = vv(j);
            x = 3*u;
            y = 2.5*v;
            z = 0.25*sin(pi*u)*cos(1.2*pi*v) + 0.05*u;
            QI(i,j,:) = [x, y, z];
        end
    end

    Sint = geom.NURBSSurface.globalInterpNet(QI, 3, 3);

    fig13 = figure('Name','13 - Global interpolation surface');
    Sint.plot(26, 26, 'ShowCP', true, 'ShowIso', false);
    hold on;
    QIf = reshape(QI, [], 3);
    plot3(QIf(:,1), QIf(:,2), QIf(:,3), 'ko', 'MarkerFaceColor','y', 'MarkerSize', 5);
    title('Interpolated rectangular data net');
    view(3);
catch ME
    fprintf('  Global interpolation failed: %s\n', ME.message);
end

%% ======================================================================
% 13. GLOBAL LEAST-SQUARES FIT OF RECTANGULAR NET
%% ======================================================================
fprintf('\n--- 13. Global LSQ fit of rectangular net ---\n');

try
    uu = linspace(0, 1, 18);
    vv = linspace(0, 1, 16);
    QF = zeros(numel(uu), numel(vv), 3);
    for i = 1:numel(uu)
        for j = 1:numel(vv)
            u = uu(i);
            v = vv(j);
            x = 3.2*u;
            y = 2.7*v;
            z = 0.22*sin(1.5*pi*u)*cos(1.1*pi*v) + 0.04*cos(3*pi*v);
            QF(i,j,:) = [x, y, z];
        end
    end

    Sfit = geom.NURBSSurface.globalLeastSquaresFitNet(QF, 3, 3, 8, 7);

    fig14 = figure('Name','14 - Global LSQ fit surface');
    Sfit.plot(28, 28, 'ShowCP', true, 'ShowIso', false);
    hold on;
    QFf = reshape(QF, [], 3);
    plot3(QFf(:,1), QFf(:,2), QFf(:,3), 'k.', 'MarkerSize', 8);
    title('Least-squares fit to rectangular data net');
    view(3);
catch ME
    fprintf('  Global LSQ fit failed: %s\n', ME.message);
end

%% ======================================================================
% 14. RULED SURFACE
%% ======================================================================
fprintf('\n--- 14. Ruled surface constructor ---\n');

try
    cP1 = [0.0 0.0 0.0;
           1.0 0.2 0.3;
           2.0 0.0 0.2;
           3.0 0.0 0.0];

    cP2 = [0.0 2.0 0.5;
           1.0 2.2 0.7;
           2.0 2.0 0.4;
           3.0 2.0 0.2];

    C1 = geom.NURBSCurve(cP1, 3);
    C2 = geom.NURBSCurve(cP2, 3);

    Srule = geom.NURBSSurface.ruled(C1, C2);

    fig15 = figure('Name','15 - Ruled surface');
    Srule.plot(24, 12, 'ShowCP', true, 'ShowIso', true);
    title('Ruled surface');
    view(3);
catch ME
    fprintf('  Ruled surface failed: %s\n', ME.message);
end

%% ======================================================================
% 15. BILINEAR PATCH
%% ======================================================================
fprintf('\n--- 15. Bilinear patch constructor ---\n');

try
    Sbil = geom.NURBSSurface.bilinearCoons( ...
        [0 0 0], [2 0 0.5], [0 2 0.2], [2 2 0.8]);

    fig16 = figure('Name','16 - Bilinear patch');
    Sbil.plot(8, 8, 'ShowCP', true, 'ShowIso', true);
    title('Bilinear patch');
    view(3);
catch ME
    fprintf('  Bilinear patch failed: %s\n', ME.message);
end

%% ======================================================================
% SUMMARY
%% ======================================================================
fprintf('\n=== Surface Kernel Demo Complete ===\n');
fprintf('Figures:\n');
fprintf('  (1) Base surface\n');
fprintf('  (2) Surface normals\n');
fprintf('  (3) Closest point\n');
fprintf('  (4) Iso-curves\n');
fprintf('  (5) Refinement\n');
fprintf('  (6) Split U\n');
fprintf('  (7) Split V\n');
fprintf('  (8) Extracted subpatch\n');
fprintf('  (9) Bezier decomposition\n');
fprintf(' (10) Degree elevation\n');
fprintf(' (11) Degree reduction\n');
fprintf(' (12) Knot removal\n');
fprintf(' (13) Global interpolation\n');
fprintf(' (14) Global LSQ fit\n');
fprintf(' (15) Ruled surface\n');
fprintf(' (16) Bilinear patch\n');