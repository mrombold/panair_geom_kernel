%% DEMO_NURBS_SURFACE_KERNEL.m
% Comprehensive demo / validation for the NURBS surface kernel.
%
% Sections:
%   1.  Base surface construction / validation
%   2.  Evaluation / derivatives / normal / curvature
%   3.  Iso-mesh / normals / closest point
%   4.  Iso-curves
%   5.  Refinement
%   6.  SplitU / SplitV
%   7.  Subpatch extraction
%   8.  Bezier patch decomposition
%   9.  Degree elevation / reduction
%  10.  Knot removal
%  11.  Rectangular interpolation
%  12.  Rectangular LSQ fit
%  13.  Ruled surface
%  14.  Bilinear patch
%  15.  Loft / skinning
%  16.  Coons patch
%  17.  Gordon surface
%  18.  Multi-Gordon blend
%  19.  Surface of revolution
%  20.  Swept surface
%  21.  Trimmed surface: UV loops / inside-outside logic
%  22.  Trimmed evaluation / trim mask
%  23.  Trim-preserving edits smoke test

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
%  1. BASE SURFACE CONSTRUCTION / VALIDATION
%% ======================================================================
fprintf('\n--- 1. Base surface construction / validation ---\n');

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

fig1 = figure('Name','1 - Base surface');
S.plot(24, 24, 'ShowCP', true, 'ShowIso', true);
title('Base surface');
view(3);

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
%  3. ISO-MESH / NORMALS / CLOSEST POINT
%% ======================================================================
fprintf('\n--- 3. Iso-mesh / normals / closest point ---\n');

try
    mesh = S.isoMesh(20, 20);
    fprintf('  mesh size = %d x %d\n', mesh.nu, mesh.nv);
    fprintf('  # quads   = %d\n', size(mesh.connectivity,1));

    fig2 = figure('Name','2 - Surface normals');
    S.plot(18, 18, 'ShowCP', false, 'ShowIso', false);
    hold on;
    S.plotNormals(8, 8, 0.18);
    title('Surface normals');
    view(3);

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
    fprintf('  Mesh / closest-point section failed: %s\n', ME.message);
end

%% ======================================================================
%  4. ISO-CURVES
%% ======================================================================
fprintf('\n--- 4. Iso-curves ---\n');

try
    Cu = S.isoCurveU(0.35);
    Cv = S.isoCurveV(0.65);

    fig4 = figure('Name','4 - Iso-curves');
    S.plot(18, 18, 'ShowCP', false, 'ShowIso', false);
    hold on;
    Cu.plot(120, 'ShowCP', true, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
    Cv.plot(120, 'ShowCP', true, 'Color', [0.1 0.1 0.9], 'LineWidth', 2.0);
    title('Iso-curves');
    view(3);
catch ME
    fprintf('  Iso-curve section failed: %s\n', ME.message);
end

%% ======================================================================
%  5. REFINEMENT
%% ======================================================================
fprintf('\n--- 5. Refinement ---\n');

try
    Sref = S.refine([0.25 0.50 0.75], [0.33 0.66]);

    fprintf('  original net size = %d x %d\n', size(S.P,1), size(S.P,2));
    fprintf('  refined  net size = %d x %d\n', size(Sref.P,1), size(Sref.P,2));

    fig5 = figure('Name','5 - Refinement');
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
%  6. SPLITU / SPLITV
%% ======================================================================
fprintf('\n--- 6. SplitU / SplitV ---\n');

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
%  7. SUBPATCH EXTRACTION
%% ======================================================================
fprintf('\n--- 7. Subpatch extraction ---\n');

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
%  8. BEZIER PATCH DECOMPOSITION
%% ======================================================================
fprintf('\n--- 8. Bezier patch decomposition ---\n');

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
%  9. DEGREE ELEVATION / REDUCTION
%% ======================================================================
fprintf('\n--- 9. Degree elevation / reduction ---\n');

try
    SeU = S.elevateU(1);
    SeV = S.elevateV(1);

    fprintf('  elevateU: p %d -> %d\n', S.p, SeU.p);
    fprintf('  elevateV: q %d -> %d\n', S.q, SeV.q);

    [SrU, errU] = SeU.reduceDegreeU(1, 1e-8);
    [SrV, errV] = SeV.reduceDegreeV(1, 1e-8);

    fprintf('  reduceDegreeU maxErr = %.3e\n', errU);
    fprintf('  reduceDegreeV maxErr = %.3e\n', errV);

    fig10 = figure('Name','10 - Degree elevation');
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
% 10. KNOT REMOVAL
%% ======================================================================
fprintf('\n--- 10. Knot removal ---\n');

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
% 11. RECTANGULAR INTERPOLATION
%% ======================================================================
fprintf('\n--- 11. Rectangular interpolation ---\n');

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
    fprintf('  interpolated net size = %d x %d\n', size(Sint.P,1), size(Sint.P,2));

    fig13 = figure('Name','13 - Global interpolation');
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
% 12. RECTANGULAR LSQ FIT
%% ======================================================================
fprintf('\n--- 12. Rectangular LSQ fit ---\n');

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
    fprintf('  LSQ fit net size = %d x %d\n', size(Sfit.P,1), size(Sfit.P,2));

    fig14 = figure('Name','14 - Global LSQ fit');
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
% 13. RULED SURFACE
%% ======================================================================
fprintf('\n--- 13. Ruled surface ---\n');

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
    fprintf('  ruled surface net size = %d x %d\n', size(Srule.P,1), size(Srule.P,2));

    fig15 = figure('Name','15 - Ruled surface');
    Srule.plot(24, 12, 'ShowCP', true, 'ShowIso', true);
    title('Ruled surface');
    view(3);
catch ME
    fprintf('  Ruled surface failed: %s\n', ME.message);
end

%% ======================================================================
% 14. BILINEAR PATCH
%% ======================================================================
fprintf('\n--- 14. Bilinear patch ---\n');

try
    Sbil = geom.NURBSSurface.bilinearCoons( ...
        [0 0 0], [2 0 0.5], [0 2 0.2], [2 2 0.8]);

    fprintf('  bilinear patch p,q = %d, %d\n', Sbil.p, Sbil.q);

    fig16 = figure('Name','16 - Bilinear patch');
    Sbil.plot(8, 8, 'ShowCP', true, 'ShowIso', true);
    title('Bilinear patch');
    view(3);
catch ME
    fprintf('  Bilinear patch failed: %s\n', ME.message);
end

%% ======================================================================
% 15. LOFT / SKINNING
%% ======================================================================
fprintf('\n--- 15. Loft / skinning ---\n');

try
    c1 = geom.NURBSCurve([0 0 0; 1 0.4 0.1; 2 0.1 0.0; 3 0 0], 3);
    c2 = geom.NURBSCurve([0 0.5 0.6; 1 0.9 0.8; 2 0.8 0.5; 3 0.6 0.3], 3);
    c3 = geom.NURBSCurve([0 1.1 1.0; 1 1.4 1.2; 2 1.5 0.7; 3 1.4 0.4], 3);
    c4 = geom.NURBSCurve([0 1.8 1.2; 1 2.0 1.0; 2 2.2 0.5; 3 2.3 0.2], 3);

    [cc, pmax, Ucommon] = geom.NURBSSurface.makeCompatibleCurves({c1,c2,c3,c4});
    fprintf('  compatible curves: pmax = %d, knot count = %d\n', pmax, numel(Ucommon)); %#ok<NASGU>

    Sloft = geom.NURBSSurface.loft({c1,c2,c3,c4}, 3, 'centripetal');
    fprintf('  loft net size = %d x %d\n', size(Sloft.P,1), size(Sloft.P,2));

    fig17 = figure('Name','17 - Lofted surface');
    Sloft.plot(24, 24, 'ShowCP', true, 'ShowIso', true);
    hold on;
    c1.plot(120, 'ShowCP', false, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
    c2.plot(120, 'ShowCP', false, 'Color', [0.1 0.5 0.9], 'LineWidth', 2.0);
    c3.plot(120, 'ShowCP', false, 'Color', [0.1 0.7 0.2], 'LineWidth', 2.0);
    c4.plot(120, 'ShowCP', false, 'Color', [0.8 0.4 0.1], 'LineWidth', 2.0);
    title('Loft / skinning');
    view(3);
catch ME
    fprintf('  Loft section failed: %s\n', ME.message);
end

%% ======================================================================
% 16. COONS PATCH
%% ======================================================================
fprintf('\n--- 16. Coons patch ---\n');

try
    Cu0 = geom.NURBSCurve([0 0 0; 1 0.0 0.3; 2 0.0 0.2; 3 0 0], 3);
    Cu1 = geom.NURBSCurve([0 2 0.2; 1 2.0 0.8; 2 2.0 0.7; 3 2 0.1], 3);
    Cv0 = geom.NURBSCurve([0 0 0; 0 0.8 0.4; 0 1.4 0.5; 0 2 0.2], 3);
    Cv1 = geom.NURBSCurve([3 0 0; 3 0.7 0.2; 3 1.4 0.4; 3 2 0.1], 3);

    Scoons = geom.NURBSSurface.coons(Cu0, Cu1, Cv0, Cv1, 3, 3, 21, 21);
    fprintf('  Coons net size = %d x %d\n', size(Scoons.P,1), size(Scoons.P,2));

    fig18 = figure('Name','18 - Coons patch');
    Scoons.plot(24, 24, 'ShowCP', true, 'ShowIso', true);
    hold on;
    Cu0.plot(100, 'ShowCP', false, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
    Cu1.plot(100, 'ShowCP', false, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
    Cv0.plot(100, 'ShowCP', false, 'Color', [0.1 0.1 0.9], 'LineWidth', 2.0);
    Cv1.plot(100, 'ShowCP', false, 'Color', [0.1 0.1 0.9], 'LineWidth', 2.0);
    title('Coons patch');
    view(3);
catch ME
    fprintf('  Coons section failed: %s\n', ME.message);
end

%% ======================================================================
% 17. GORDON SURFACE / NETWORK COMPATIBILITY
%% ======================================================================
fprintf('\n--- 17. Gordon surface / network compatibility ---\n');

try
    p1 = geom.NURBSCurve([0 0 0; 1 0.0 0.4; 2 0.0 0.3; 3 0 0], 3);
    p2 = geom.NURBSCurve([0 0.8 0.3; 1 0.8 0.8; 2 0.8 0.6; 3 0.8 0.2], 3);
    p3 = geom.NURBSCurve([0 1.6 0.1; 1 1.6 0.7; 2 1.6 0.9; 3 1.6 0.3], 3);
    p4 = geom.NURBSCurve([0 2.4 0.0; 1 2.4 0.2; 2 2.4 0.5; 3 2.4 0.1], 3);

    g1 = geom.NURBSCurve([0 0 0; 0 0.8 0.4; 0 1.6 0.1; 0 2.4 0], 3);
    g2 = geom.NURBSCurve([1 0.0 0.4; 1 0.8 0.8; 1 1.6 0.7; 1 2.4 0.2], 3);
    g3 = geom.NURBSCurve([2 0.0 0.3; 2 0.8 0.6; 2 1.6 0.9; 2 2.4 0.5], 3);
    g4 = geom.NURBSCurve([3 0.0 0; 3 0.8 0.2; 3 1.6 0.3; 3 2.4 0.1], 3);

    [profilesC, guidesC, uPar, vPar] = geom.NURBSSurface.makeCompatibleNetwork( ...
        {p1,p2,p3,p4}, {g1,g2,g3,g4});
    fprintf('  network compatible: %d profiles, %d guides\n', numel(profilesC), numel(guidesC));
    fprintf('  uPar = %s\n', mat2str(uPar.',3));
    fprintf('  vPar = %s\n', mat2str(vPar.',3));

    Sgordon = geom.NURBSSurface.gordon({p1,p2,p3,p4}, {g1,g2,g3,g4}, 3, 3, 25, 25);
    fprintf('  Gordon net size = %d x %d\n', size(Sgordon.P,1), size(Sgordon.P,2));

    fig19 = figure('Name','19 - Gordon surface');
    Sgordon.plot(28, 28, 'ShowCP', true, 'ShowIso', true);
    hold on;
    p1.plot(100, 'ShowCP', false, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
    p2.plot(100, 'ShowCP', false, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
    p3.plot(100, 'ShowCP', false, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
    p4.plot(100, 'ShowCP', false, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
    g1.plot(100, 'ShowCP', false, 'Color', [0.1 0.1 0.9], 'LineWidth', 2.0);
    g2.plot(100, 'ShowCP', false, 'Color', [0.1 0.1 0.9], 'LineWidth', 2.0);
    g3.plot(100, 'ShowCP', false, 'Color', [0.1 0.1 0.9], 'LineWidth', 2.0);
    g4.plot(100, 'ShowCP', false, 'Color', [0.1 0.1 0.9], 'LineWidth', 2.0);
    title('Gordon surface');
    view(3);
catch ME
    fprintf('  Gordon section failed: %s\n', ME.message);
end

%% ======================================================================
% 18. MULTI-GORDON BLEND
%% ======================================================================
fprintf('\n--- 18. Multi-Gordon blend ---\n');

try
    p1b = geom.NURBSCurve([0 0 0; 1 0.1 0.2; 2 0.0 0.4; 3 0 0.1], 3);
    p2b = geom.NURBSCurve([0 0.8 0.1; 1 0.8 0.5; 2 0.8 0.7; 3 0.8 0.3], 3);
    p3b = geom.NURBSCurve([0 1.6 0.2; 1 1.6 0.4; 2 1.6 0.8; 3 1.6 0.2], 3);
    p4b = geom.NURBSCurve([0 2.4 0.1; 1 2.4 0.2; 2 2.4 0.4; 3 2.4 0.0], 3);

    g1b = geom.NURBSCurve([0 0 0; 0 0.8 0.1; 0 1.6 0.2; 0 2.4 0.1], 3);
    g2b = geom.NURBSCurve([1 0.1 0.2; 1 0.8 0.5; 1 1.6 0.4; 1 2.4 0.2], 3);
    g3b = geom.NURBSCurve([2 0.0 0.4; 2 0.8 0.7; 2 1.6 0.8; 2 2.4 0.4], 3);
    g4b = geom.NURBSCurve([3 0.1 0.1; 3 0.8 0.3; 3 1.6 0.2; 3 2.4 0.0], 3);

    Smg = geom.NURBSSurface.multiGordon( ...
        {{p1,p2,p3,p4}, {p1b,p2b,p3b,p4b}}, ...
        {{g1,g2,g3,g4}, {g1b,g2b,g3b,g4b}}, ...
        [0.65; 0.35], 3, 3, 25, 25);

    fprintf('  multi-Gordon net size = %d x %d\n', size(Smg.P,1), size(Smg.P,2));

    fig20 = figure('Name','20 - Multi-Gordon');
    Smg.plot(28, 28, 'ShowCP', true, 'ShowIso', true);
    title('Multi-Gordon blended surface');
    view(3);
catch ME
    fprintf('  Multi-Gordon section failed: %s\n', ME.message);
end

%% ======================================================================
% 19. SURFACE OF REVOLUTION
%% ======================================================================
fprintf('\n--- 19. Surface of revolution ---\n');

try
    Cprof = geom.NURBSCurve([1.0 0.0 0.0;
                             1.2 0.0 0.4;
                             1.4 0.0 0.8;
                             1.6 0.0 1.2], 3);

    Srev = geom.NURBSSurface.revolve(Cprof, [0 0 0], [0 0 1], 2*pi, 3, 9);
    fprintf('  revolution net size = %d x %d\n', size(Srev.P,1), size(Srev.P,2));

    fig21 = figure('Name','21 - Surface of revolution');
    Srev.plot(30, 28, 'ShowCP', true, 'ShowIso', true);
    title('Surface of revolution');
    view(3);
catch ME
    fprintf('  Revolution section failed: %s\n', ME.message);
end

%% ======================================================================
% 20. SWEPT SURFACE
%% ======================================================================
fprintf('\n--- 20. Swept surface ---\n');

try
    Cprof2 = geom.NURBSCurve([0  -0.4  0.0;
                              0  -0.1  0.2;
                              0   0.2  0.2;
                              0   0.5  0.0], 3);

    Cspine = geom.NURBSCurve([0.0 0.0 0.0;
                              1.0 0.4 0.5;
                              2.2 0.2 1.1;
                              3.2 0.8 1.6], 3);

    Sswp = geom.NURBSSurface.sweep(Cprof2, Cspine, 3, 10, [0 0 1]);
    fprintf('  sweep net size = %d x %d\n', size(Sswp.P,1), size(Sswp.P,2));

    fig22 = figure('Name','22 - Swept surface');
    Sswp.plot(28, 22, 'ShowCP', true, 'ShowIso', true);
    hold on;
    Cspine.plot(200, 'ShowCP', false, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.0);
    title('Swept surface');
    view(3);
catch ME
    fprintf('  Sweep section failed: %s\n', ME.message);
end

%% ======================================================================
% 21. TRIMMED SURFACE: UV LOOPS / INSIDE-OUTSIDE LOGIC
%% ======================================================================
fprintf('\n--- 21. Trimmed surface: UV loops / inside-outside logic ---\n');

try
    th = linspace(0, 2*pi, 180).';
    outer = [0.5 + 0.42*cos(th), 0.5 + 0.33*sin(th)];
    inner = [0.55 + 0.12*cos(th), 0.48 + 0.09*sin(th)];

    Strim = S.setTrims({outer}, {inner});
    fprintf('  isTrimmed() = %d\n', Strim.isTrimmed());
    fprintf('  inside trim at (0.50,0.50) = %d\n', Strim.isInsideTrim(0.50,0.50));
    fprintf('  inside trim at (0.10,0.10) = %d\n', Strim.isInsideTrim(0.10,0.10));
    fprintf('  inside trim at (0.90,0.50) = %d\n', Strim.isInsideTrim(0.90,0.50));

    fig23 = figure('Name','23 - Trim loops in UV');
    Strim.plotTrimUV();
    view(2);

    fig24 = figure('Name','24 - Trimmed surface');
    Strim.plot(32, 32, 'ShowCP', false, 'ShowIso', false, 'ShowTrims', true);
    title('Trimmed surface with outer loop and inner hole');
    view(3);
catch ME
    fprintf('  Trim logic section failed: %s\n', ME.message);
end

%% ======================================================================
% 22. TRIMMED EVALUATION / TRIM MASK
%% ======================================================================
fprintf('\n--- 22. Trimmed evaluation / trim mask ---\n');

try
    p_in = Strim.evaluateTrimmed(0.20, 0.50);
    fprintf('  evaluateTrimmed(0.20,0.50) = [%.6f %.6f %.6f]\n', p_in);

    try
        Strim.evaluateTrimmed(0.50, 0.48); % likely inside hole
        fprintf('  hole test unexpectedly evaluated.\n');
    catch MEhole
        fprintf('  hole test blocked as expected: %s\n', MEhole.message);
    end

    mg = Strim.isoMesh(25, 25, 'RespectTrim', true);
    kept = sum(mg.trimMask(:));
    total = numel(mg.trimMask);
    fprintf('  trim mask kept %d / %d mesh nodes\n', kept, total);

    fig25 = figure('Name','25 - Trim mask visualization');
    imagesc(mg.trimMask.');
    axis equal tight;
    title('Trim mask on mesh nodes');
    xlabel('u-index');
    ylabel('v-index');
    colorbar;
catch ME
    fprintf('  Trim evaluation section failed: %s\n', ME.message);
end

%% ======================================================================
% 23. TRIM-PRESERVING EDITS SMOKE TEST
%% ======================================================================
fprintf('\n--- 23. Trim-preserving edits smoke test ---\n');

try
    Srt = Strim.refine([0.3 0.6], [0.4 0.7]);
    Seu = Strim.elevateU(1);
    Sev = Strim.elevateV(1);

    fprintf('  trimmed refined net = %d x %d\n', size(Srt.P,1), size(Srt.P,2));
    fprintf('  trimmed elevateU: p %d -> %d\n', Strim.p, Seu.p);
    fprintf('  trimmed elevateV: q %d -> %d\n', Strim.q, Sev.q);
    fprintf('  refined still trimmed = %d\n', Srt.isTrimmed());
    fprintf('  elevateU still trimmed = %d\n', Seu.isTrimmed());
    fprintf('  elevateV still trimmed = %d\n', Sev.isTrimmed());

    fig26 = figure('Name','26 - Trim-preserving refine');
    Srt.plot(34, 34, 'ShowCP', false, 'ShowIso', false, 'ShowTrims', true);
    title('Trimmed + refined');
    view(3);
catch ME
    fprintf('  Trim-preserving edit section failed: %s\n', ME.message);
end

%% ======================================================================
% SUMMARY
%% ======================================================================
fprintf('\n=== Surface Kernel Demo Complete ===\n');
fprintf('Figures:\n');
fprintf('  (1)  Base surface\n');
fprintf('  (2)  Surface normals\n');
fprintf('  (3)  Closest point\n');
fprintf('  (4)  Iso-curves\n');
fprintf('  (5)  Refinement\n');
fprintf('  (6)  Split U\n');
fprintf('  (7)  Split V\n');
fprintf('  (8)  Extracted subpatch\n');
fprintf('  (9)  Bezier decomposition\n');
fprintf('  (10) Degree elevation\n');
fprintf('  (11) Degree reduction\n');
fprintf('  (12) Knot removal\n');
fprintf('  (13) Global interpolation\n');
fprintf('  (14) Global LSQ fit\n');
fprintf('  (15) Ruled surface\n');
fprintf('  (16) Bilinear patch\n');
fprintf('  (17) Loft / skinning\n');
fprintf('  (18) Coons patch\n');
fprintf('  (19) Gordon surface\n');
fprintf('  (20) Multi-Gordon\n');
fprintf('  (21) Surface of revolution\n');
fprintf('  (22) Swept surface\n');
fprintf('  (23) Trim loops in UV\n');
fprintf('  (24) Trimmed surface\n');
fprintf('  (25) Trim mask visualization\n');
fprintf('  (26) Trim-preserving refine\n');