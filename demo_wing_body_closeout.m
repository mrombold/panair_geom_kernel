function demo_wing_body_closeout()
% DEMO_WING_BODY_CLOSEOUT
% Build a fuselage-like body and a half-wing, compute the intersection,
% split both into useful local regions, and create simple upper/lower
% closeout patches between the wing-root side and the fuselage side.
%
% This demo is intended as the next step after demo_fuselage_six_patch_split.
%
% What it does:
%   1) builds fuselage and wing surfaces
%   2) computes a sampled/refined wing-body intersection
%   3) splits fuselage into forward/mid/aft and upper/lower
%   4) splits wing into root and outboard bands
%   5) builds simple ruled closeout surfaces (upper and lower)
%   6) meshes everything and optionally packages WGS networks
%
% Notes:
%   - This is a first-pass structured workflow, not a full CAD trim/B-rep.
%   - The closeouts here are ruled surfaces between extracted boundary curves.
%   - For higher quality junctions, replace ruled closeouts with Coons or
%     tangent-matched patches later.
%
% Requirements:
%   +geom/Aircraft.m
%   +geom/PatchOps.m
%   +geom/SurfaceSurfaceIntersect.m
%   +geom/PlotHelpers.m
%   +geom/MeshWriter.m   (optional)

    close all;
    clc;

    fprintf('=== Wing / Body Closeout Demo ===\n');

    haveSSI = exist('geom.SurfaceSurfaceIntersect', 'class') == 8;
    havePatchOps = exist('geom.PatchOps', 'class') == 8;
    haveMeshWriter = exist('geom.MeshWriter', 'class') == 8;
    if ~haveSSI || ~havePatchOps
        error(['This demo requires geom.SurfaceSurfaceIntersect and ' ...
               'geom.PatchOps on the MATLAB path.']);
    end

    %% 1) Build fuselage surface
    fprintf('\n--- 1. Build fuselage surface ---\n');

    xprof = linspace(0, 12, 13)';
    rprof = [0.05; 0.35; 0.70; 1.00; 1.18; 1.22; 1.18; 1.05; 0.88; 0.65; 0.40; 0.18; 0.03];
    profile_pts = [xprof, rprof, zeros(size(xprof))];

    S_body = geom.Aircraft.revolutionSurface(profile_pts, [0 0 0], [1 0 0], 0, 360);
    fprintf(' Body control net: %d x %d\n', size(S_body.P,1), size(S_body.P,2));

    %% 2) Build half-wing
    fprintf('\n--- 2. Build half-wing surface ---\n');

    [x_af, y_af] = geom.Aircraft.naca5Coords(23012, 60);
    airfoil_xy = [x_af, y_af];

    spans     = [0; 1.2; 2.6; 4.2; 6.0];
    chords    = [3.0; 2.7; 2.2; 1.7; 1.2];
    sweeps    = [1.0; 1.15; 1.45; 1.9; 2.5];
    twists    = [0; -1.0; -2.0; -3.0; -4.0];
    dihedrals = [0.15; 0.25; 0.38; 0.56; 0.80];

    S_wing = geom.Aircraft.wingSurface(airfoil_xy, spans, chords, sweeps, twists, dihedrals);
    if isa(S_wing, 'geom.LoftedSurface')
        S_wing = S_wing.surface;
    end
    fprintf(' Wing control net: %d x %d\n', size(S_wing.P,1), size(S_wing.P,2));

    %% 3) Intersect wing and body
    fprintf('\n--- 3. Surface-surface intersection ---\n');

    hits = geom.SurfaceSurfaceIntersect.intersect(S_wing, S_body, ...
        'Nu', 18, 'Nv', 12, 'TolPoint', 1e-5, 'Refine', true);

    nh = numel(hits);
    fprintf(' Intersection hits found: %d\n', nh);
    if nh < 4
        error(['Intersection produced too few points. Increase sampling ' ...
               'density or adjust geometry.']);
    end

    pts = reshape([hits.pt], 3, []).';
    uw = [hits.u1].';
    vw = [hits.v1].';
    ub = [hits.u2].';
    vb = [hits.v2].';

    residuals = [hits.residual];
    fprintf(' Best residual: %.3e\n', min(residuals));
    fprintf(' Worst residual: %.3e\n', max(residuals));

    %% 4) Split fuselage into six patches as before
    fprintf('\n--- 4. Split fuselage into six patches ---\n');

    U0b = S_body.domainU(1); U1b = S_body.domainU(2);
    V0b = S_body.domainV(1); V1b = S_body.domainV(2);

    u_fwd = min(ub);
    u_aft = max(ub);
    du_margin = 0.01 * max(U1b - U0b, 1.0);
    u_fwd = max(U0b + 1e-8, u_fwd - du_margin);
    u_aft = min(U1b - 1e-8, u_aft + du_margin);

    % For this revolution surface, v_mid ~ half-wrap gives upper/lower split
    v_mid_body = 0.5 * (V0b + V1b);

    fus = struct();
    fus.fwd_upper = geom.PatchOps.extractUV(S_body, [U0b, u_fwd], [V0b, v_mid_body]);
    fus.fwd_lower = geom.PatchOps.extractUV(S_body, [U0b, u_fwd], [v_mid_body, V1b]);
    fus.mid_upper = geom.PatchOps.extractUV(S_body, [u_fwd, u_aft], [V0b, v_mid_body]);
    fus.mid_lower = geom.PatchOps.extractUV(S_body, [u_fwd, u_aft], [v_mid_body, V1b]);
    fus.aft_upper = geom.PatchOps.extractUV(S_body, [u_aft, U1b], [V0b, v_mid_body]);
    fus.aft_lower = geom.PatchOps.extractUV(S_body, [u_aft, U1b], [v_mid_body, V1b]);

    %% 5) Split wing into root and outboard patches
    fprintf('\n--- 5. Split wing into root/outboard patches ---\n');

    U0w = S_wing.domainU(1); U1w = S_wing.domainU(2);
    V0w = S_wing.domainV(1); V1w = S_wing.domainV(2);

    % Wing-body intersection typically occupies a small root band in the loft direction.
    % We use the extrema of the wing-side intersection parameters to define that band.
    uw_min = min(uw);
    uw_max = max(uw);
    duw_margin = 0.01 * max(U1w - U0w, 1.0);

    uw_root0 = max(U0w + 1e-8, uw_min - duw_margin);
    uw_root1 = min(U1w - 1e-8, uw_max + duw_margin);

    % Split upper/lower using a representative wing v mid.
    % This is only a first-pass convenience split.
    v_mid_wing = 0.5 * (V0w + V1w);

    wing = struct();
    wing.root_upper = geom.PatchOps.extractUV(S_wing, [uw_root0, uw_root1], [V0w, v_mid_wing]);
    wing.root_lower = geom.PatchOps.extractUV(S_wing, [uw_root0, uw_root1], [v_mid_wing, V1w]);

    % Optional outboard bands
    wing.out_upper = geom.PatchOps.extractUV(S_wing, [uw_root1, U1w], [V0w, v_mid_wing]);
    wing.out_lower = geom.PatchOps.extractUV(S_wing, [uw_root1, U1w], [v_mid_wing, V1w]);

    fprintf(' Wing root band u=[%.6f, %.6f]\n', uw_root0, uw_root1);

    %% 6) Build simple boundary curves for closeouts
    fprintf('\n--- 6. Build upper/lower closeout curves ---\n');

    % Use the ordered intersection points and split them by z sign.
    upper_idx = pts(:,3) >= 0;
    lower_idx = pts(:,3) < 0;

    pts_upper = pts(upper_idx, :);
    pts_lower = pts(lower_idx, :);

    % Ensure each branch is sorted roughly streamwise.
    if ~isempty(pts_upper)
        [~, idx] = sort(pts_upper(:,1));
        pts_upper = pts_upper(idx,:);
    end
    if ~isempty(pts_lower)
        [~, idx] = sort(pts_lower(:,1));
        pts_lower = pts_lower(idx,:);
    end

    % Build interpolated intersection curves where possible.
    C_int_upper = [];
    C_int_lower = [];
    try
        if size(pts_upper,1) >= 4
            C_int_upper = geom.LoftedSurface.globalCurveInterp(pts_upper, 3, size(pts_upper,1));
        elseif size(pts_upper,1) >= 2
            C_int_upper = geom.NURBSCurve(pts_upper, 1);
        end
    catch
    end

    try
        if size(pts_lower,1) >= 4
            C_int_lower = geom.LoftedSurface.globalCurveInterp(pts_lower, 3, size(pts_lower,1));
        elseif size(pts_lower,1) >= 2
            C_int_lower = geom.NURBSCurve(pts_lower, 1);
        end
    catch
    end

    % Build a matching wing-side "inner" curve by sampling the root patches near the inboard side.
    % This is a practical first-pass curve for ruled closeout generation.
    C_wing_upper = geom.PatchOps.extractInnerSpanCurve(wing.root_upper);
    C_wing_lower = geom.PatchOps.extractInnerSpanCurve(wing.root_lower);

    % Build a matching fuselage-side curve by sampling the mid patches near the intersection side.
    C_body_upper = geom.PatchOps.extractInnerCircCurve(fus.mid_upper);
    C_body_lower = geom.PatchOps.extractInnerCircCurve(fus.mid_lower);

    % Prefer actual intersection-derived curves when available on the fuselage side.
    if isempty(C_int_upper)
        C_close_upper_body = C_body_upper;
    else
        C_close_upper_body = C_int_upper;
    end
    if isempty(C_int_lower)
        C_close_lower_body = C_body_lower;
    else
        C_close_lower_body = C_int_lower;
    end

    %% 7) Build upper and lower closeout surfaces
    fprintf('\n--- 7. Build ruled closeout surfaces ---\n');

    % Harmonize and create ruled surfaces.
    S_close_upper = geom.Aircraft.ruledSurface(C_close_upper_body, C_wing_upper);
    S_close_lower = geom.Aircraft.ruledSurface(C_close_lower_body, C_wing_lower);

    fprintf(' Upper closeout net: %d x %d\n', size(S_close_upper.P,1), size(S_close_upper.P,2));
    fprintf(' Lower closeout net: %d x %d\n', size(S_close_lower.P,1), size(S_close_lower.P,2));

    %% 8) Mesh everything
    fprintf('\n--- 8. Mesh patches and closeouts ---\n');

    meshes = struct();

    meshes.fus_fwd_upper = fus.fwd_upper.isoMesh(18, 16, 'SpacingU', 'cosine', 'SpacingV', 'uniform');
    meshes.fus_fwd_lower = fus.fwd_lower.isoMesh(18, 16, 'SpacingU', 'cosine', 'SpacingV', 'uniform');
    meshes.fus_mid_upper = fus.mid_upper.isoMesh(24, 16, 'SpacingU', 'uniform', 'SpacingV', 'uniform');
    meshes.fus_mid_lower = fus.mid_lower.isoMesh(24, 16, 'SpacingU', 'uniform', 'SpacingV', 'uniform');
    meshes.fus_aft_upper = fus.aft_upper.isoMesh(16, 16, 'SpacingU', 'cosine', 'SpacingV', 'uniform');
    meshes.fus_aft_lower = fus.aft_lower.isoMesh(16, 16, 'SpacingU', 'cosine', 'SpacingV', 'uniform');

    meshes.wing_root_upper = wing.root_upper.isoMesh(18, 14, 'SpacingU', 'uniform', 'SpacingV', 'uniform');
    meshes.wing_root_lower = wing.root_lower.isoMesh(18, 14, 'SpacingU', 'uniform', 'SpacingV', 'uniform');
    meshes.wing_out_upper  = wing.out_upper.isoMesh(24, 14, 'SpacingU', 'uniform', 'SpacingV', 'uniform');
    meshes.wing_out_lower  = wing.out_lower.isoMesh(24, 14, 'SpacingU', 'uniform', 'SpacingV', 'uniform');

    meshes.close_upper = S_close_upper.isoMesh(20, 10, 'SpacingU', 'uniform', 'SpacingV', 'uniform');
    meshes.close_lower = S_close_lower.isoMesh(20, 10, 'SpacingU', 'uniform', 'SpacingV', 'uniform');

    %% 9) Optional network packaging
    fprintf('\n--- 9. Package WGS-style networks (optional) ---\n');
    networks = {};
    if haveMeshWriter
        labels = { ...
            'FusFwdUpper','FusFwdLower','FusMidUpper','FusMidLower','FusAftUpper','FusAftLower', ...
            'WingRootUpper','WingRootLower','WingOutUpper','WingOutLower', ...
            'CloseUpper','CloseLower'};
        meshNames = fieldnames(meshes);
        for k = 1:numel(meshNames)
            networks{k} = geom.MeshWriter.meshToNetwork(meshes.(meshNames{k}), labels{k}, 0, 0); %#ok<AGROW>
        end
        fprintf(' Created %d network structs.\n', numel(networks));
    else
        fprintf(' geom.MeshWriter not found; skipping network packaging.\n');
    end

    %% 10) Plot geometry decomposition
    fprintf('\n--- 10. Plot geometry and closeouts ---\n');

    figure('Name', 'Wing / Body Closeout Geometry');
    hold on;

    % Fuselage patches
    fusColors = [ ...
        0.90 0.75 0.30;
        0.95 0.60 0.30;
        0.30 0.80 0.45;
        0.20 0.65 0.40;
        0.75 0.55 0.90;
        0.60 0.45 0.80];
    fusNames = {'fwd_upper','fwd_lower','mid_upper','mid_lower','aft_upper','aft_lower'};
    for k = 1:numel(fusNames)
        fus.(fusNames{k}).plot(16, 12, 'ShowCP', false, ...
            'Alpha', 0.75, 'EdgeAlpha', 0.15, 'FaceColor', fusColors(k,:));
    end

    % Wing patches
    wing.root_upper.plot(22, 12, 'ShowCP', false, 'Alpha', 0.85, 'EdgeAlpha', 0.15, 'FaceColor', [0.35 0.55 0.95]);
    wing.root_lower.plot(22, 12, 'ShowCP', false, 'Alpha', 0.85, 'EdgeAlpha', 0.15, 'FaceColor', [0.25 0.45 0.85]);
    wing.out_upper.plot(22, 12, 'ShowCP', false, 'Alpha', 0.45, 'EdgeAlpha', 0.10, 'FaceColor', [0.55 0.70 1.00]);
    wing.out_lower.plot(22, 12, 'ShowCP', false, 'Alpha', 0.45, 'EdgeAlpha', 0.10, 'FaceColor', [0.45 0.60 0.92]);

    % Closeouts
    S_close_upper.plot(16, 8, 'ShowCP', false, 'Alpha', 0.95, 'EdgeAlpha', 0.25, 'FaceColor', [0.95 0.20 0.20]);
    S_close_lower.plot(16, 8, 'ShowCP', false, 'Alpha', 0.95, 'EdgeAlpha', 0.25, 'FaceColor', [0.80 0.10 0.10]);

    % Intersection points
    plot3(pts(:,1), pts(:,2), pts(:,3), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 5);

    title('Wing / Body closeout workflow');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal;
    grid on;
    view(35, 22);

    %% 11) Plot structured meshes
    fprintf('\n--- 11. Plot structured meshes ---\n');

    figure('Name', 'Wing / Body Closeout Meshes');
    hold on;

    meshNames = fieldnames(meshes);
    meshColors = [
        0.90 0.75 0.30;
        0.95 0.60 0.30;
        0.30 0.80 0.45;
        0.20 0.65 0.40;
        0.75 0.55 0.90;
        0.60 0.45 0.80;
        0.35 0.55 0.95;
        0.25 0.45 0.85;
        0.55 0.70 1.00;
        0.45 0.60 0.92;
        0.95 0.20 0.20;
        0.80 0.10 0.10];
    for k = 1:numel(meshNames)
        geom.PlotHelpers.meshStruct(meshes.(meshNames{k}), meshColors(k,:));
    end
    plot3(pts(:,1), pts(:,2), pts(:,3), 'k.', 'MarkerSize', 12);

    title('Structured meshes on fuselage, wing, and closeout patches');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal;
    grid on;
    view(35, 22);

    %% 12) Summary
    fprintf('\n=== Demo Summary ===\n');
    fprintf('Built:\n');
    fprintf('  - six fuselage patches\n');
    fprintf('  - wing root and outboard patches\n');
    fprintf('  - upper and lower ruled closeout surfaces\n');
    fprintf('\nThis is a practical first-pass wing/body patch network.\n');
    fprintf('Next refinements would be:\n');
    fprintf('  1) replace ruled closeouts with Coons/tangent-matched patches,\n');
    fprintf('  2) use true branch-separated intersection curves for upper/lower,\n');
    fprintf('  3) enforce matching edge point distributions before WGS export.\n');

    if haveMeshWriter
        fprintf('\nYou can export the networks with something like:\n');
        fprintf('  geom.MeshWriter.toWGS(''wing_body_closeout.wgs'', networks)\n');
    end
end
