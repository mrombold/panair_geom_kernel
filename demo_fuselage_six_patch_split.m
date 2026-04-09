function demo_fuselage_six_patch_split()
% DEMO_FUSELAGE_SIX_PATCH_SPLIT
% Build a fuselage-like revolution surface and a half-wing, compute the
% wing/body intersection, then decompose the fuselage into six structured
% patches:
%
%   fwd upper, fwd lower,
%   mid upper, mid lower,
%   aft upper, aft lower
%
% The decomposition is done by exact parameter-space splitting/extraction
% rather than general CAD trimming, which is better suited to structured
% PanAir-style patch meshing.
%
% Requirements:
%   +geom/NURBSSurface.m         (patched splitU/splitV recommended)
%   +geom/PatchOps.m
%   +geom/SurfaceSurfaceIntersect.m
%   +geom/Aircraft.m
%   +geom/MeshWriter.m          (optional, only for WGS network packaging)

    close all;
    clc;

    fprintf('=== Fuselage Six-Patch Split Demo ===\n');

    haveSSI = exist('geom.SurfaceSurfaceIntersect', 'class') == 8;
    havePatchOps = exist('geom.PatchOps', 'class') == 8;
    haveMeshWriter = exist('geom.MeshWriter', 'class') == 8;

    if ~haveSSI || ~havePatchOps
        error(['This demo requires geom.SurfaceSurfaceIntersect and ' ...
               'geom.PatchOps on the MATLAB path.']);
    end

    %% 1) Build fuselage-like revolution surface
    fprintf('\n--- 1. Build fuselage surface ---\n');

    xprof = linspace(0, 12, 13)';
    rprof = [0.05; 0.35; 0.70; 1.00; 1.18; 1.22; 1.18; 1.05; 0.88; 0.65; 0.40; 0.18; 0.03];
    profile_pts = [xprof, rprof, zeros(size(xprof))];

    S_body = geom.Aircraft.revolutionSurface(profile_pts, [0 0 0], [1 0 0], 0, 360);
    fprintf(' Body control net: %d x %d\n', size(S_body.P,1), size(S_body.P,2));
    fprintf(' Body domains: U=[%.4f, %.4f], V=[%.4f, %.4f]\n', ...
        S_body.domainU(1), S_body.domainU(2), S_body.domainV(1), S_body.domainV(2));

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
    if nh < 2
        error(['Intersection produced too few points. Try increasing ' ...
               'sampling density or adjusting geometry.']);
    end

    residuals = [hits.residual];
    fprintf(' Best residual: %.3e\n', min(residuals));
    fprintf(' Worst residual: %.3e\n', max(residuals));

    pts = reshape([hits.pt], 3, []).';
    ub = [hits.u2].';
    vb = [hits.v2].';

    % For the revolution surface built here:
    %   U ~ streamwise/axial profile parameter
    %   V ~ circumferential wrap parameter
    %
    % Split forward/aft at the extremal axial parameters of the join.
    u_fwd = min(ub);
    u_aft = max(ub);

    fprintf(' Raw fuselage split params from intersection: u_fwd=%.6f, u_aft=%.6f\n', u_fwd, u_aft);

    % Add a small buffer so the mid patch cleanly contains the join curve.
    du_margin = 0.01 * max(S_body.domainU(2) - S_body.domainU(1), 1.0);
    u_fwd = max(S_body.domainU(1) + 1e-8, u_fwd - du_margin);
    u_aft = min(S_body.domainU(2) - 1e-8, u_aft + du_margin);

    fprintf(' Buffered fuselage split params: u_fwd=%.6f, u_aft=%.6f\n', u_fwd, u_aft);

    % Upper/lower split for a 0..360 revolution:
    % angle 0   -> starboard (z=0)
    % angle pi  -> port     (z=0)
    % so upper half is approximately V in [0, 0.5]
    v_mid = 0.5 * (S_body.domainV(1) + S_body.domainV(2));
    fprintf(' Circumferential split param (upper/lower): v_mid=%.6f\n', v_mid);

    %% 4) Extract the six fuselage patches
    fprintf('\n--- 4. Extract six fuselage patches ---\n');

    U0 = S_body.domainU(1);
    U1 = S_body.domainU(2);
    V0 = S_body.domainV(1);
    V1 = S_body.domainV(2);

    patches = struct();

    % Forward band
    patches.fwd_upper = geom.PatchOps.extractUV(S_body, [U0, u_fwd], [V0, v_mid]);
    patches.fwd_lower = geom.PatchOps.extractUV(S_body, [U0, u_fwd], [v_mid, V1]);

    % Mid band
    patches.mid_upper = geom.PatchOps.extractUV(S_body, [u_fwd, u_aft], [V0, v_mid]);
    patches.mid_lower = geom.PatchOps.extractUV(S_body, [u_fwd, u_aft], [v_mid, V1]);

    % Aft band
    patches.aft_upper = geom.PatchOps.extractUV(S_body, [u_aft, U1], [V0, v_mid]);
    patches.aft_lower = geom.PatchOps.extractUV(S_body, [u_aft, U1], [v_mid, V1]);

    names = fieldnames(patches);
    for k = 1:numel(names)
        S = patches.(names{k});
        fprintf(' %-10s : [%2dx%2d] control net\n', names{k}, size(S.P,1), size(S.P,2));
    end

    %% 5) Mesh each patch
    fprintf('\n--- 5. Mesh the six patches ---\n');

    % Suggestion:
    %   nu = streamwise density for each band
    %   nv = circumferential density within upper/lower half
    %
    % You can tune these separately per patch later.
    nu_fwd = 18;
    nu_mid = 24;
    nu_aft = 16;
    nv_half = 16;

    meshes = struct();

    meshes.fwd_upper = patches.fwd_upper.isoMesh(nu_fwd, nv_half, ...
        'SpacingU', 'cosine', 'SpacingV', 'uniform');
    meshes.fwd_lower = patches.fwd_lower.isoMesh(nu_fwd, nv_half, ...
        'SpacingU', 'cosine', 'SpacingV', 'uniform');

    meshes.mid_upper = patches.mid_upper.isoMesh(nu_mid, nv_half, ...
        'SpacingU', 'uniform', 'SpacingV', 'uniform');
    meshes.mid_lower = patches.mid_lower.isoMesh(nu_mid, nv_half, ...
        'SpacingU', 'uniform', 'SpacingV', 'uniform');

    meshes.aft_upper = patches.aft_upper.isoMesh(nu_aft, nv_half, ...
        'SpacingU', 'cosine', 'SpacingV', 'uniform');
    meshes.aft_lower = patches.aft_lower.isoMesh(nu_aft, nv_half, ...
        'SpacingU', 'cosine', 'SpacingV', 'uniform');

    for k = 1:numel(names)
        M = meshes.(names{k});
        fprintf(' %-10s : mesh %2d x %2d\n', names{k}, M.nu, M.nv);
    end

    %% 6) Optional network packaging for WGS export
    fprintf('\n--- 6. Package as WGS-style networks (optional) ---\n');
    networks = {};
    if haveMeshWriter
        labels = {'FwdUpper','FwdLower','MidUpper','MidLower','AftUpper','AftLower'};
        for k = 1:numel(names)
            networks{k} = geom.MeshWriter.meshToNetwork(meshes.(names{k}), labels{k}, 0, 0); %#ok<AGROW>
        end
        fprintf(' Created %d network structs for MeshWriter.toWGS().\n', numel(networks));
    else
        fprintf(' geom.MeshWriter not found; skipping network packaging.\n');
    end

    %% 7) Plot result
    fprintf('\n--- 7. Plot decomposition and meshes ---\n');

    figure('Name', 'Fuselage Six-Patch Split');
    hold on;

    % Plot full wing semi-transparent
    S_wing.plot(28, 12, 'ShowCP', false, 'Alpha', 0.78, ...
        'EdgeAlpha', 0.10, 'FaceColor', [0.35 0.55 0.95]);

    % Plot six fuselage patches with different colors
    patchColors = [ ...
        0.90 0.75 0.30;   % fwd upper
        0.95 0.60 0.30;   % fwd lower
        0.30 0.80 0.45;   % mid upper
        0.20 0.65 0.40;   % mid lower
        0.75 0.55 0.90;   % aft upper
        0.60 0.45 0.80];  % aft lower

    for k = 1:numel(names)
        patches.(names{k}).plot(18, 12, 'ShowCP', false, ...
            'Alpha', 0.85, 'EdgeAlpha', 0.22, 'FaceColor', patchColors(k,:));
    end

    % Plot intersection points
    plot3(pts(:,1), pts(:,2), pts(:,3), 'ko', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 5);

    title('Fuselage split into forward/mid/aft and upper/lower patches');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal;
    grid on;
    view(35, 22);

    %% 8) Plot the actual structured meshes
    figure('Name', 'Fuselage Six-Patch Meshes');
    hold on;

    geom.PlotHelpers.meshStruct(meshes.fwd_upper, patchColors(1,:));
    geom.PlotHelpers.meshStruct(meshes.fwd_lower, patchColors(2,:));
    geom.PlotHelpers.meshStruct(meshes.mid_upper, patchColors(3,:));
    geom.PlotHelpers.meshStruct(meshes.mid_lower, patchColors(4,:));
    geom.PlotHelpers.meshStruct(meshes.aft_upper, patchColors(5,:));
    geom.PlotHelpers.meshStruct(meshes.aft_lower, patchColors(6,:));

    plot3(pts(:,1), pts(:,2), pts(:,3), 'k.', 'MarkerSize', 12);

    title('Structured meshes on six fuselage patches');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal;
    grid on;
    view(35, 22);

    %% 9) Summary
    fprintf('\n=== Demo Summary ===\n');
    fprintf('Created fuselage patches:\n');
    fprintf('  fwd_upper, fwd_lower,\n');
    fprintf('  mid_upper, mid_lower,\n');
    fprintf('  aft_upper, aft_lower\n');
    fprintf('\nThis is a split-based decomposition, not general CAD trimming.\n');
    fprintf('That is usually the better first move for PanAir-style structured grids.\n');

    if haveMeshWriter
        fprintf('\nYou can now export the networks with something like:\n');
        fprintf('  geom.MeshWriter.toWGS(''fuselage_six_patch.wgs'', networks)\n');
    end

end
