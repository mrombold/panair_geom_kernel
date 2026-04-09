function demo_simple_algorithm_validation()
% DEMO_SIMPLE_ALGORITHM_VALIDATION
% A simpler, more controlled validation of the new geometry algorithms.
%
% This demo checks three things independently:
%   1) exact surface splitting / extraction on a simple fuselage-like surface
%   2) surface-surface intersection on a known planar cut through a body
%   3) ruled closeout generation between two simple, known curves
%
% The goal is to validate the algorithms before using them in the more
% complicated wing/body junction workflow.

    close all;
    clc;

    fprintf('=== Simple Algorithm Validation ===\n');

    if exist('geom.PatchOps', 'class') ~= 8 || ...
       exist('geom.SurfaceSurfaceIntersect', 'class') ~= 8
        error('This demo requires geom.PatchOps and geom.SurfaceSurfaceIntersect.');
    end

    %% --------------------------------------------------------------------
    % 1) Build a simple body of revolution
    % ---------------------------------------------------------------------
    fprintf('\n--- 1. Simple body surface ---\n');

    xprof = linspace(0, 10, 11)';
    rprof = [0.05; 0.50; 0.90; 1.10; 1.20; 1.20; 1.10; 0.90; 0.60; 0.25; 0.05];
    profile_pts = [xprof, rprof, zeros(size(xprof))];

    S_body = geom.Aircraft.revolutionSurface(profile_pts, [0 0 0], [1 0 0], 0, 360);
    fprintf(' Body net: %d x %d\n', size(S_body.P,1), size(S_body.P,2));

    %% --------------------------------------------------------------------
    % 2) Validate exact splitting / extraction with known parameters
    % ---------------------------------------------------------------------
    fprintf('\n--- 2. Exact split / extract validation ---\n');

    u_fwd = 0.30;
    u_aft = 0.70;
    v_mid = 0.50;

    S_fwd_upper = geom.PatchOps.extractUV(S_body, [0.00, u_fwd], [0.00, v_mid]);
    S_mid_upper = geom.PatchOps.extractUV(S_body, [u_fwd, u_aft], [0.00, v_mid]);
    S_aft_upper = geom.PatchOps.extractUV(S_body, [u_aft, 1.00], [0.00, v_mid]);

    fprintf(' fwd upper net: %d x %d\n', size(S_fwd_upper.P,1), size(S_fwd_upper.P,2));
    fprintf(' mid upper net: %d x %d\n', size(S_mid_upper.P,1), size(S_mid_upper.P,2));
    fprintf(' aft upper net: %d x %d\n', size(S_aft_upper.P,1), size(S_aft_upper.P,2));

    % Check that extracted geometry matches parent surface at mapped sample points
    % We'll sample the mid-upper patch.
    uv_sub = [ ...
        0.15 0.10;
        0.35 0.25;
        0.65 0.40;
        0.85 0.15];

    err_extract = zeros(size(uv_sub,1),1);
    for k = 1:size(uv_sub,1)
        us = uv_sub(k,1);
        vs = uv_sub(k,2);

        % Map local subpatch coordinates back to parent body coordinates
        up = u_fwd + us * (u_aft - u_fwd);
        vp = 0.00 + vs * (v_mid - 0.00);

        p_parent = S_body.evaluate(up, vp);
        p_sub    = S_mid_upper.evaluate(us, vs);

        err_extract(k) = norm(p_parent - p_sub);
    end

    fprintf(' max extractUV geometry error: %.3e\n', max(err_extract));

    %% --------------------------------------------------------------------
    % 3) Validate surface-surface intersection with a very simple cut surface
    % ---------------------------------------------------------------------
    fprintf('\n--- 3. Surface-surface intersection with simple planar cut ---\n');

    % Build a planar rectangular NURBS surface at y = +0.6 spanning x,z.
    % This cuts the body in a predictable vertical section.
    x0 = 2.0; x1 = 8.0;
    z0 = -1.8; z1 = 1.8;
    ycut = 0.60;

    Pplane = zeros(2,2,3);
    Pplane(1,1,:) = [x0, ycut, z0];
    Pplane(2,1,:) = [x1, ycut, z0];
    Pplane(1,2,:) = [x0, ycut, z1];
    Pplane(2,2,:) = [x1, ycut, z1];

    S_cut = geom.NURBSSurface(Pplane, 1, 1, [0 0 1 1], [0 0 1 1], ones(2,2));

    hits = geom.SurfaceSurfaceIntersect.intersect(S_cut, S_body, ...
        'Nu', 12, 'Nv', 8, 'TolPoint', 1e-6, 'Refine', true);

    fprintf(' intersection hits found: %d\n', numel(hits));
    if isempty(hits)
        warning('No hits found in simple validation case.');
    else
        pts = reshape([hits.pt], 3, []).';
        y_err = max(abs(pts(:,2) - ycut));
        x_min = min(pts(:,1));
        x_max = max(pts(:,1));

        fprintf(' max |y - ycut| on returned points: %.3e\n', y_err);
        fprintf(' x-range of intersection points: [%.4f, %.4f]\n', x_min, x_max);
        fprintf(' best residual: %.3e\n', min([hits.residual]));
        fprintf(' worst residual: %.3e\n', max([hits.residual]));
    end

    %% --------------------------------------------------------------------
    % 4) Validate ruled closeout on two simple known curves
    % ---------------------------------------------------------------------
    fprintf('\n--- 4. Ruled closeout validation with simple curves ---\n');

    % Two simple curves with same parameter direction:
    %   lower curve: straight line
    %   upper curve: same x-y trace but raised and bowed in z
    t = linspace(0,1,7)';
    C1_pts = [2 + 4*t, 0.8*ones(size(t)), zeros(size(t))];
    C2_pts = [2 + 4*t, 1.6*ones(size(t)), 0.3*sin(pi*t)];

    C1 = geom.LoftedSurface.globalCurveInterp(C1_pts, 3, size(C1_pts,1));
    C2 = geom.LoftedSurface.globalCurveInterp(C2_pts, 3, size(C2_pts,1));

    S_close = geom.Aircraft.ruledSurface(C1, C2);
    fprintf(' ruled closeout net: %d x %d\n', size(S_close.P,1), size(S_close.P,2));

    % Check boundary reproduction
    svals = linspace(0,1,9)';
    err_bnd_1 = zeros(numel(svals),1);
    err_bnd_2 = zeros(numel(svals),1);
    for k = 1:numel(svals)
        pk1 = C1.evaluate(svals(k));
        ps1 = S_close.evaluate(svals(k), 0.0);
        err_bnd_1(k) = norm(pk1 - ps1);

        pk2 = C2.evaluate(svals(k));
        ps2 = S_close.evaluate(svals(k), 1.0);
        err_bnd_2(k) = norm(pk2 - ps2);
    end
    fprintf(' max boundary error at v=0: %.3e\n', max(err_bnd_1));
    fprintf(' max boundary error at v=1: %.3e\n', max(err_bnd_2));

    %% --------------------------------------------------------------------
    % 5) Plot everything
    % ---------------------------------------------------------------------
    fprintf('\n--- 5. Plot validation cases ---\n');

    figure('Name','Simple Validation: Body split and planar cut');
    hold on;
    S_body.plot(30, 20, 'ShowCP', false, 'Alpha', 0.25, 'EdgeAlpha', 0.08, 'FaceColor', [0.7 0.75 0.85]);
    S_fwd_upper.plot(14, 10, 'ShowCP', false, 'Alpha', 0.90, 'EdgeAlpha', 0.18, 'FaceColor', [0.95 0.70 0.30]);
    S_mid_upper.plot(14, 10, 'ShowCP', false, 'Alpha', 0.90, 'EdgeAlpha', 0.18, 'FaceColor', [0.30 0.85 0.40]);
    S_aft_upper.plot(14, 10, 'ShowCP', false, 'Alpha', 0.90, 'EdgeAlpha', 0.18, 'FaceColor', [0.75 0.55 0.90]);
    S_cut.plot(8, 8, 'ShowCP', false, 'Alpha', 0.35, 'EdgeAlpha', 0.25, 'FaceColor', [0.35 0.55 0.95]);

    if ~isempty(hits)
        pts = reshape([hits.pt], 3, []).';
        plot3(pts(:,1), pts(:,2), pts(:,3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
        plot3(pts(:,1), pts(:,2), pts(:,3), 'r-', 'LineWidth', 1.5);
    end

    title('Simple split validation and planar body cut');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; grid on; view(35,22);

    figure('Name','Simple Validation: Ruled closeout');
    hold on;
    S_close.plot(24, 8, 'ShowCP', false, 'Alpha', 0.85, 'EdgeAlpha', 0.20, 'FaceColor', [0.95 0.25 0.25]);
    C1.plot(80, 'ShowCP', true, 'Color', [0 0 0], 'LineWidth', 1.5);
    C2.plot(80, 'ShowCP', true, 'Color', [0 0 0], 'LineWidth', 1.5);

    title('Ruled closeout between two known curves');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; grid on; view(35,22);

    %% --------------------------------------------------------------------
    % 6) Summary
    % ---------------------------------------------------------------------
    fprintf('\n=== Validation Summary ===\n');
    fprintf('What to look for:\n');
    fprintf('  - extractUV error should be near machine precision\n');
    fprintf('  - intersection points should satisfy y = %.3f to tight tolerance\n', ycut);
    fprintf('  - ruled closeout should reproduce the two input curves exactly on its boundaries\n');
    fprintf('\nIf these all look good, then the core algorithms are behaving,\n');
    fprintf('and any issue in the wing/body demo is likely in the higher-level setup,\n');
    fprintf('not in split/extract, simple intersection, or ruled-surface creation.\n');
end
