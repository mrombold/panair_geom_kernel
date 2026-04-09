function demo_wing_body_workflow()
% DEMO_WING_BODY_WORKFLOW
% End-to-end demonstration for:
%   1) building a fuselage-like revolution surface
%   2) building a half-wing surface
%   3) intersecting wing and body
%   4) extracting a local body patch around the junction
%   5) plotting the result
%
% Intended as a practical validation script for the PanAir geometry kernel.

    close all;
    clc;

    fprintf('=== Wing/Body Workflow Demo ===\n');

    if exist('+geom', 'dir')
        fprintf('geom package detected in current working tree.\n');
    end

    %% 1) Build fuselage-like body from a meridian and revolution
    fprintf('\n--- 1. Build fuselage surface ---\n');

    xprof = linspace(0, 12, 13)';
    rprof = [0.05; 0.35; 0.70; 1.00; 1.18; 1.22; 1.18; 1.05; 0.88; 0.65; 0.40; 0.18; 0.03];
    profile_pts = [xprof, rprof, zeros(size(xprof))];  % lies in XY plane, revolve about X axis

    S_body = geom.Aircraft.revolutionSurface(profile_pts, [0 0 0], [1 0 0], 0, 360);
    fprintf(' Body control net: %d x %d\n', size(S_body.P,1), size(S_body.P,2));

    %% 2) Build simple half-wing
    fprintf('\n--- 2. Build half-wing surface ---\n');

    % Use a simple NACA 23012 section from the aircraft helpers
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

    %% 3) Surface-surface intersection
    fprintf('\n--- 3. Surface-surface intersection ---\n');
    hits = geom.SurfaceSurfaceIntersect.intersect(S_wing, S_body, ...
        'Nu', 18, 'Nv', 12, 'TolPoint', 1e-5, 'Refine', true);

    fprintf(' Hits found: %d\n', numel(hits));
    if isempty(hits)
        warning('No intersections found. Try increasing sampling density or changing geometry.');
        return;
    end

    residuals = [hits.residual];
    fprintf(' Best residual: %.3e\n', min(residuals));
    fprintf(' Worst residual: %.3e\n', max(residuals));

    uvWing = [[hits.u1].', [hits.v1].'];
    uvBody = [[hits.u2].', [hits.v2].'];

    fprintf(' Wing parameter box: u=[%.4f, %.4f], v=[%.4f, %.4f]\n', ...
        min(uvWing(:,1)), max(uvWing(:,1)), min(uvWing(:,2)), max(uvWing(:,2)));
    fprintf(' Body parameter box: u=[%.4f, %.4f], v=[%.4f, %.4f]\n', ...
        min(uvBody(:,1)), max(uvBody(:,1)), min(uvBody(:,2)), max(uvBody(:,2)));

    %% 4) Extract a local body patch around the intersection
    fprintf('\n--- 4. Extract local body patch around junction ---\n');

    du = 0.03 * max(S_body.domainU(2)-S_body.domainU(1), 1);
    dv = 0.06 * max(S_body.domainV(2)-S_body.domainV(1), 1);

    ur = [max(S_body.domainU(1), min(uvBody(:,1)) - du), ...
          min(S_body.domainU(2), max(uvBody(:,1)) + du)];
    vr = [max(S_body.domainV(1), min(uvBody(:,2)) - dv), ...
          min(S_body.domainV(2), max(uvBody(:,2)) + dv)];

    % PatchOps expects normalized parameter behavior because splitU/splitV renormalize
    % child patches. Body revolution surface is built on normalized domains, so this is OK.
    S_body_local = geom.PatchOps.extractUV(S_body, ur, vr);
    fprintf(' Local body patch net: %d x %d\n', size(S_body_local.P,1), size(S_body_local.P,2));

    %% 5) Optional smooth intersection curve
    fprintf('\n--- 5. Build interpolated intersection curve ---\n');
    C_int = [];
    try
        C_int = geom.SurfaceSurfaceIntersect.toCompositeCurve(hits);
        fprintf(' Composite intersection curve built.\n');
    catch ME
        fprintf(' Could not build smooth intersection curve: %s\n', ME.message);
    end

    %% 6) Plot
    fprintf('\n--- 6. Plot workflow ---\n');
    figure('Name', 'Wing / Body Workflow Demo');
    hold on;

    % Full body and wing
    S_body.plot(36, 24, 'ShowCP', false, 'Alpha', 0.25, 'EdgeAlpha', 0.08, 'FaceColor', [0.65 0.70 0.80]);
    S_wing.plot(28, 12, 'ShowCP', false, 'Alpha', 0.78, 'EdgeAlpha', 0.12, 'FaceColor', [0.35 0.55 0.95]);

    % Local extracted body patch
    S_body_local.plot(20, 14, 'ShowCP', false, 'Alpha', 0.9, 'EdgeAlpha', 0.2, 'FaceColor', [0.95 0.75 0.35]);

    % Intersection points
    pts = reshape([hits.pt], 3, []).';
    plot3(pts(:,1), pts(:,2), pts(:,3), 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 5);

    % Smooth curve if available
    if ~isempty(C_int)
        try
            C_int.plot(120, 'Color', [0 0 0], 'LineWidth', 2.0, 'ShowCP', false);
        catch
            % CompositeCurve plot signature may vary by repo version
            xyz = pts;
            plot3(xyz(:,1), xyz(:,2), xyz(:,3), 'k-', 'LineWidth', 2.0);
        end
    else
        plot3(pts(:,1), pts(:,2), pts(:,3), 'k-', 'LineWidth', 1.5);
    end

    title('Wing / Body Surface Intersection Workflow');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal;
    grid on;
    view(35, 22);

    %% 7) Summary
    fprintf('\n=== Demo Summary ===\n');
    fprintf('This demo validates the basic wing/body workflow:\n');
    fprintf('  - surface generation\n');
    fprintf('  - surface/surface intersection sampling + refinement\n');
    fprintf('  - local patch extraction around the junction\n');
    fprintf('  - plotting of the resulting cut region\n');
    fprintf('\nSuggested next steps:\n');
    fprintf('  1) add marching/tracing so the intersection curve is not only sample-based,\n');
    fprintf('  2) extract matching wing-side local patches,\n');
    fprintf('  3) add 4-sided closeout patch constructors,\n');
    fprintf('  4) export the resulting patch network to WGS.\n');
end
