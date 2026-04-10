function demo_wing_upper_lower()
% DEMO_WING_UPPER_LOWER
% Validate separate upper/lower wing surface generation from Selig-format airfoil data.

    close all;
    clc;

    script_dir = fileparts(mfilename('fullpath'));
    if ~isempty(script_dir)
        addpath(script_dir);
    end

    fprintf('=== Wing Upper/Lower Surface Demo ===\n');

    if exist('geom.wingSurfaceUpperLower', 'file') ~= 2
        error('Need geom.wingSurfaceUpperLower on the MATLAB path.');
    end
    if exist('geom.splitAirfoilSelig', 'file') ~= 2
        error('Need geom.splitAirfoilSelig on the MATLAB path.');
    end

    [x_af, y_af] = geom.Aircraft.naca5Coords(23012, 80);
    airfoil_xy = [x_af, y_af];

    spans     = [0; 1.2; 2.6; 4.2; 6.0];
    chords    = [3.0; 2.7; 2.2; 1.7; 1.2];
    sweeps    = [1.0; 1.15; 1.45; 1.9; 2.5];
    twists    = [0; -1.0; -2.0; -3.0; -4.0];
    dihedrals = [0.15; 0.25; 0.38; 0.56; 0.80];

    [S_up, S_lo, details] = geom.wingSurfaceUpperLower( ...
        airfoil_xy, spans, chords, sweeps, twists, dihedrals, ...
        'CloseTE', true, 'TEPoint', 'average', 'NCP', 20);

    if isa(S_up, 'geom.LoftedSurface'), S_up_plot = S_up.surface; else, S_up_plot = S_up; end
    if isa(S_lo, 'geom.LoftedSurface'), S_lo_plot = S_lo.surface; else, S_lo_plot = S_lo; end

    fprintf(' Upper wing net: %d x %d\n', size(S_up_plot.P,1), size(S_up_plot.P,2));
    fprintf(' Lower wing net: %d x %d\n', size(S_lo_plot.P,1), size(S_lo_plot.P,2));

    xy_up = details.upper_xy{1};
    xy_lo = details.lower_xy{1};
    te_gap = norm(xy_up(1,:) - xy_lo(end,:));
    fprintf(' Root section TE gap after closure: %.3e\n', te_gap);

    figure('Name', 'Wing upper/lower surfaces');
    hold on;
    S_up_plot.plot(32, 12, 'ShowCP', false, 'Alpha', 0.85, 'EdgeAlpha', 0.12, 'FaceColor', [0.35 0.55 0.95]);
    S_lo_plot.plot(32, 12, 'ShowCP', false, 'Alpha', 0.85, 'EdgeAlpha', 0.12, 'FaceColor', [0.90 0.40 0.40]);
    title('Separate upper and lower wing surfaces');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; grid on; view(35,22);

    figure('Name', 'Root airfoil branches');
    hold on;
    plot(xy_up(:,1), xy_up(:,2), '-o', 'LineWidth', 1.5);
    plot(xy_lo(:,1), xy_lo(:,2), '-o', 'LineWidth', 1.5);
    legend('Upper: TE->LE', 'Lower: LE->TE', 'Location', 'Best');
    title('Split Selig airfoil branches');
    xlabel('x/c'); ylabel('z/c');
    axis equal; grid on;

    fprintf('\n=== Summary ===\n');
    fprintf('This builds the wing as two separate surfaces from Selig airfoil input.\n');
    fprintf('That makes wing-body intersection and fuselage splitting much easier.\n');
end
