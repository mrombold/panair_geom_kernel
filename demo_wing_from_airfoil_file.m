function demo_wing_from_airfoil_file()
% DEMO_WING_FROM_AIRFOIL_FILE
% Read an airfoil file, build separate upper/lower wing surfaces with a
% sharp trailing edge, and return the 3D TE curve for wake-shedding use.
%
% Before running:
%   set airfoil_file below to a valid .dat file on your machine.

    close all;
    clc;

    fprintf('=== Wing From Airfoil File Demo ===\n');

    rehash
    if isempty(which('geom.readAirfoilFile'))
        error('Need geom.readAirfoilFile on the MATLAB path.');
    end
    if isempty(which('geom.wingSurfaceUpperLowerFromFile'))
        error('Need geom.wingSurfaceUpperLowerFromFile on the MATLAB path.');
    end

    % ---- USER: set this to a real airfoil file path ----
    airfoil_file = 'airfoil.dat';

    if ~isfile(airfoil_file)
        error(['Set airfoil_file inside demo_wing_from_airfoil_file.m to a real file. ' ...
               'Current value not found: %s'], airfoil_file);
    end

    [xy, meta] = geom.readAirfoilFile(airfoil_file, 'NormalizeChord', true, 'SharpTE', true);
    fprintf(' Read file: %s\n', meta.filename);
    fprintf(' Points: %d\n', meta.npts);
    fprintf(' TE endpoint gap after sharpening: %.3e\n', norm(xy(1,:) - xy(end,:)));

    spans     = [0; 1.2; 2.6; 4.2; 6.0];
    chords    = [3.0; 2.7; 2.2; 1.7; 1.2];
    sweeps    = [1.0; 1.15; 1.45; 1.9; 2.5];
    twists    = [0; -1.0; -2.0; -3.0; -4.0];
    dihedrals = [0.15; 0.25; 0.38; 0.56; 0.80];

    [S_up, S_lo, C_te, details] = geom.wingSurfaceUpperLowerFromFile( ...
        airfoil_file, spans, chords, sweeps, twists, dihedrals, ...
        'NormalizeChord', true, 'SharpTE', true, 'TEPoint', 'average');

    if isa(S_up, 'geom.LoftedSurface'), S_up_plot = S_up.surface; else, S_up_plot = S_up; end
    if isa(S_lo, 'geom.LoftedSurface'), S_lo_plot = S_lo.surface; else, S_lo_plot = S_lo; end

    fprintf(' Upper wing net: %d x %d\n', size(S_up_plot.P,1), size(S_up_plot.P,2));
    fprintf(' Lower wing net: %d x %d\n', size(S_lo_plot.P,1), size(S_lo_plot.P,2));

    te_root = details.te_pts(1,:);
    te_tip  = details.te_pts(end,:);
    fprintf(' TE root point: [%.6f %.6f %.6f]\n', te_root(1), te_root(2), te_root(3));
    fprintf(' TE tip  point: [%.6f %.6f %.6f]\n', te_tip(1), te_tip(2), te_tip(3));

    figure('Name', 'Airfoil file branches');
    hold on;
    plot(xy(:,1), xy(:,2), 'k.-', 'LineWidth', 1.0);
    xy_up = details.upper_xy{1};
    xy_lo = details.lower_xy{1};
    plot(xy_up(:,1), xy_up(:,2), '-o', 'LineWidth', 1.5);
    plot(xy_lo(:,1), xy_lo(:,2), '-o', 'LineWidth', 1.5);
    legend('Input / sharpened file', 'Upper: TE->LE', 'Lower: LE->TE', 'Location', 'Best');
    xlabel('x/c'); ylabel('z/c');
    title('Airfoil file split into upper/lower branches');
    axis equal; grid on;

    fprintf('\n=== Summary ===\n');
    fprintf('This workflow reads a real airfoil file, enforces a sharp TE, builds\n');
    fprintf('separate upper/lower wing surfaces, and returns an explicit 3D TE curve.\n');
    fprintf('That TE curve is the right place to attach a wake surface later.\n');
end
