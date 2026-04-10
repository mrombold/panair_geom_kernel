function [S_upper, S_lower, C_te, details] = wingSurfaceUpperLowerFromFile(airfoil_files, spans, chords, sweeps, twists, dihedrals, varargin)
% WINGSURFACEUPPERLOWERFROMFILE  Build separate upper/lower wing surfaces from airfoil files.
%
%   [S_upper, S_lower, C_te, details] = geom.wingSurfaceUpperLowerFromFile( ...
%       airfoil_files, spans, chords, sweeps, twists, dihedrals, ...)
%
% Inputs:
%   airfoil_files - string path to one airfoil file used for every section,
%                   or cell array of per-section file paths
%   spans       - [nsec x 1] spanwise stations (Y)
%   chords      - [nsec x 1] chord lengths
%   sweeps      - [nsec x 1] leading-edge X offsets (optional)
%   twists      - [nsec x 1] twist angles in deg, positive nose-up (optional)
%   dihedrals   - [nsec x 1] Z offsets (optional)
%
% Name-value options:
%   'Degree'          - curve degree for branch fitting (default 3)
%   'NormalizeChord'  - normalize file coordinates to unit chord (default true)
%   'SharpTE'         - force exact trailing-edge coincidence (default true)
%   'TEPoint'         - 'average' | 'first' | 'last' (default 'average')
%   'Method'          - loft method passed to geom.LoftedSurface (default 'chord')
%   'QuarterChordX'   - twist axis x/c location (default 0.25)
%
% Outputs:
%   S_upper, S_lower - separate lofted wing surfaces
%   C_te             - 3D trailing-edge curve through section TE points
%   details          - struct with file metadata, split info, and section curves

    pa = inputParser;
    addParameter(pa, 'Degree', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(pa, 'NormalizeChord', true, @(x)islogical(x) || isnumeric(x));
    addParameter(pa, 'SharpTE', true, @(x)islogical(x) || isnumeric(x));
    addParameter(pa, 'TEPoint', 'average', @(s)ischar(s) || isstring(s));
    addParameter(pa, 'Method', 'chord', @(s)ischar(s) || isstring(s));
    addParameter(pa, 'QuarterChordX', 0.25, @(x)isnumeric(x)&&isscalar(x));
    parse(pa, varargin{:});
    opts = pa.Results;

    nsec = numel(spans);
    spans = spans(:);
    chords = chords(:);

    if nargin < 4 || isempty(sweeps),    sweeps    = zeros(nsec,1); else, sweeps = sweeps(:); end
    if nargin < 5 || isempty(twists),    twists    = zeros(nsec,1); else, twists = twists(:); end
    if nargin < 6 || isempty(dihedrals), dihedrals = zeros(nsec,1); else, dihedrals = dihedrals(:); end

    if ~iscell(airfoil_files)
        f0 = airfoil_files;
        airfoil_files = repmat({f0}, nsec, 1);
    end

    if numel(airfoil_files) ~= nsec
        error('geom.wingSurfaceUpperLowerFromFile: number of airfoil files must match spans.');
    end
    if numel(chords) ~= nsec || numel(sweeps) ~= nsec || numel(twists) ~= nsec || numel(dihedrals) ~= nsec
        error('geom.wingSurfaceUpperLowerFromFile: spans/chords/sweeps/twists/dihedrals must all have nsec entries.');
    end

    p = opts.Degree;

    upper_curves = cell(nsec,1);
    lower_curves = cell(nsec,1);
    file_meta    = cell(nsec,1);
    split_info   = cell(nsec,1);
    upper_xy_all = cell(nsec,1);
    lower_xy_all = cell(nsec,1);
    te_pts       = zeros(nsec, 3);

    for k = 1:nsec
        [xy, meta] = geom.readAirfoilFile(airfoil_files{k}, ...
            'NormalizeChord', opts.NormalizeChord, ...
            'SharpTE', opts.SharpTE, ...
            'TEPoint', opts.TEPoint);

        [xy_up, xy_lo, info] = geom.splitAirfoilSelig(xy, ...
            'CloseTE', opts.SharpTE, 'TEPoint', opts.TEPoint);

        file_meta{k}  = meta;
        split_info{k} = info;
        upper_xy_all{k} = xy_up;
        lower_xy_all{k} = xy_lo;

        P_up = [xy_up(:,1), zeros(size(xy_up,1),1), xy_up(:,2)];
        P_lo = [xy_lo(:,1), zeros(size(xy_lo,1),1), xy_lo(:,2)];

        Cup = geom.LoftedSurface.globalCurveInterp(P_up, p, size(P_up,1));
        Clo = geom.LoftedSurface.globalCurveInterp(P_lo, p, size(P_lo,1));

        chord_k    = chords(k);
        twist_k    = twists(k) * pi/180;
        span_k     = spans(k);
        sweep_k    = sweeps(k);
        dihedral_k = dihedrals(k);
        xqc        = opts.QuarterChordX * chord_k;

        Cup = localTransformCurve(Cup, chord_k, twist_k, xqc, span_k, sweep_k, dihedral_k);
        Clo = localTransformCurve(Clo, chord_k, twist_k, xqc, span_k, sweep_k, dihedral_k);

        upper_curves{k} = Cup;
        lower_curves{k} = Clo;

        % Shared trailing-edge point for this section
        te_local = [info.te_used(1), 0, info.te_used(2)];
        te_local(:,1) = te_local(:,1) * chord_k;
        te_local(:,3) = te_local(:,3) * chord_k;

        te_local(:,1) = te_local(:,1) - xqc;
        x_rot  = te_local(:,1) * cos(twist_k) - te_local(:,3) * sin(twist_k);
        z_rot  = te_local(:,1) * sin(twist_k) + te_local(:,3) * cos(twist_k);
        te_local(:,1) = x_rot + xqc + sweep_k;
        te_local(:,3) = z_rot + dihedral_k;
        te_local(:,2) = span_k;

        te_pts(k,:) = te_local;
    end

    S_upper = geom.LoftedSurface(upper_curves, 'method', char(opts.Method));
    S_lower = geom.LoftedSurface(lower_curves, 'method', char(opts.Method));
    C_te    = geom.LoftedSurface.globalCurveInterp(te_pts, min(3, max(1, size(te_pts,1)-1)), size(te_pts,1));

    details = struct();
    details.upper_curves = upper_curves;
    details.lower_curves = lower_curves;
    details.file_meta = file_meta;
    details.upper_xy = upper_xy_all;
    details.lower_xy = lower_xy_all;
    details.split_info = split_info;
    details.te_pts = te_pts;
    details.options = opts;
end

function C2 = localTransformCurve(C0, chord_k, twist_k, xqc, span_k, sweep_k, dihedral_k)
    P = C0.P;

    P(:,1) = P(:,1) * chord_k;
    P(:,3) = P(:,3) * chord_k;

    P(:,1) = P(:,1) - xqc;
    x_rot  = P(:,1) * cos(twist_k) - P(:,3) * sin(twist_k);
    z_rot  = P(:,1) * sin(twist_k) + P(:,3) * cos(twist_k);
    P(:,1) = x_rot + xqc + sweep_k;
    P(:,3) = z_rot + dihedral_k;

    P(:,2) = span_k;

    C2 = geom.NURBSCurve(P, C0.p, C0.U, C0.W);
end
