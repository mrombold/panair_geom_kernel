function [S_upper, S_lower, details] = wingSurfaceUpperLower(airfoil_xy, spans, chords, sweeps, twists, dihedrals, varargin)
% WINGSURFACEUPPERLOWER  Build separate upper/lower wing surfaces from Selig airfoil data.
%
%   [S_upper, S_lower, details] = geom.wingSurfaceUpperLower( ...
%       airfoil_xy, spans, chords, sweeps, twists, dihedrals, ...)
%
% Inputs:
%   airfoil_xy  - [N x 2] Selig-format airfoil coordinates for one section,
%                 or cell array of per-section Selig airfoils
%   spans       - [nsec x 1] spanwise stations (Y)
%   chords      - [nsec x 1] chord lengths
%   sweeps      - [nsec x 1] leading-edge X offsets (optional)
%   twists      - [nsec x 1] twist angles in deg, positive nose-up (optional)
%   dihedrals   - [nsec x 1] Z offsets (optional)
%
% Name-value options:
%   'Degree'          - curve degree for branch fitting (default 3)
%   'NCP'             - control points per section branch (default 18)
%   'CloseTE'         - force upper/lower TE coincidence (default true)
%   'TEPoint'         - 'average' | 'first' | 'last' (default 'average')
%   'Method'          - loft method passed to geom.LoftedSurface (default 'chord')
%   'QuarterChordX'   - twist axis x/c location (default 0.25)
%
% Outputs:
%   S_upper, S_lower  - lofted upper and lower surfaces
%   details           - struct with section curves and split metadata

    pa = inputParser;
    addParameter(pa, 'Degree', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(pa, 'NCP', 18, @(x)isnumeric(x)&&isscalar(x)&&x>=4);
    addParameter(pa, 'CloseTE', true, @(x)islogical(x) || isnumeric(x));
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

    if ~iscell(airfoil_xy)
        af_base = airfoil_xy;
        airfoil_xy = repmat({af_base}, nsec, 1);
    end

    if numel(airfoil_xy) ~= nsec
        error('geom.wingSurfaceUpperLower: number of airfoil sections must match spans.');
    end
    if numel(chords) ~= nsec || numel(sweeps) ~= nsec || numel(twists) ~= nsec || numel(dihedrals) ~= nsec
        error('geom.wingSurfaceUpperLower: spans/chords/sweeps/twists/dihedrals must all have nsec entries.');
    end

    p = opts.Degree;
    n_cp = opts.NCP;

    upper_curves = cell(nsec,1);
    lower_curves = cell(nsec,1);
    split_info   = cell(nsec,1);
    upper_xy_all = cell(nsec,1);
    lower_xy_all = cell(nsec,1);

    for k = 1:nsec
        xy = airfoil_xy{k};
        [xy_up, xy_lo, info] = geom.splitAirfoilSelig(xy, ...
            'CloseTE', opts.CloseTE, 'TEPoint', opts.TEPoint);

        upper_xy_all{k} = xy_up;
        lower_xy_all{k} = xy_lo;
        split_info{k} = info;

        P_up = [xy_up(:,1), zeros(size(xy_up,1),1), xy_up(:,2)];
        P_lo = [xy_lo(:,1), zeros(size(xy_lo,1),1), xy_lo(:,2)];

        Cup = geom.LoftedSurface.globalCurveInterp(P_up, p, min(n_cp, size(P_up,1)));
        Clo = geom.LoftedSurface.globalCurveInterp(P_lo, p, min(n_cp, size(P_lo,1)));

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
    end

    S_upper = geom.LoftedSurface(upper_curves, 'method', char(opts.Method));
    S_lower = geom.LoftedSurface(lower_curves, 'method', char(opts.Method));

    details = struct();
    details.upper_curves = upper_curves;
    details.lower_curves = lower_curves;
    details.upper_xy = upper_xy_all;
    details.lower_xy = lower_xy_all;
    details.split_info = split_info;
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
