function [xy_upper, xy_lower, info] = splitAirfoilSelig(xy, varargin)
% SPLITAIRFOILSELIG  Split a Selig-format airfoil into upper and lower branches.
%
%   [xy_upper, xy_lower, info] = geom.splitAirfoilSelig(xy)
%   [xy_upper, xy_lower, info] = geom.splitAirfoilSelig(xy, 'CloseTE', true)
%
% Inputs:
%   xy  - [N x 2] airfoil coordinates in Selig order:
%         TE -> upper -> LE -> lower -> TE
%
% Name-value options:
%   'CloseTE'       - logical, force upper/lower TE points to coincide (default true)
%   'TEPoint'       - 'average' | 'first' | 'last'  (default 'average')
%   'RemoveDupLE'   - logical, remove duplicate LE point from lower branch (default true)
%   'Tol'           - tolerance for repeated points / LE detection (default 1e-10)
%
% Outputs:
%   xy_upper  - [Nu x 2] upper branch, ordered TE -> LE
%   xy_lower  - [Nl x 2] lower branch, ordered LE -> TE
%   info      - struct with fields:
%                 .iLE
%                 .te_first
%                 .te_last
%                 .te_used

    pa = inputParser;
    addParameter(pa, 'CloseTE', true, @(x)islogical(x) || isnumeric(x));
    addParameter(pa, 'TEPoint', 'average', @(s)ischar(s) || isstring(s));
    addParameter(pa, 'RemoveDupLE', true, @(x)islogical(x) || isnumeric(x));
    addParameter(pa, 'Tol', 1e-10, @(x)isnumeric(x) && isscalar(x) && x > 0);
    parse(pa, varargin{:});
    opts = pa.Results;

    if size(xy,2) ~= 2
        error('geom.splitAirfoilSelig: xy must be [N x 2].');
    end
    if size(xy,1) < 5
        error('geom.splitAirfoilSelig: need at least 5 points.');
    end

    xy = double(xy);

    xmin = min(xy(:,1));
    i_candidates = find(abs(xy(:,1) - xmin) <= opts.Tol);
    iLE = i_candidates(ceil(numel(i_candidates)/2));

    xy_upper = xy(1:iLE, :);      % TE -> upper -> LE
    xy_lower = xy(iLE:end, :);    % LE -> lower -> TE

    te_first = xy_upper(1,:);
    te_last  = xy_lower(end,:);

    if opts.CloseTE
        switch lower(string(opts.TEPoint))
            case "first"
                te_used = te_first;
            case "last"
                te_used = te_last;
            otherwise
                te_used = 0.5*(te_first + te_last);
        end
        xy_upper(1,:) = te_used;
        xy_lower(end,:) = te_used;
    else
        te_used = 0.5*(te_first + te_last);
    end

    info = struct();
    info.iLE = iLE;
    info.te_first = te_first;
    info.te_last  = te_last;
    info.te_used  = te_used;
end
