function [xy, meta] = readAirfoilFile(filename, varargin)
% READAIRFOILFILE  Read a standard airfoil coordinate file.
%
%   [xy, meta] = geom.readAirfoilFile(filename)
%   [xy, meta] = geom.readAirfoilFile(filename, 'NormalizeChord', true, 'SharpTE', true)
%
% Inputs:
%   filename - path to an airfoil coordinate file (typically .dat)
%
% Name-value options:
%   'NormalizeChord' - logical, normalize x to [0,1] using min/max x (default true)
%   'SharpTE'        - logical, force first and last points to share one TE point (default false)
%   'TEPoint'        - 'average' | 'first' | 'last' (default 'average')
%   'Tol'            - parsing/duplicate tolerance (default 1e-10)
%
% Outputs:
%   xy   - [N x 2] numeric coordinates
%   meta - struct with fields:
%            .filename
%            .header
%            .npts
%            .xmin .xmax
%            .te_first .te_last .te_used
%
% Notes:
%   - This reader skips blank lines and non-numeric header lines.
%   - It preserves input ordering.
%   - Most curated airfoil files are already in Selig order:
%       TE -> upper -> LE -> lower -> TE

    pa = inputParser;
    addParameter(pa, 'NormalizeChord', true, @(x)islogical(x) || isnumeric(x));
    addParameter(pa, 'SharpTE', false, @(x)islogical(x) || isnumeric(x));
    addParameter(pa, 'TEPoint', 'average', @(s)ischar(s) || isstring(s));
    addParameter(pa, 'Tol', 1e-10, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    parse(pa, varargin{:});
    opts = pa.Results;

    if ~(ischar(filename) || isstring(filename))
        error('geom.readAirfoilFile: filename must be a string.');
    end
    filename = char(filename);

    fid = fopen(filename, 'r');
    if fid < 0
        error('geom.readAirfoilFile: could not open file "%s".', filename);
    end
    cleanup = onCleanup(@() fclose(fid));

    lines = {};
    while true
        tline = fgetl(fid);
        if ~ischar(tline), break; end
        lines{end+1,1} = tline; %#ok<AGROW>
    end

    xy = zeros(0,2);
    header = '';

    for k = 1:numel(lines)
        s = strtrim(lines{k});
        if isempty(s)
            continue;
        end

        vals = sscanf(s, '%f %f');
        if numel(vals) >= 2
            xy(end+1,:) = vals(1:2).'; %#ok<AGROW>
        else
            if isempty(header)
                header = s;
            end
        end
    end

    if size(xy,1) < 5
        error('geom.readAirfoilFile: file "%s" did not contain enough coordinate rows.', filename);
    end

    % Remove consecutive duplicate points
    keep = true(size(xy,1),1);
    for k = 2:size(xy,1)
        if norm(xy(k,:) - xy(k-1,:)) <= opts.Tol
            keep(k) = false;
        end
    end
    xy = xy(keep,:);

    xmin = min(xy(:,1));
    xmax = max(xy(:,1));

    if opts.NormalizeChord
        if abs(xmax - xmin) < opts.Tol
            error('geom.readAirfoilFile: degenerate x range.');
        end
        xy(:,1) = (xy(:,1) - xmin) / (xmax - xmin);
        xy(:,2) = xy(:,2) / (xmax - xmin);
        xmin = min(xy(:,1));
        xmax = max(xy(:,1));
    end

    te_first = xy(1,:);
    te_last  = xy(end,:);

    if opts.SharpTE
        switch lower(string(opts.TEPoint))
            case "first"
                te_used = te_first;
            case "last"
                te_used = te_last;
            otherwise
                te_used = 0.5 * (te_first + te_last);
        end
        xy(1,:)   = te_used;
        xy(end,:) = te_used;
    else
        te_used = 0.5 * (te_first + te_last);
    end

    meta = struct();
    meta.filename = filename;
    meta.header = header;
    meta.npts = size(xy,1);
    meta.xmin = xmin;
    meta.xmax = xmax;
    meta.te_first = te_first;
    meta.te_last  = te_last;
    meta.te_used  = te_used;
end
