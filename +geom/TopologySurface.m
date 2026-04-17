classdef TopologySurface < handle
    properties
        name
        S   % geom.NURBSSurface
    end

    methods
        function obj = TopologySurface(name, S)
            if nargin < 2
                error('TopologySurface requires name and surface.');
            end
            obj.name = char(name);
            obj.S = S;
        end
    end
end