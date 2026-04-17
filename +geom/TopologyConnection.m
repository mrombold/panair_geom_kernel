classdef TopologyConnection
    properties
        surfA
        edgeA
        surfB
        edgeB
        tag
    end

    methods
        function obj = TopologyConnection(surfA, edgeA, surfB, edgeB, tag)
            if nargin < 4
                error('TopologyConnection requires surfA, edgeA, surfB, edgeB.');
            end
            if nargin < 5
                tag = '';
            end

            obj.surfA = surfA;
            obj.edgeA = char(edgeA);
            obj.surfB = surfB;
            obj.edgeB = char(edgeB);
            obj.tag   = char(tag);
        end
    end
end