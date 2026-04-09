classdef PlotHelpers
% PLOTHELPERS  Small plotting helpers for demo scripts.

    methods (Static)
        function meshStruct(M, faceColor)
        % MESHSTRUCT  Plot a mesh struct produced by isoMesh().
            surf(M.X, M.Y, M.Z, ...
                'FaceColor', faceColor, ...
                'FaceAlpha', 0.55, ...
                'EdgeColor', [0.15 0.15 0.15], ...
                'EdgeAlpha', 0.55);
        end
    end
end
