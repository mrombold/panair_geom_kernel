classdef TopologyMesh < handle
    properties
        surfaces = {}
        connections = {}
        meshes = struct()
    end

    methods
        function addSurface(obj, topoSurf)
            if ~isa(topoSurf, 'geom.TopologySurface')
                error('TopologyMesh.addSurface expects geom.TopologySurface.');
            end
            obj.surfaces{end+1} = topoSurf;
        end

        function connect(obj, surfA, edgeA, surfB, edgeB, tag)
            if nargin < 6
                tag = '';
            end
            obj.connections{end+1} = ...
                geom.TopologyConnection(surfA, edgeA, surfB, edgeB, tag);
        end

        function buildMeshes(obj, resolution)
            obj.meshes = struct();

            for k = 1:numel(obj.surfaces)
                s = obj.surfaces{k};

                if ~isfield(resolution, s.name)
                    error('TopologyMesh.buildMeshes: no resolution specified for "%s".', s.name);
                end

                r = resolution.(s.name);
                if numel(r) ~= 2
                    error('TopologyMesh.buildMeshes: resolution for "%s" must be [nu, nv].', s.name);
                end

                obj.meshes.(s.name) = s.S.isoMesh(r(1), r(2));

                fprintf(' Meshed %-10s : [%d x %d]\n', ...
                    s.name, ...
                    size(obj.meshes.(s.name).X,1), ...
                    size(obj.meshes.(s.name).X,2));
            end
        end

        function report = validate(obj, tolerance, verbose)
            if nargin < 2 || isempty(tolerance)
                tolerance = 1e-10;
            end
            if nargin < 3
                verbose = true;
            end

            report = struct();
            report.connections = cell(1, numel(obj.connections));
            report.ok = true;

            if verbose
                fprintf('Topology validation (tol = %.3e)\n', tolerance);
            end

            for k = 1:numel(obj.connections)
                c = obj.connections{k};

                if ~isfield(obj.meshes, c.surfA.name)
                    error('TopologyMesh.validate: missing mesh for surface "%s".', c.surfA.name);
                end
                if ~isfield(obj.meshes, c.surfB.name)
                    error('TopologyMesh.validate: missing mesh for surface "%s".', c.surfB.name);
                end

                mA = obj.meshes.(c.surfA.name);
                mB = obj.meshes.(c.surfB.name);

                PA = geom.getEdge(mA, c.edgeA);
                PB = geom.getEdge(mB, c.edgeB);

                cmp = geom.compareEdgePair(PA, PB);

                item = struct();
                item.surfA = c.surfA.name;
                item.edgeA = c.edgeA;
                item.surfB = c.surfB.name;
                item.edgeB = c.edgeB;
                item.tag   = c.tag;
                item.sameCount = cmp.sameCount;
                item.sameError = cmp.sameError;
                item.revError  = cmp.revError;
                item.bestError = cmp.bestError;
                item.mode      = cmp.mode;
                item.withinTolerance = cmp.sameCount && (cmp.bestError <= tolerance);

                report.connections{k} = item;

                if ~item.withinTolerance
                    report.ok = false;
                end

                if verbose
                    if ~cmp.sameCount
                        fprintf('  %-10s %-3s <-> %-10s %-3s : node count mismatch\n', ...
                            item.surfA, item.edgeA, item.surfB, item.edgeB);
                    else
                        fprintf('  %-10s %-3s <-> %-10s %-3s : best=%.3e mode=%s %s\n', ...
                            item.surfA, item.edgeA, ...
                            item.surfB, item.edgeB, ...
                            item.bestError, item.mode, ...
                            ternaryLocal(item.withinTolerance, '[ok]', '[FAIL]'));
                    end
                end
            end
        end

        function list = neighbors(obj, surfName)
            surfName = char(surfName);
            list = {};

            for k = 1:numel(obj.connections)
                c = obj.connections{k};

                if strcmpi(c.surfA.name, surfName)
                    list(end+1,:) = {c.surfA.name, c.edgeA, c.surfB.name, c.edgeB}; %#ok<AGROW>
                elseif strcmpi(c.surfB.name, surfName)
                    list(end+1,:) = {c.surfB.name, c.edgeB, c.surfA.name, c.edgeA}; %#ok<AGROW>
                end
            end
        end
    end
end

function s = ternaryLocal(cond, a, b)
if cond
    s = a;
else
    s = b;
end
end