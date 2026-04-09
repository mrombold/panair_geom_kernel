classdef MeshWriter
% MESHWRITER  Export isoparametric mesh structs to various formats.
%
%   Converts the mesh struct produced by NURBSSurface.isoMesh() or
%   LoftedSurface.isoMesh() into file formats for downstream tools.
%
%   Static methods:
%     MeshWriter.toWGS(networks, filename)    % PanAir / LaWGS format
%     MeshWriter.toVTK(mesh, filename)        % VTK legacy unstructured grid
%     MeshWriter.toSTL(mesh, filename)        % ASCII STL
%     MeshWriter.toCSV(mesh, filename)        % flat CSV of node coords
%
%   'networks' is a struct array or cell array of mesh structs, each
%   with fields:  .X .Y .Z .name (optional)
%
%   'mesh' is a single mesh struct from isoMesh().
%
% PanAir WGS reference:
%   Carmichael, R.L. & Erickson, L.L. (1981). PAN AIR -- A Higher Order
%   Panel Method for Predicting Subsonic or Supersonic Linear Potential
%   Flows about Arbitrary Configurations. AIAA-81-1255.

    methods (Static)

        % ============================================================== %
        %  PanAir WGS (LaWGS) format
        % ============================================================== %

        function toWGS(networks, filename)
        % TOWGS  Write PanAir LaWGS format grid file.
        %
        %   toWGS(networks, filename)
        %
        %   networks  - struct array or cell array of mesh structs, each with:
        %               .X, .Y, .Z  [nRows x nCols] coordinate arrays
        %               .name       string label (optional, default 'Network_k')
        %               .iSymX      x-symmetry: 0=none, 1=reflect in XZ plane (y→-y)
        %               .iSymZ      z-symmetry: 0=none, 1=reflect in XY plane (z→-z)
        %   filename  - output path, e.g. 'wing.wgs'
        %
        %   WGS format (LaWGS):
        %     Line 1: number of networks
        %     Per network:
        %       Name(1-32)  iSymX  iSymZ  nRows  nCols
        %       x(row1,col1) x(row1,col2) ... [5 per line, E12.5]
        %       x(row2,...)
        %       ...
        %       y(row1,col1) ...
        %       z(row1,col1) ...

            if ~iscell(networks)
                % Convert struct array to cell
                tmp = networks;
                networks = cell(1, numel(tmp));
                for k = 1:numel(tmp)
                    networks{k} = tmp(k);
                end
            end
            n_net = numel(networks);

            fid = fopen(filename, 'w');
            if fid < 0
                error('MeshWriter.toWGS: cannot open "%s" for writing.', filename);
            end

            % Line 1: number of networks
            fprintf(fid, '%6d\n', n_net);

            for k = 1:n_net
                net = networks{k};

                % Defaults
                if ~isfield(net, 'name'),  net.name  = sprintf('Network_%02d', k); end
                if ~isfield(net, 'iSymX'), net.iSymX = 0; end
                if ~isfield(net, 'iSymZ'), net.iSymZ = 0; end

                nRows = size(net.X, 1);
                nCols = size(net.X, 2);

                % Header line: name padded to 32 chars, symmetry, dimensions
                name_padded = sprintf('%-32s', net.name(1:min(end,32)));
                fprintf(fid, ' %s %2d %2d %4d %4d\n', ...
                        name_padded, net.iSymX, net.iSymZ, nRows, nCols);

                % Write coordinate blocks: X then Y then Z, row by row, 5 per line
                for coord_block = 1:3
                    switch coord_block
                        case 1, data = net.X;
                        case 2, data = net.Y;
                        case 3, data = net.Z;
                    end
                    n = 0;
                    for row = 1:nRows
                        for col = 1:nCols
                            fprintf(fid, ' %12.5E', data(row,col));
                            n = n + 1;
                            if mod(n, 5) == 0
                                fprintf(fid, '\n');
                            end
                        end
                    end
                    if mod(n, 5) ~= 0
                        fprintf(fid, '\n');
                    end
                end
            end

            fclose(fid);
            fprintf('MeshWriter: wrote %d network(s) to "%s"\n', n_net, filename);
        end

        % ============================================================== %
        %  VTK Legacy Unstructured Grid
        % ============================================================== %

        function toVTK(mesh_or_networks, filename, varargin)
        % TOVTK  Write VTK legacy format unstructured grid.
        %
        %   toVTK(mesh, filename)
        %   toVTK(mesh, filename, 'IncludeNormals', true)
        %
        %   Writes quad cells for direct loading in ParaView.
        %   Multiple meshes can be passed as a cell array - they are
        %   merged into one unstructured grid.
        %
        %   Output: VTK legacy ASCII unstructured grid with:
        %     - POINTS: all node coordinates
        %     - CELLS:  quad connectivity
        %     - POINT_DATA: surface normals (if present and requested)

            pa = inputParser;
            addParameter(pa, 'IncludeNormals', true);
            addParameter(pa, 'BinaryMode',     false);
            parse(pa, varargin{:});
            opts = pa.Results;

            % Normalize to cell array
            if ~iscell(mesh_or_networks)
                meshes = {mesh_or_networks};
            else
                meshes = mesh_or_networks;
            end

            % Collect all points and cells across meshes
            all_pts  = zeros(0, 3);
            all_quads= zeros(0, 4);
            all_nrms = zeros(0, 3);
            has_normals = true;
            offset = 0;

            for k = 1:numel(meshes)
                m = meshes{k};
                nu = m.nu; nv = m.nv;

                % Flatten node grid to list
                pts = [m.X(:), m.Y(:), m.Z(:)];
                all_pts = [all_pts; pts]; %#ok<AGROW>

                if has_normals && isfield(m, 'normals')
                    nrm = reshape(m.normals, [], 3);
                    all_nrms = [all_nrms; nrm]; %#ok<AGROW>
                else
                    has_normals = false;
                end

                % Quad connectivity (0-based for VTK)
                for i = 1:nu-1
                    for j = 1:nv-1
                        n1 = offset + (i-1)*nv + (j-1);
                        n2 = offset + (i  )*nv + (j-1);
                        n3 = offset + (i  )*nv + (j  );
                        n4 = offset + (i-1)*nv + (j  );
                        all_quads(end+1,:) = [n1, n2, n3, n4]; %#ok<AGROW>
                    end
                end
                offset = offset + nu*nv;
            end

            n_pts   = size(all_pts,  1);
            n_cells = size(all_quads,1);

            fid = fopen(filename, 'w');
            if fid < 0
                error('MeshWriter.toVTK: cannot open "%s" for writing.', filename);
            end

            % VTK header
            fprintf(fid, '# vtk DataFile Version 3.0\n');
            fprintf(fid, 'Aircraft Geometry Kernel surface mesh\n');
            fprintf(fid, 'ASCII\n');
            fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n\n');

            % Points
            fprintf(fid, 'POINTS %d float\n', n_pts);
            fprintf(fid, '%f %f %f\n', all_pts');
            fprintf(fid, '\n');

            % Cells (VTK_QUAD = type 9, 4 nodes each)
            fprintf(fid, 'CELLS %d %d\n', n_cells, n_cells*5);
            for q = 1:n_cells
                fprintf(fid, '4 %d %d %d %d\n', all_quads(q,:));
            end
            fprintf(fid, '\n');

            % Cell types
            fprintf(fid, 'CELL_TYPES %d\n', n_cells);
            fprintf(fid, repmat('9\n', 1, n_cells));
            fprintf(fid, '\n');

            % Point data: normals
            if has_normals && opts.IncludeNormals && ~isempty(all_nrms)
                fprintf(fid, 'POINT_DATA %d\n', n_pts);
                fprintf(fid, 'NORMALS SurfaceNormals float\n');
                fprintf(fid, '%f %f %f\n', all_nrms');
            end

            fclose(fid);
            fprintf('MeshWriter: wrote %d pts, %d quads to "%s"\n', ...
                    n_pts, n_cells, filename);
        end

        % ============================================================== %
        %  STL
        % ============================================================== %

        function toSTL(mesh_or_networks, filename)
        % TOSTL  Write ASCII STL file (triangulated quads).
        %
        %   Each quad is split into 2 triangles.

            if ~iscell(mesh_or_networks)
                meshes = {mesh_or_networks};
            else
                meshes = mesh_or_networks;
            end

            fid = fopen(filename, 'w');
            if fid < 0
                error('MeshWriter.toSTL: cannot open "%s" for writing.', filename);
            end

            solid_name = 'AircraftGeom';
            fprintf(fid, 'solid %s\n', solid_name);

            for k = 1:numel(meshes)
                m = meshes{k};
                nu = m.nu; nv = m.nv;

                for i = 1:nu-1
                    for j = 1:nv-1
                        % 4 corner points
                        V = zeros(4,3);
                        V(1,:) = [m.X(i,j),   m.Y(i,j),   m.Z(i,j)  ];
                        V(2,:) = [m.X(i+1,j), m.Y(i+1,j), m.Z(i+1,j)];
                        V(3,:) = [m.X(i+1,j+1),m.Y(i+1,j+1),m.Z(i+1,j+1)];
                        V(4,:) = [m.X(i,j+1), m.Y(i,j+1), m.Z(i,j+1)];

                        % Triangle 1: V1 V2 V3
                        n1 = geom.MeshWriter.triNormal(V(1,:),V(2,:),V(3,:));
                        fprintf(fid, '  facet normal %f %f %f\n', n1);
                        fprintf(fid, '    outer loop\n');
                        fprintf(fid, '      vertex %f %f %f\n', V(1,:));
                        fprintf(fid, '      vertex %f %f %f\n', V(2,:));
                        fprintf(fid, '      vertex %f %f %f\n', V(3,:));
                        fprintf(fid, '    endloop\n');
                        fprintf(fid, '  endfacet\n');

                        % Triangle 2: V1 V3 V4
                        n2 = geom.MeshWriter.triNormal(V(1,:),V(3,:),V(4,:));
                        fprintf(fid, '  facet normal %f %f %f\n', n2);
                        fprintf(fid, '    outer loop\n');
                        fprintf(fid, '      vertex %f %f %f\n', V(1,:));
                        fprintf(fid, '      vertex %f %f %f\n', V(3,:));
                        fprintf(fid, '      vertex %f %f %f\n', V(4,:));
                        fprintf(fid, '    endloop\n');
                        fprintf(fid, '  endfacet\n');
                    end
                end
            end

            fprintf(fid, 'endsolid %s\n', solid_name);
            fclose(fid);

            n_tris = sum(cellfun(@(m) 2*(m.nu-1)*(m.nv-1), meshes));
            fprintf('MeshWriter: wrote %d triangles to "%s"\n', n_tris, filename);
        end

        % ============================================================== %
        %  CSV
        % ============================================================== %

        function toCSV(mesh, filename)
        % TOCSV  Write mesh nodes as CSV: x,y,z,nx,ny,nz,u,v

            fid = fopen(filename, 'w');
            if fid < 0
                error('MeshWriter.toCSV: cannot open "%s".', filename);
            end

            has_normals = isfield(mesh, 'normals');
            has_params  = isfield(mesh, 'u');

            if has_normals && has_params
                fprintf(fid, 'x,y,z,nx,ny,nz,u,v\n');
            elseif has_normals
                fprintf(fid, 'x,y,z,nx,ny,nz\n');
            else
                fprintf(fid, 'x,y,z\n');
            end

            nu = mesh.nu; nv = mesh.nv;
            for i = 1:nu
                for j = 1:nv
                    if has_normals && has_params
                        fprintf(fid, '%.8f,%.8f,%.8f,%.6f,%.6f,%.6f,%.6f,%.6f\n', ...
                                mesh.X(i,j), mesh.Y(i,j), mesh.Z(i,j), ...
                                mesh.normals(i,j,1), mesh.normals(i,j,2), mesh.normals(i,j,3), ...
                                mesh.u(i,j), mesh.v(i,j));
                    elseif has_normals
                        fprintf(fid, '%.8f,%.8f,%.8f,%.6f,%.6f,%.6f\n', ...
                                mesh.X(i,j), mesh.Y(i,j), mesh.Z(i,j), ...
                                mesh.normals(i,j,1), mesh.normals(i,j,2), mesh.normals(i,j,3));
                    else
                        fprintf(fid, '%.8f,%.8f,%.8f\n', ...
                                mesh.X(i,j), mesh.Y(i,j), mesh.Z(i,j));
                    end
                end
            end

            fclose(fid);
            fprintf('MeshWriter: wrote %dx%d = %d nodes to "%s"\n', ...
                    nu, nv, nu*nv, filename);
        end

        % ============================================================== %
        %  Helper: build WGS network struct from isoMesh output
        % ============================================================== %

        function net = meshToNetwork(mesh, name, iSymX, iSymZ)
        % MESHTONETWORK  Package an isoMesh struct as a WGS network.
        %
        %   net = meshToNetwork(mesh, name, iSymX, iSymZ)
        %
        %   PanAir expects the mesh ordered so that the normal computed
        %   from (row2-row1) x (col2-col1) points into the flow domain
        %   (typically outward from the body).
        %
        %   iSymX: 0=no sym, 1=reflect y->-y (XZ plane symmetry)
        %   iSymZ: 0=no sym, 1=reflect z->-z (XY plane symmetry)

            if nargin < 2, name  = 'Network'; end
            if nargin < 3, iSymX = 0; end
            if nargin < 4, iSymZ = 0; end

            net.X     = mesh.X;
            net.Y     = mesh.Y;
            net.Z     = mesh.Z;
            net.name  = name;
            net.iSymX = iSymX;
            net.iSymZ = iSymZ;
        end

        function checkNormals(mesh)
        % CHECKNORMALS  Report normal direction statistics (diagnostic).
        %   Useful for verifying PanAir outward normal convention.

            nrm = reshape(mesh.normals, [], 3);
            fprintf('Normal direction statistics:\n');
            fprintf('  Mean normal: [%.3f %.3f %.3f]\n', mean(nrm));
            fprintf('  Z-component: min=%.3f  max=%.3f\n', ...
                    min(nrm(:,3)), max(nrm(:,3)));

            % For a wing upper surface, Z should be mostly positive
            z_pos_frac = sum(nrm(:,3) > 0) / size(nrm,1);
            fprintf('  Fraction pointing +Z: %.1f%%\n', 100*z_pos_frac);
            if z_pos_frac < 0.5
                fprintf('  WARNING: normals may need flipping (call mesh.Z = flipud(mesh.Z) etc.)\n');
            end
        end

    end  % static methods

    methods (Static, Access = private)

        function n = triNormal(V1, V2, V3)
            e1 = V2 - V1; e2 = V3 - V1;
            n  = cross(e1, e2);
            nrm = norm(n);
            if nrm > eps, n = n/nrm; end
        end

    end

end  % classdef MeshWriter
