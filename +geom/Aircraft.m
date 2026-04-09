classdef Aircraft
% AIRCRAFT  Factory methods for common aircraft surface geometry.
%
%   Provides convenience constructors for wings, fuselages, and other
%   aircraft-specific primitives, returning geom.NURBSSurface or
%   geom.LoftedSurface objects ready for meshing and analysis.
%
%   Static factory methods:
%
%   Curves:
%     C = geom.Aircraft.conicArc(P0, P2, T0, T2, rho)
%         Roy Liming conic arc through endpoints with tangent directions.
%
%     C = geom.Aircraft.airfoilFromCoords(xy, p)
%         NURBS curve fit to airfoil coordinate array [x, y].
%
%   Surfaces:
%     S = geom.Aircraft.ruledSurface(C1, C2, n)
%         Ruled surface between two compatible NURBSCurves.
%
%     S = geom.Aircraft.wingSurface(sections, spans, chords, twists)
%         Trapezoidal wing from section airfoils.
%
%     S = geom.Aircraft.revolutionSurface(profile, axis, angles)
%         Surface of revolution (e.g., axisymmetric fuselage).
%
%     S = geom.Aircraft.fuselageFromSections(sections)
%         Fuselage OML from cross-section NURBSCurves.

    methods (Static)

        % ============================================================== %
        %  CONIC ARC  (Roy Liming / NURBS weight method)
        % ============================================================== %

        function C = conicArc(P0, P2, T0, T2, rho)
        % CONICARC  Rational quadratic Bezier conic arc.
        %
        %   C = conicArc(P0, P2, T0, T2, rho)
        %
        %   P0, P2  - start and end points  [1x3]
        %   T0, T2  - tangent directions at P0, P2 (need not be unit)
        %   rho     - shoulder parameter (0 < rho < 1)
        %             rho = 0.5 -> parabola
        %             rho > 0.5 -> ellipse
        %             rho < 0.5 -> hyperbola (rarely used)
        %
        %   Returns a degree-2 NURBS curve (rational quadratic Bezier).
        %   This is the Roy Liming method used in classic aircraft lofting.

            P0 = P0(:)'; P2 = P2(:)'; T0 = T0(:)'; T2 = T2(:)';
            T0 = T0/norm(T0); T2 = T2/norm(T2);

            % Intersection of tangent lines to find P1 (shoulder point)
            % Solve:  P0 + s*T0 = P2 + t*T2
            % [T0, -T2] * [s; t] = P2 - P0
            A = [T0(1), -T2(1); T0(2), -T2(2)];
            b = (P2(1:2) - P0(1:2))';
            if abs(det(A)) < 1e-12
                % Use 3D least-squares
                A3 = [T0(:), -T2(:)];
                st = A3 \ (P2-P0)';
            else
                st = A \ b;
            end
            P1 = P0 + st(1)*T0;

            % Weight for middle control point (shoulder parameter -> weight)
            % w1 = rho / (1 - rho)
            w1 = rho / (1 - rho);

            P  = [P0; P1; P2];
            W  = [1; w1; 1];
            U  = [0 0 0 1 1 1];  % clamped quadratic
            C  = geom.NURBSCurve(P, 2, U, W);
        end

        function C = conicArcThroughPoint(P0, P2, T0, T2, Q)
            % CONICARCTHROUGHPOINT
            % Rational quadratic Bezier conic through:
            %   - endpoints P0, P2
            %   - endpoint tangents T0, T2
            %   - pass-through point Q
            
                P0 = P0(:)'; P2 = P2(:)'; T0 = T0(:)'; T2 = T2(:)'; Q = Q(:)';
            
                T0 = T0 / norm(T0);
                T2 = T2 / norm(T2);
            
                % Tangent-line intersection V
                A = [T0(:), -T2(:)];
                st = A \ (P2 - P0)';
                V = P0 + st(1) * T0;
            
                % Solve barycentric coordinates of Q in triangle (P0, V, P2)
                M = [P0(:), V(:), P2(:);
                     1,     1,    1];
                rhs = [Q(:); 1];
                lam = M \ rhs;
            
                l0 = lam(1);
                l1 = lam(2);
                l2 = lam(3);
            
                % Validity checks
                if abs(l0*l2) < 1e-14
                    error('Pass-through point too close to an endpoint or degenerate configuration.');
                end
            
                % Middle weight from Liming/projective conic relation
                w1 = l1 / (2 * sqrt(l0 * l2));
            
                P = [P0; V; P2];
                W = [1; w1; 1];
                U = [0 0 0 1 1 1];
            
                C = geom.NURBSCurve(P, 2, U, W);
            end

        function C = conicSection(P_top, P_bot, P_fwd, P_aft, rho_fw, rho_bk, rho_lr)
        % CONICSECTION  Full closed conic cross-section.
        %
        %   Assembles a closed fuselage-style cross section from four
        %   quadrant conic arcs.  Returns a composite NURBSCurve.
        %   For now returns the 4 arcs as a cell array.
        %
        %   This approximates the classic Liming conic fuselage section.

            % Default rho
            if nargin < 5, rho_fw = 0.5; end
            if nargin < 6, rho_bk = 0.5; end
            if nargin < 7, rho_lr = 0.5; end

            % Tangents at cardinal points (tangent = perpendicular to radius)
            % Top/bottom tangent in Y->X direction, Fwd/Aft tangent X->Y
            T_top = [1 0 0]; T_bot = [-1 0 0];
            T_fwd = [0  1 0]; T_aft = [0 -1 0];

            % Four arcs: top-fwd, fwd-bot, bot-aft, aft-top
            arcs{1} = geom.Aircraft.conicArc(P_top, P_fwd, T_top,  T_fwd,  rho_fw);
            arcs{2} = geom.Aircraft.conicArc(P_fwd, P_bot, T_fwd,  T_bot,  rho_lr);
            arcs{3} = geom.Aircraft.conicArc(P_bot, P_aft, T_bot,  T_aft,  rho_bk);
            arcs{4} = geom.Aircraft.conicArc(P_aft, P_top, T_aft,  T_top,  rho_lr);

            C = arcs;  % return cell; use CompositeCurve to join
        end

        % ============================================================== %
        %  AIRFOIL FROM COORDINATES
        % ============================================================== %

        function C = airfoilFromCoords(xy, p, n_cp)
        % AIRFOILFROMCOORDS  Fit a NURBS curve to airfoil coordinate data.
        %
        %   C = airfoilFromCoords(xy, p, n_cp)
        %
        %   xy    - [N x 2] airfoil coordinates (Selig format: TE->upper->LE->lower->TE)
        %   p     - degree (default 3)
        %   n_cp  - number of control points (default: match data = exact interp)
        %
        %   Returns geom.NURBSCurve representing the closed airfoil contour.

            if nargin < 2, p = 3; end

            % Ensure closed
            if norm(xy(1,:) - xy(end,:)) > 1e-6
                xy = [xy; xy(1,:)];
            end

            % Default: exact interpolation through every data point
            if nargin < 3 || isempty(n_cp)
                n_cp = size(xy, 1);
            end

            pts3 = [xy(:,1), zeros(size(xy,1), 1), xy(:,2)];

            C = geom.LoftedSurface.globalCurveInterp(pts3, p, n_cp);
        end

        % ============================================================== %
        %  RULED SURFACE
        % ============================================================== %

        function S = ruledSurface(C1, C2, n_pts)
        % RULEDSURFACE  Linear ruled surface between two NURBSCurves.
        %
        %   S = ruledSurface(C1, C2)
        %
        %   C1, C2 must have compatible parameterization.
        %   Returns a geom.NURBSSurface (degree 1 in v).

            if nargin < 3, n_pts = []; end

            % Harmonize if needed
            cells = geom.LoftedSurface.harmonizeCurves({C1, C2});
            C1 = cells{1}; C2 = cells{2};

            n   = C1.n;
            p   = C1.p;
            U   = C1.U;
            V   = [0 0 1 1];  % linear in v

            P = zeros(n+1, 2, 3);
            W = zeros(n+1, 2);

            P(:,1,:) = C1.P;
            P(:,2,:) = C2.P;
            W(:,1)   = C1.W;
            W(:,2)   = C2.W;

            S = geom.NURBSSurface(P, p, 1, U, V, W);
        end

        % ============================================================== %
        %  WING SURFACE
        % ============================================================== %

        function S = wingSurface(airfoil_xy, spans, chords, sweeps, twists, dihedrals)
        % WINGSURFACE  Build a lofted wing surface from section airfoil data.
        %
        %   S = wingSurface(airfoil_xy, spans, chords, sweeps, twists, dihedrals)
        %
        %   airfoil_xy  - [N x 2] base airfoil coordinates (normalized, unit chord)
        %                 OR cell array of per-section {[N x 2], ...}
        %   spans       - [nsec x 1] spanwise stations (Y-axis)
        %   chords      - [nsec x 1] chord lengths at each station
        %   sweeps      - [nsec x 1] leading-edge sweep (X offset from root LE)
        %   twists      - [nsec x 1] twist angles (deg, positive nose up)
        %   dihedrals   - [nsec x 1] dihedral offset (Z offset from root)
        %
        %   Returns geom.LoftedSurface.

            nsec = numel(spans);
            if nargin < 4 || isempty(sweeps),    sweeps    = zeros(nsec,1); end
            if nargin < 5 || isempty(twists),    twists    = zeros(nsec,1); end
            if nargin < 6 || isempty(dihedrals), dihedrals = zeros(nsec,1); end

            % Use same airfoil for all sections if single array given
            if ~iscell(airfoil_xy)
                af_base = airfoil_xy;
                airfoil_xy = repmat({af_base}, nsec, 1);
            end

            n_cp = 24;   % control points per section
            p    =  3;   % degree

            section_curves = cell(nsec, 1);

            for k = 1:nsec
                % Build normalized airfoil curve
                xy = airfoil_xy{k};
                C0 = geom.Aircraft.airfoilFromCoords(xy, p, n_cp);

                % Scale by chord, apply twist, translate to span station
                chord_k  = chords(k);
                twist_k  = twists(k) * pi/180;
                span_k   = spans(k);
                sweep_k  = sweeps(k);
                dihedral_k = dihedrals(k);

                % Build affine transform: scale -> twist -> translate
                P = C0.P;

                % Scale (chord)
                P(:,1) = P(:,1) * chord_k;
                P(:,3) = P(:,3) * chord_k;

                % Twist about quarter-chord (Y-axis rotation in XZ plane)
                xqc    = 0.25 * chord_k;
                P(:,1) = P(:,1) - xqc;
                x_rot  = P(:,1)*cos(twist_k) - P(:,3)*sin(twist_k);
                z_rot  = P(:,1)*sin(twist_k) + P(:,3)*cos(twist_k);
                P(:,1) = x_rot + xqc + sweep_k;
                P(:,3) = z_rot + dihedral_k;

                % Span station (Y)
                P(:,2) = span_k;

                section_curves{k} = geom.NURBSCurve(P, p, C0.U, C0.W);
            end

            S = geom.LoftedSurface(section_curves, 'method', 'chord');
        end

        % ============================================================== %
        %  SURFACE OF REVOLUTION
        % ============================================================== %

        function S = revolutionSurface(profile_pts, axis_origin, axis_dir, theta_start, theta_end)
        % REVOLUTIONSURFACE  Surface of revolution from a profile curve.
        %
        %   S = revolutionSurface(profile_pts, axis_origin, axis_dir, theta_start, theta_end)
        %
        %   profile_pts  - [N x 3] points defining the profile (meridian)
        %   axis_origin  - [1 x 3] point on axis of revolution
        %   axis_dir     - [1 x 3] axis direction (unit vector)
        %   theta_start  - start angle (deg, default 0)
        %   theta_end    - end angle   (deg, default 360)
        %
        %   Uses the exact NURBS circle construction (rational quadratics).
        %   Returns geom.NURBSSurface.

            if nargin < 4, theta_start = 0;   end
            if nargin < 5, theta_end   = 360; end

            theta  = (theta_end - theta_start) * pi/180;
            axis_dir = axis_dir(:)' / norm(axis_dir);
            axis_origin = axis_origin(:)';

            npts  = size(profile_pts, 1);

            % Build NURBS circle for each profile point
            % Assemble into surface control net
            [P_arcs, W_arcs, V_circ, n_circ] = ...
                geom.Aircraft.buildRevolutionNet(profile_pts, axis_origin, axis_dir, ...
                                                  theta_start*pi/180, theta);

            % U knot: chord-length from profile
            U = geom.BasisFunctions.ChordLengthKnotVector(profile_pts, 2);

            S = geom.NURBSSurface(P_arcs, 2, 2, U, V_circ, W_arcs);
        end

        % ============================================================== %
        %  FUSELAGE FROM SECTIONS
        % ============================================================== %

        function S = fuselageFromSections(sections, x_stations)
        % FUSELAGEFROMSECTIONS  Loft fuselage OML through cross-section curves.
        %
        %   sections    - cell array of NURBSCurve cross sections
        %   x_stations  - [nsec x 1] X (axial) positions of each section
        %
        %   Returns geom.LoftedSurface.

            nsec = numel(sections);
            if nargin < 2
                x_stations = linspace(0, 1, nsec)';
            end
            x_stations = x_stations(:);

            % Ensure all sections lie at their x_station
            for k = 1:nsec
                P = sections{k}.P;
                P(:,1) = x_stations(k);  % force x coordinate
                sections{k} = geom.NURBSCurve(P, sections{k}.p, sections{k}.U, sections{k}.W);
            end

            % Sort by x station
            [~, idx] = sort(x_stations);
            sections = sections(idx);

            S = geom.LoftedSurface(sections, 'method', 'chord');
        end

        % ============================================================== %
        %  SYMMETRIC WING  (mirror half-wing across XZ plane)
        % ============================================================== %

        function [S_left, S_right] = symmetricWing(S_half)
        % SYMMETRICWING  Mirror a half-wing surface to create full-span wing.
        %
        %   [S_left, S_right] = symmetricWing(S_half)
        %
        %   S_half   - LoftedSurface or NURBSSurface of the starboard (right) wing,
        %              built from span Y >= 0.
        %
        %   Returns S_right (original) and S_left (mirror in XZ plane, Y -> -Y).
        %   Use with PanAir iSymX=1 flag to avoid meshing both sides.

            S_right = S_half;

            % Mirror: flip Y control points
            if isa(S_half, 'geom.LoftedSurface')
                surf = S_half.surface;
            else
                surf = S_half;
            end

            P2 = surf.P;
            P2(:,:,2) = -P2(:,:,2);   % negate Y

            % Flip u-direction so normals still point outward
            P2 = flipud(P2);
            W2 = flipud(surf.W);
            U2 = 1 - fliplr(surf.U);

            S_left_surf = geom.NURBSSurface(P2, surf.p, surf.q, U2, surf.V, W2);
            S_left = S_left_surf;
        end

        % ============================================================== %
        %  NACA 5-SERIES AIRFOIL
        % ============================================================== %

        function [x, y] = naca5Coords(des, N)
        % NACA5COORDS  NACA 5-digit airfoil coordinates.
        %
        %   [x, y] = naca5Coords(des, N)
        %
        %   des  - 5-digit designator as integer (e.g. 23012) or string '23012'
        %   N    - points per surface (default 60)
        %
        %   Returns [x,y] in Selig format (unit chord).
        %
        %   Camber line: reflexed (2x3) and non-reflexed (21x, 22x, 23x, 24x, 25x)
        %   Thickness:   standard NACA symmetric (same as 4-series)

            if nargin < 2, N = 60; end
            if isnumeric(des), des = sprintf('%05d', des); end

            % The k1 table (Ladson et al. NASA TM-4741) is calibrated for
            % first digit = 2.  For other first digits scale proportionally.
            % Do NOT multiply by cl directly - that would double-apply scaling.
            first_digit = str2double(des(1));
            camber_scale = first_digit / 2;      % 1.0 for 2xxxx, 0.5 for 1xxxx
            p_pos  = str2double(des(2)) / 20;    % max camber position (0.05..0.25)
            reflex = str2double(des(3));          % 0 = non-reflexed, 1 = reflexed
            t      = str2double(des(4:5)) / 100; % thickness fraction

            % Lookup table for camber line coefficients
            % Non-reflexed (Ladson et al.)
            p_table    = [0.05 0.10 0.15 0.20 0.25];
            m_table    = [0.0580 0.1260 0.2025 0.2900 0.3910];
            k1_table   = [361.40 51.640 15.957 6.643 3.230];

            [~, idx] = min(abs(p_table - p_pos));
            m  = m_table(idx);
            k1 = k1_table(idx);

            % Cosine spacing
            beta = linspace(pi, 0, N+1)';
            xc   = 0.5*(1 - cos(beta));

            % Thickness (NACA symmetric)
            yt = (t/0.2)*(0.2969*sqrt(xc) - 0.1260*xc - 0.3516*xc.^2 ...
                         + 0.2843*xc.^3 - 0.1015*xc.^4);

            % Camber line
            yc  = zeros(size(xc));
            dyc = zeros(size(xc));

            if reflex == 0
                % Standard 5-series (non-reflexed)
                fwd = xc < m;
                bwd = ~fwd;
                yc(fwd)  = (k1/6) * (xc(fwd).^3 - 3*m*xc(fwd).^2 + m^2*(3-m)*xc(fwd));
                yc(bwd)  = (k1*m^3/6) * (1 - xc(bwd));
                dyc(fwd) = (k1/6) * (3*xc(fwd).^2 - 6*m*xc(fwd) + m^2*(3-m));
                dyc(bwd) = -(k1*m^3/6) * ones(sum(bwd),1);
            else
                % Reflexed camber line (230xx series)
                k2_k1 = 0.1;  % approximate
                fwd = xc < m;
                bwd = ~fwd;
                yc(fwd)  = (k1/6) * ((xc(fwd)-m).^3 - k2_k1*(1-m)^3*xc(fwd) - m^3*(xc(fwd)) + m^3);
                yc(bwd)  = (k1/6) * (k2_k1*(xc(bwd)-m).^3 - k2_k1*(1-m)^3*xc(bwd) + k2_k1*(1-m)^3*m);
                dyc(fwd) = (k1/6) * (3*(xc(fwd)-m).^2 - k2_k1*(1-m)^3);
                dyc(bwd) = (k1/6) * (3*k2_k1*(xc(bwd)-m).^2 - k2_k1*(1-m)^3);
            end

            % Apply camber scale (relative to digit=2 which the k1 table encodes)
            yc  = yc  * camber_scale;
            dyc = dyc * camber_scale;

            theta = atan(dyc);
            xu = xc - yt.*sin(theta);
            yu = yc + yt.*cos(theta);
            xl = xc + yt.*sin(theta);
            yl = yc - yt.*cos(theta);

            x = [xu(end:-1:1); xl(2:end)];
            y = [yu(end:-1:1); yl(2:end)];
        end

        % ============================================================== %
        %  CLOSED FUSELAGE SECTION  (via CompositeCurve)
        % ============================================================== %

        function CC = closedFuselageSection(x_sta, width, height, rho_top, rho_bot, rho_side)
        % CLOSEDFUSELAGESECTION  Closed conic cross-section at axial station x.
        %
        %   CC = closedFuselageSection(x_sta, width, height, rho_top, rho_bot, rho_side)
        %
        %   x_sta   - axial (X) position of the section
        %   width   - half-width in Y direction
        %   height  - half-height: [z_top, z_bot] or scalar (symmetric)
        %   rho_top - shoulder parameter for top arc (default 0.5)
        %   rho_bot - shoulder parameter for bottom arc (default 0.5)
        %   rho_side- shoulder parameter for side arcs (default 0.5)
        %
        %   Returns geom.CompositeCurve (4 conic arcs, closed).
        %
        %   Cardinal points:
        %     P_top = [x, 0,      z_top]
        %     P_bot = [x, 0,      z_bot]
        %     P_fwd = [x, +width, 0    ]   (starboard)
        %     P_aft = [x, -width, 0    ]   (port)

            if nargin < 4, rho_top  = 0.5; end
            if nargin < 5, rho_bot  = 0.5; end
            if nargin < 6, rho_side = 0.5; end

            if isscalar(height)
                z_top = height; z_bot = -height;
            else
                z_top = height(1); z_bot = height(2);
            end

            P_top = [x_sta,   0,      z_top];
            P_bot = [x_sta,   0,      z_bot];
            P_fwd = [x_sta,  width,   0    ];
            P_aft = [x_sta, -width,   0    ];

            rhos = [rho_top, rho_side, rho_bot, rho_side];
            CC   = geom.CompositeCurve.fromConicSection(P_top, P_bot, P_fwd, P_aft, rhos);
        end

        % ============================================================== %
        %  WING + FUSELAGE ARRANGEMENT
        % ============================================================== %

        function surfaces = wingBodyLayout(wing_half, fuselage, varargin)
        % WINGBODYLAYOUT  Package wing and fuselage as named WGS network array.
        %
        %   surfaces = wingBodyLayout(wing_half, fuselage)
        %   surfaces = wingBodyLayout(wing_half, fuselage, 'nu', 30, 'nv', 10)
        %
        %   Returns cell array of network structs ready for MeshWriter.toWGS().
        %   Wing is stored with iSymX=1 (XZ symmetry) so PanAir mirrors it.

            pa = inputParser;
            addParameter(pa, 'nu_wing', 25);
            addParameter(pa, 'nv_wing', 10);
            addParameter(pa, 'nu_fus',  20);
            addParameter(pa, 'nv_fus',  16);
            addParameter(pa, 'SpacingChord', 'cosine');
            parse(pa, varargin{:});
            opts = pa.Results;

            surfaces = {};

            if ~isempty(wing_half)
                mesh_w = wing_half.isoMesh(opts.nu_wing, opts.nv_wing, ...
                    'SpacingU', opts.SpacingChord);
                net_w  = geom.MeshWriter.meshToNetwork(mesh_w, 'Wing', 1, 0);
                surfaces{end+1} = net_w;
            end

            if ~isempty(fuselage)
                mesh_f = fuselage.isoMesh(opts.nu_fus, opts.nv_fus);
                net_f  = geom.MeshWriter.meshToNetwork(mesh_f, 'Fuselage', 0, 0);
                surfaces{end+1} = net_f;
            end
        end

    end  % static methods

    % ------------------------------------------------------------------ %
    %  Private helpers
    % ------------------------------------------------------------------ %
    methods (Static, Access = private)

        function [P_net, W_net, V, n_arc] = buildRevolutionNet(profile_pts, O, axis, ...
                                                                 theta0, theta)
        % Build NURBS surface of revolution control net.
        % Uses exact rational quadratic circle arcs.

            % Number of circular arcs needed
            if theta <= pi/2
                narcs = 1;
            elseif theta <= pi
                narcs = 2;
            elseif theta <= 3*pi/2
                narcs = 3;
            else
                narcs = 4;
            end

            n_arc = 2*narcs + 1;  % total control points in v
            dtheta = theta / narcs;
            w1     = cos(dtheta/2);

            % V knot vector for circle
            switch narcs
                case 1
                    V = [0 0 0 1 1 1];
                case 2
                    V = [0 0 0 0.5 0.5 1 1 1];
                case 3
                    V = [0 0 0 1/3 1/3 2/3 2/3 1 1 1];
                case 4
                    V = [0 0 0 0.25 0.25 0.5 0.5 0.75 0.75 1 1 1];
            end

            npts = size(profile_pts, 1);
            P_net= zeros(npts, n_arc, 3);
            W_net= zeros(npts, n_arc);

            for ip = 1:npts
                Pk = profile_pts(ip,:);
                % Project to axis, get radius vector
                OP  = Pk - O;
                axproj = dot(OP, axis) * axis;
                X0  = OP - axproj;   % radial component
                r   = norm(X0);

                if r < 1e-10
                    % On axis: degenerate
                    for j = 1:n_arc
                        P_net(ip,j,:) = Pk;
                        W_net(ip,j)   = 1;
                    end
                    continue;
                end

                X0 = X0 / r;  % radial unit vector
                Y0 = cross(axis, X0);  % tangential unit vector

                % Sweep arcs
                angle = theta0;
                P0 = O + axproj + r*cos(angle)*X0 + r*sin(angle)*Y0;
                j  = 1;
                P_net(ip,j,:) = P0;
                W_net(ip,j)   = 1;

                for a = 1:narcs
                    angle = angle + dtheta;
                    P2  = O + axproj + r*cos(angle)*X0 + r*sin(angle)*Y0;
                    T0  = -sin(angle-dtheta)*X0 + cos(angle-dtheta)*Y0;
                    T2  = -sin(angle)*X0        + cos(angle)*Y0;
                    % Shoulder point P1 = P0 + s*T0 intersect P2 + t*T2
                    A3  = [T0(:), -T2(:)];
                    dp  = (P2-P0)';
                    st  = A3 \ dp;
                    P1  = P0 + st(1)*T0;

                    P_net(ip,j+1,:) = P1;
                    P_net(ip,j+2,:) = P2;
                    W_net(ip,j+1)   = w1;
                    W_net(ip,j+2)   = 1;
                    j  = j + 2;
                    P0 = P2;
                end
            end
        end

    end  % private static

end  % classdef Aircraft