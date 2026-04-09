classdef Demo
% DEMO  Helper functions for the geometry kernel demo script.

    methods (Static)

        function [x, y] = naca4Coords(m, p, t, N)
        % NACA4COORDS  Generate NACA 4-digit airfoil coordinates.
        %
        %   [x, y] = naca4Coords(m, p, t, N)
        %
        %   m  - max camber as fraction of chord (e.g. 0.02 for NACA 2xxx)
        %   p  - location of max camber as fraction (e.g. 0.4 for NACAxxx)
        %   t  - max thickness as fraction of chord (e.g. 0.12 for NACA xx12)
        %   N  - number of points on each surface (default 50)
        %
        %   Returns [x,y] in Selig format (TE -> upper -> LE -> lower -> TE),
        %   unit chord, x in [0,1].

            if nargin < 4, N = 50; end

            % Cosine spacing (denser near LE/TE)
            beta = linspace(0, pi, N+1)';
            xc   = 0.5*(1 - cos(beta));   % chord stations

            % Thickness distribution (NACA symmetric formula)
            yt = (t/0.2) * (0.2969*sqrt(xc) ...
                          - 0.1260*xc ...
                          - 0.3516*xc.^2 ...
                          + 0.2843*xc.^3 ...
                          - 0.1015*xc.^4);

            % Camber line and slope
            yc    = zeros(size(xc));
            dyc   = zeros(size(xc));

            if abs(m) > 1e-6 && abs(p) > 1e-6
                fwd = xc <= p;
                bwd = ~fwd;
                yc(fwd)  = (m/p^2)     * (2*p*xc(fwd)  - xc(fwd).^2);
                yc(bwd)  = (m/(1-p)^2) * ((1-2*p) + 2*p*xc(bwd) - xc(bwd).^2);
                dyc(fwd) = (2*m/p^2)   * (p - xc(fwd));
                dyc(bwd) = (2*m/(1-p)^2) * (p - xc(bwd));
            end

            theta = atan(dyc);

            % Upper and lower surface
            xu = xc - yt.*sin(theta);
            yu = yc + yt.*cos(theta);
            xl = xc + yt.*sin(theta);
            yl = yc - yt.*cos(theta);

            % Assemble Selig format: upper surface TE->LE, then lower LE->TE
            x = [xu(end:-1:1); xl(2:end)];
            y = [yu(end:-1:1); yl(2:end)];
        end

        function plotAirfoil(m, p, t, N)
        % PLOTAIRFOIL  Quick plot of a NACA 4-digit airfoil.

            if nargin < 4, N = 100; end
            [x,y] = geom.Demo.naca4Coords(m, p, t, N);
            figure;
            plot(x, y, 'b-', 'LineWidth', 1.5); hold on;
            fill(x, y, [0.8 0.85 0.95], 'FaceAlpha', 0.5, 'EdgeColor','none');
            axis equal; grid on;
            title(sprintf('NACA %d%d%02d', round(m*100), round(p*10), round(t*100)));
            xlabel('x/c'); ylabel('y/c');
        end

    end

end  % classdef Demo
