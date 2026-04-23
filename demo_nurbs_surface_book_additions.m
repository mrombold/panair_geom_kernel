function demo_nurbs_surface_book_additions()
%DEMO_NURBS_SURFACE_BOOK_ADDITIONS
% Exercise added NURBSSurface utilities:
%   * reverseU / reverseV
%   * reduceU / reduceV aliases

clear; clc;
fprintf('=== NURBSSurface book-additions demo ===\n');

%% 1) Build a simple bicubic tensor-product surface
[xg, yg] = ndgrid(linspace(0, 1.5, 4), linspace(0, 1.0, 4));
zg = 0.15 * sin(pi*xg/1.5) .* cos(pi*yg/1.0) + 0.05*xg.^2 - 0.08*yg;

P = zeros(4,4,3);
P(:,:,1) = xg;
P(:,:,2) = yg;
P(:,:,3) = zg;

W = ones(4,4);
p = 3; q = 3;
U = [0 0 0 0 1 1 1 1];
V = [0 0 0 0 1 1 1 1];

S = geom.NURBSSurface(P, p, q, U, V, W);

%% 2) Add a rectangular trim loop in UV if setTrims exists
if any(strcmp(methods(S), 'setTrims'))
    outer = [ ...
        0.15 0.15
        0.85 0.15
        0.85 0.85
        0.15 0.85
        0.15 0.15];
    S = S.setTrims({outer}, {});
end

%% 3) reverseU invariance check
Su = S.reverseU();

uv = [ ...
    0.17 0.22
    0.30 0.70
    0.61 0.41
    0.92 0.18];

u = uv(:,1);
v = uv(:,2);

a = S.domainU(1);
b = S.domainU(2);

X0 = S.evaluate(a + b - u, v);
X1 = Su.evaluate(u, v);
errU = max(vecnorm(X0 - X1, 2, 2));

fprintf('reverseU max geometry mismatch = %.3e\n', errU);

if any(strcmp(methods(S), 'isInsideTrim'))
    trimU0 = arrayfun(@(k) S.isInsideTrim(a+b-u(k), v(k)), 1:numel(u)).';
    trimU1 = arrayfun(@(k) Su.isInsideTrim(u(k), v(k)), 1:numel(u)).';
    fprintf('reverseU trim-mask agreement   = %d / %d\n', nnz(trimU0 == trimU1), numel(u));
end

%% 4) reverseV invariance check
Sv = S.reverseV();

c = S.domainV(1);
d = S.domainV(2);

X0 = S.evaluate(u, c + d - v);
X1 = Sv.evaluate(u, v);
errV = max(vecnorm(X0 - X1, 2, 2));

fprintf('reverseV max geometry mismatch = %.3e\n', errV);

if any(strcmp(methods(S), 'isInsideTrim'))
    trimV0 = arrayfun(@(k) S.isInsideTrim(u(k), c+d-v(k)), 1:numel(u)).';
    trimV1 = arrayfun(@(k) Sv.isInsideTrim(u(k), v(k)), 1:numel(u)).';
    fprintf('reverseV trim-mask agreement   = %d / %d\n', nnz(trimV0 == trimV1), numel(u));
end

%% 5) Double reversal should recover the original surface
Suu = Su.reverseU();
Svv = Sv.reverseV();

uv2 = [ ...
    0.05 0.10
    0.25 0.55
    0.50 0.50
    0.73 0.81
    0.97 0.25];

Xref = S.evaluate(uv2(:,1), uv2(:,2));
Xuu  = Suu.evaluate(uv2(:,1), uv2(:,2));
Xvv  = Svv.evaluate(uv2(:,1), uv2(:,2));

errUU = max(vecnorm(Xref - Xuu, 2, 2));
errVV = max(vecnorm(Xref - Xvv, 2, 2));

fprintf('reverseU(reverseU(S)) mismatch = %.3e\n', errUU);
fprintf('reverseV(reverseV(S)) mismatch = %.3e\n', errVV);

%% 6) Exercise the reduceU/reduceV aliases
SeU = S.elevateU(1);
[SrU, errRedU] = SeU.reduceU(1, inf, 200);

SeV = S.elevateV(1);
[SrV, errRedV] = SeV.reduceV(1, inf, 200);

Xeu = SeU.evaluate(uv2(:,1), uv2(:,2));
Xru = SrU.evaluate(uv2(:,1), uv2(:,2));
XuMismatch = max(vecnorm(Xeu - Xru, 2, 2));

Xev = SeV.evaluate(uv2(:,1), uv2(:,2));
Xrv = SrV.evaluate(uv2(:,1), uv2(:,2));
XvMismatch = max(vecnorm(Xev - Xrv, 2, 2));

fprintf('reduceU alias reported err     = %.3e\n', errRedU);
fprintf('reduceU sampled mismatch       = %.3e\n', XuMismatch);
fprintf('reduceV alias reported err     = %.3e\n', errRedV);
fprintf('reduceV sampled mismatch       = %.3e\n', XvMismatch);

%% 7) Quick visualization, if plot exists
if any(strcmp(methods(S), 'plot'))
    figure('Color', 'w', 'Name', 'NURBSSurface book additions');
    tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    try
        S.plot(35, 35, 'ShowCP', true, 'ShowTrims', true);
    catch
        S.plot(35, 35);
    end
    title('Original S');
    axis equal; view(3);

    nexttile;
    try
        Su.plot(35, 35, 'ShowCP', true, 'ShowTrims', true);
    catch
        Su.plot(35, 35);
    end
    title('reverseU(S)');
    axis equal; view(3);

    nexttile;
    try
        Sv.plot(35, 35, 'ShowCP', true, 'ShowTrims', true);
    catch
        Sv.plot(35, 35);
    end
    title('reverseV(S)');
    axis equal; view(3);
end

fprintf('Demo complete.\n');
end
