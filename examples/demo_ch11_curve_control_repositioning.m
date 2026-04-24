function demo_ch11_curve_control_repositioning()
%DEMO_CH11_CURVE_CONTROL_REPOSITIONING
clear; clc;
fprintf('=== Chapter 11.2 exact curve control-point repositioning demo ===\n');

P = [ ...
    0.00 0.00 0.0
    0.10 0.00 0.0
    0.18 0.55 0.0
    0.36 0.78 0.0
    0.52 0.92 0.0
    0.68 0.82 0.0
    0.84 0.58 0.0
    0.92 0.00 0.0
    1.00 0.00 0.0 ];
p = 3;
U = [0 0 0 0 0.20 0.35 0.52 0.69 0.84 1 1 1 1];
W = ones(size(P,1),1);

C = geom.NURBSCurve(P, p, U, W);
nodes = geom.curveControlNodes(C);

ubar = 0.43;
deltaP = [0.0, -0.18, 0.0];

P0 = C.evaluate(ubar);

[C1, info1] = geom.curveRepositionOneControl(C, ubar, deltaP);
P1 = C1.evaluate(ubar);

[C2, info2] = geom.curveRepositionInsertedControl(C, ubar, deltaP);
P2 = C2.evaluate(ubar);

[C3, info3] = geom.curveRepositionTwoControls(C, ubar, deltaP);
P3 = C3.evaluate(ubar);

fprintf('Picked parameter ubar                  = %.6f\n', ubar);
fprintf('Original point C(ubar)                = [% .6f % .6f]\n', P0(1), P0(2));
fprintf('Desired translated point              = [% .6f % .6f]\n', P0(1)+deltaP(1), P0(2)+deltaP(2));
fprintf('Single-control exact error            = %.3e\n', norm(P1 - (P0 + deltaP)));
fprintf('Insert-and-move exact error           = %.3e\n', norm(P2 - (P0 + deltaP)));
fprintf('Two-control exact error               = %.3e\n', norm(P3 - (P0 + deltaP)));
fprintf('Single-control chosen k               = %d\n', info1.k);
fprintf('Inserted-knot value uhat              = %.6f\n', info2.insert.uhat);
fprintf('Two-control chosen k, gamma           = %d, %.6f\n', info3.k, info3.gamma);

figure('Color','w','Name','Chapter 11.2 exact curve repositioning');
tiledlayout(2,2, 'Padding','compact', 'TileSpacing','compact');

nexttile;
hold on;
localPlotCurveAndPolygon(C, [0 0 0], [0.3 0.3 0.3]);
localMarkNodes(C, nodes);
plot(P0(1), P0(2), 'ko', 'MarkerFaceColor','y', 'MarkerSize',7);
plot(P0(1)+deltaP(1), P0(2)+deltaP(2), 'rx', 'LineWidth',2, 'MarkerSize',10);
title({'Original curve', 'yellow = picked point, red x = target point'});
axis equal; grid on;

nexttile;
hold on;
localPlotCurveAndPolygon(C, [0.65 0.65 0.65], [0.8 0.8 0.8]);
localPlotCurveAndPolygon(C1, [0 0.3 0.9], [0 0.3 0.9]);
plot(P0(1), P0(2), 'ko', 'MarkerFaceColor','y', 'MarkerSize',7);
plot(P1(1), P1(2), 'bo', 'MarkerFaceColor','c', 'MarkerSize',7);
plot(P0(1)+deltaP(1), P0(2)+deltaP(2), 'rx', 'LineWidth',2, 'MarkerSize',10);
plot(C.P(info1.k,1),  C.P(info1.k,2),  'ks', 'MarkerFaceColor','k', 'MarkerSize',7);
plot(C1.P(info1.k,1), C1.P(info1.k,2), 'bs', 'MarkerFaceColor','b', 'MarkerSize',7);
title({'One control point moved', sprintf('Eq. 11.3, k = %d', info1.k)});
axis equal; grid on;

nexttile;
hold on;
localPlotCurveAndPolygon(C, [0.65 0.65 0.65], [0.8 0.8 0.8]);
localPlotCurveAndPolygon(C2, [0.7 0 0.7], [0.7 0 0.7]);
plot(P0(1), P0(2), 'ko', 'MarkerFaceColor','y', 'MarkerSize',7);
plot(P2(1), P2(2), 'mo', 'MarkerFaceColor','m', 'MarkerSize',7);
plot(P0(1)+deltaP(1), P0(2)+deltaP(2), 'rx', 'LineWidth',2, 'MarkerSize',10);
q = info2.insert.qIndex;
plot(C2.P(q,1), C2.P(q,2), 'ms', 'MarkerFaceColor','m', 'MarkerSize',7);
title({'Insert knot, then move new control point', sprintf('Eq. 11.6, q = %d', q)});
axis equal; grid on;

nexttile;
hold on;
localPlotCurveAndPolygon(C, [0.65 0.65 0.65], [0.8 0.8 0.8]);
localPlotCurveAndPolygon(C3, [0.85 0 0], [0.85 0 0]);
plot(P0(1), P0(2), 'ko', 'MarkerFaceColor','y', 'MarkerSize',7);
plot(P3(1), P3(2), 'ro', 'MarkerFaceColor','r', 'MarkerSize',7);
plot(P0(1)+deltaP(1), P0(2)+deltaP(2), 'rx', 'LineWidth',2, 'MarkerSize',10);
plot(C3.P(info3.k,1),   C3.P(info3.k,2),   'rs', 'MarkerFaceColor','r', 'MarkerSize',7);
plot(C3.P(info3.k+1,1), C3.P(info3.k+1,2), 'rs', 'MarkerFaceColor','r', 'MarkerSize',7);
title({'Two neighboring control points moved', sprintf('Eq. 11.7, gamma = %.3f', info3.gamma)});
axis equal; grid on;

fprintf('Demo complete.\n');
end

function localPlotCurveAndPolygon(C, curveColor, polyColor)
u = linspace(C.U(C.p+1), C.U(end-C.p), 400).';
X = zeros(numel(u), 3);
for i = 1:numel(u)
    X(i,:) = C.evaluate(u(i));
end
plot(X(:,1), X(:,2), '-', 'Color', curveColor, 'LineWidth', 2);
plot(C.P(:,1), C.P(:,2), ':', 'Color', polyColor, 'LineWidth', 1.2);
plot(C.P(:,1), C.P(:,2), 'o', 'Color', polyColor, 'MarkerSize', 4, ...
    'MarkerFaceColor', [0.95 0.95 0.95]);
end

function localMarkNodes(C, nodes)
for i = 1:numel(nodes)
    Pi = C.evaluate(nodes(i));
    text(Pi(1), Pi(2), sprintf('  t_%d', i-1), 'FontSize', 8);
end
end
