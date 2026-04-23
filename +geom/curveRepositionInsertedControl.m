function [C2, info] = curveRepositionInsertedControl(C, ubar, deltaP, k)
% geom.curveRepositionInsertedControl
% Knot-insertion refinement path from The NURBS Book Sec. 11.2
%
% If k is omitted, it is chosen inside curveInsertKnotForNode.

if nargin < 4 || isempty(k)
    [Cref, insInfo] = geom.curveInsertKnotForNode(C, ubar);
else
    [Cref, insInfo] = geom.curveInsertKnotForNode(C, ubar, k);
end

q = insInfo.qIndex;
[C2, repInfo] = geom.curveRepositionOneControl(Cref, ubar, deltaP, q);

info = struct();
info.insert = insInfo;
info.reposition = repInfo;
info.equation = 'Implements the knot-insertion refinement path of Sec. 11.2.';
end
