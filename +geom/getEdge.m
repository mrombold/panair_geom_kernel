function P = getEdge(mesh, name)
switch lower(strtrim(name))
    case 'i1'
        P = [mesh.X(1,:).',   mesh.Y(1,:).',   mesh.Z(1,:).'];
    case {'in','iN'}
        P = [mesh.X(end,:).', mesh.Y(end,:).', mesh.Z(end,:).'];
    case 'j1'
        P = [mesh.X(:,1),     mesh.Y(:,1),     mesh.Z(:,1)];
    case {'jn','jN'}
        P = [mesh.X(:,end),   mesh.Y(:,end),   mesh.Z(:,end)];
    otherwise
        error('geom.getEdge: unknown edge name "%s".', name);
end
end