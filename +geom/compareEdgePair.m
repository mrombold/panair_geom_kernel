function cmp = compareEdgePair(PA, PB)
cmp = struct( ...
    'sameError', inf, ...
    'revError',  inf, ...
    'bestError', inf, ...
    'mode',      'none', ...
    'sameCount', false);

if size(PA,2) ~= 3 || size(PB,2) ~= 3
    return;
end

if size(PA,1) ~= size(PB,1)
    cmp.sameCount = false;
    return;
end

cmp.sameCount = true;

cmp.sameError = max(vecnorm(PA - PB, 2, 2));
cmp.revError  = max(vecnorm(PA - flipud(PB), 2, 2));

if cmp.sameError <= cmp.revError
    cmp.bestError = cmp.sameError;
    cmp.mode = 'same';
else
    cmp.bestError = cmp.revError;
    cmp.mode = 'reversed';
end
end