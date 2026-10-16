function [r,d,pathLength,bounces,flip,paths] = sbr_bounce(r0,din,a,recordPaths)
% Trace rays inside x,y,z>=0, x+y+z<=a. The x+y+z=a face is OPEN.
% Each metallic face is a finite right-isosceles triangle, not an infinite
% reflector. Exit is an explicit competing intersection at every bounce.
% Planar PEC reflections preserve ray-tube cross-sectional area.
if nargin<4, recordPaths=false; end
n = size(r0,1);
r = r0;
d = repmat(din,n,1);
flip = ones(n,3);
pathLength = zeros(n,1);
bounces = zeros(n,1);
active = true(n,1);
tol = a*1e-11;
paths = [];
if recordPaths
    paths = nan(n,3,5);
    paths(:,:,1) = r0;
end
for step=1:4
    ii = find(active);
    if isempty(ii), break; end
    rr = r(ii,:); dd = d(ii,:);
    wallTime = inf(numel(ii),3);
    towardWall = dd < -1e-13;
    wallTime(towardWall) = -rr(towardWall)./dd(towardWall);
    wallTime(wallTime<0 & wallTime>-tol) = 0;
    assert(all(wallTime(:)>=0),'Ray left the physical cavity.');
    [tw,wall] = min(wallTime,[],2);
    outSpeed = sum(dd,2);
    te = inf(numel(ii),1);
    outward = outSpeed>1e-13;
    te(outward) = (a-sum(rr(outward,:),2))./outSpeed(outward);
    te(te<0 & te>-tol) = 0;
    assert(all(te>=0),'Negative aperture exit distance.');
    % Resolve a boundary tie by exit. Such edge rays have measure zero;
    % the independent Monte Carlo run avoids deterministic grid alignment.
    exitNow = te<=tw+tol;
    travel = tw;
    travel(exitNow) = te(exitNow);
    assert(all(isfinite(travel)),'Ray has no next intersection.');
    r(ii,:) = rr + travel.*dd;
    pathLength(ii) = pathLength(ii)+travel;
    if recordPaths, paths(ii,:,step+1) = r(ii,:); end
    active(ii(exitNow)) = false;
    jj = ii(~exitNow);
    ww = wall(~exitNow);
    if ~isempty(jj)
        linearIndex = sub2ind([n,3],jj,ww);
        r(linearIndex) = 0;
        d(linearIndex) = -d(linearIndex);
        flip(linearIndex) = -flip(linearIndex);
        bounces(jj) = bounces(jj)+1;
    end
end
assert(~any(active),'A ray did not exit after at most three reflections.');
assert(all(bounces>=1 & bounces<=3),'Unexpected reflection count.');
assert(max(abs(sum(r,2)-a))<a*1e-9,'Exit is not on the aperture.');
assert(min(r(:))>=-tol,'Exit is outside the triangular aperture.');
end
