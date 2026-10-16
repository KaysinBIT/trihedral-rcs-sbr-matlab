function out = sbr_pattern_point(a,f,az,el,N)
% Same ray/exit-aperture model as sbr_trihedral, with both polarizations.
% Planar-ray reflection states share one polarization matrix and outgoing
% direction. Grouping their complex aperture integrals accelerates sweeps;
% no angular RCS formula replaces the numerical ray integration.
persistent cachedN cachedUV
if isempty(cachedN) || cachedN~=N
    cachedN=N;
    cachedUV=zeros(N^2,2); index=0;
    for i=0:N-1
        j=(0:N-1-i)'; n=numel(j);
        cachedUV(index+(1:n),:)=[(i+1/3)*ones(n,1),j+1/3]/N; index=index+n;
        j=(0:N-2-i)'; n=numel(j);
        cachedUV(index+(1:n),:)=[(i+2/3)*ones(n,1),j+2/3]/N; index=index+n;
    end
end
lambda=299792458/f; k=2*pi/lambda;
s=[cosd(el)*cosd(az),cosd(el)*sind(az),sind(el)];
assert(all(s>0),'Only open-first-octant directions are supported.');
din=-s; normal=[1,1,1]/sqrt(3);
h=cross(din,[0,0,1]); h=h/norm(h);
v=cross(h,din); v=v/norm(v);
cosIn=dot(s,normal); dAin=sqrt(3)*a^2/(2*N^2);
rayIntegral=complex(zeros(7,1)); counts=zeros(7,1);
for first=1:65536:N^2
    last=min(first+65535,N^2);
    q=cachedUV(first:last,:); r0=a*[q,1-sum(q,2)];
    [r,d,L,~,flip]=sbr_bounce_fast(r0,din,a);
    state=round(((1-flip)/2)*[1;2;4]);
    assert(all(state>=1 & state<=7));
    phaseLength=r0*din'+L-r*s';
    weight=dAin*cosIn./(d*normal');
    rayIntegral=rayIntegral+accumarray(state,weight.*cos(k*phaseLength),[7,1],@sum,0)...
        -1i*accumarray(state,weight.*sin(k*phaseLength),[7,1],@sum,0);
    counts=counts+accumarray(state,1,[7,1],@sum,0);
end
amplitude=complex(zeros(2,2)); triple=complex(zeros(2,2));
for state=1:7
    if counts(state)==0, continue; end
    bits=double(bitget(state,1:3)); flips=1-2*bits;
    direction=din.*flips;
    reflectedE=((-1)^sum(bits))*[h;v].*flips;
    kernel=sbr_aperture_kernel(reflectedE,repmat(direction,2,1),normal,s);
    % rows = transmitted H,V; columns = received H,V.
    component=(kernel*[h;v]')*rayIntegral(state);
    amplitude=amplitude+component;
    if state==7, triple=component; end
end
factor=4*pi/lambda^2;
out=struct('HH',factor*abs(amplitude(1,1))^2,...
    'VV',factor*abs(amplitude(2,2))^2,...
    'HV',factor*abs(amplitude(1,2))^2,...
    'VH',factor*abs(amplitude(2,1))^2,...
    'tripleHH',factor*abs(triple(1,1))^2,...
    'amplitude',amplitude,'rays',N^2,'counts',counts);
end
