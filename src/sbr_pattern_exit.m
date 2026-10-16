function out=sbr_pattern_exit(a,f,az,el,N)
% SAME PEC GO/Huygens model, integrated on uniformly sampled EXIT aperture.
% Avoids the 1/(n.d_out) quadrature weight singularity of entrance sampling.
% Every possible reflected plane wave is tested by reverse ray tracing.
% Orthogonal planes permit an exact event test: reverse rays must reach all
% selected walls before re-exiting the aperture, t_last <= t_exit.
persistent cachedN unitPoints
if isempty(cachedN) || cachedN~=N
    cachedN=N; uv=zeros(N^2,2); index=0;
    for i=0:N-1
        j=(0:N-1-i)'; n=numel(j);
        uv(index+(1:n),:)=[(i+1/3)*ones(n,1),j+1/3]/N; index=index+n;
        j=(0:N-2-i)'; n=numel(j);
        uv(index+(1:n),:)=[(i+2/3)*ones(n,1),j+2/3]/N; index=index+n;
    end
    unitPoints=[uv,1-sum(uv,2)];
end
s=[cosd(el)*cosd(az),cosd(el)*sind(az),sind(el)];
assert(all(s>0)); din=-s; normal=[1,1,1]/sqrt(3);
h=cross(din,[0,0,1]); h=h/norm(h); v=cross(h,din); v=v/norm(v);
lambda=299792458/f; k=2*pi/lambda; S=sum(s);
dA=sqrt(3)*a^2/(2*N^2);
amplitude=complex(zeros(2)); triple=complex(zeros(2)); counts=zeros(7,1);
for state=1:7
    bits=logical(bitget(state,1:3)); flips=1-2*double(bits); order=nnz(bits);
    direction=din.*flips;
    if dot(direction,normal)<=1e-13, continue; end
    % Reverse propagation starts along -direction. Only selected components
    % are negative; after those reflections its direction is +s.
    L=2*sum(unitPoints(:,bits),2)/S;
    lastHit=max(unitPoints(:,bits)./s(bits),[],2);
    valid=lastHit<L-1e-11; % exclude measure-zero ray/edge ties
    counts(state)=nnz(valid);
    if ~any(valid), continue; end
    % Planes pass through the origin, so d_in.r_in + path = d_out.r_out.
    % This exact plane-wave phase avoids long-path cancellation numerics.
    phase=(k*a)*(unitPoints(valid,:)*(direction-s)');
    integral=dA*sum(exp(-1i*phase));
    reflectedE=((-1)^order)*[h;v].*flips;
    kernel=sbr_aperture_kernel(reflectedE,repmat(direction,2,1),normal,s);
    component=(kernel*[h;v]')*integral;
    amplitude=amplitude+component;
    if state==7, triple=component; end
end
factor=4*pi/lambda^2;
out=struct('HH',factor*abs(amplitude(1,1))^2,'VV',factor*abs(amplitude(2,2))^2,...
    'HV',factor*abs(amplitude(1,2))^2,'VH',factor*abs(amplitude(2,1))^2,...
    'tripleHH',factor*abs(triple(1,1))^2,'amplitude',amplitude,...
    'rays',N^2,'counts',counts);
end
