function result = sbr_trihedral(a,f,az,el,N,pol,options)
% Custom PEC SBR + vector Huygens exit-aperture integration.
% NOT MathWorks raytrace(), NOT a full-wave solver, NO fitted RCS constants.
% For broad angular sweeps use sbr_pattern_exit/run_pattern_sweep. Entrance
% sampling here can become ill-conditioned for grazing outgoing branches;
% the uniform exit-aperture integrator removes that sampling singularity.
% a: inner leg [m], f: [Hz], az/el: radar direction [deg], pol: HH or VV.
% Centroid quadrature uses N^2 equal-area triangular ray tubes. Monte Carlo
% uses N independent uniform aperture samples, with a reproducible seed.
arguments
    a (1,1) double {mustBePositive}
    f (1,1) double {mustBePositive}
    az (1,1) double
    el (1,1) double
    N (1,1) double {mustBeInteger,mustBePositive}
    pol (1,1) string {mustBeMember(pol,["HH","VV"])} = "HH"
    options.Sampling (1,1) string {mustBeMember(options.Sampling,["centroid","montecarlo"])} = "centroid"
    options.Seed (1,1) double = 340436
    options.BlockSize (1,1) double = 32768
end
started = tic;
lambda = 299792458/f;
k = 2*pi/lambda;
s = [cosd(el)*cosd(az),cosd(el)*sind(az),sind(el)];
assert(all(s>0),'This aperture model is limited to incidence in the open first octant.');
din = -s;
normal = [1,1,1]/sqrt(3);
h = cross(din,[0,0,1]); h = h/norm(h);
v = cross(h,din); v = v/norm(v);
if pol=="HH", e0=h; ex=v; else, e0=v; ex=h; end
cosIn = dot(s,normal);
apertureArea = sqrt(3)*a^2/2;
if options.Sampling=="centroid"
    nRays = N^2;
    uv = zeros(nRays,2);
    index = 0;
    for i=0:N-1
        jj = (0:N-1-i)'; n1 = numel(jj);
        uv(index+(1:n1),:) = [(i+1/3)*ones(n1,1),jj+1/3]/N;
        index = index+n1;
        jj = (0:N-2-i)'; n2 = numel(jj);
        uv(index+(1:n2),:) = [(i+2/3)*ones(n2,1),jj+2/3]/N;
        index = index+n2;
    end
    assert(index==nRays);
else
    nRays = N;
    stream = RandStream('mt19937ar','Seed',options.Seed);
    uv = rand(stream,nRays,2);
    folded = sum(uv,2)>1;
    uv(folded,:) = 1-uv(folded,:);
end
dAin = apertureArea/nRays;
amplitude = complex(zeros(3,3));
counts = zeros(3,1);
maxTransverse = 0;
maxTriplePhaseResidual = 0;
maxTripleDirectionError = 0;
for first=1:options.BlockSize:nRays
    last = min(first+options.BlockSize-1,nRays);
    q = uv(first:last,:);
    r0 = a*[q,1-sum(q,2)];
    [r,d,L,b,flip] = sbr_bounce(r0,din,a);
    % PEC: E_reflected = -E + 2*n*(n.E); product of planar reflections.
    E = ((-1).^b).*flip.*e0;
    maxTransverse = max(maxTransverse,max(abs(sum(E.*d,2))));
    cosOut = d*normal';
    assert(all(cosOut>0),'Exited ray is not outgoing.');
    % Ray-tube flux conservation: dAout*cosOut = dAin*cosIn.
    dAout = dAin*cosIn./cosOut;
    opticalPhaseLength = r0*din' + L - r*s';
    F = sbr_aperture_kernel(E,d,normal,s);
    integrand = F.*(dAout.*exp(-1i*k*opticalPhaseLength));
    for order=1:3
        take = b==order;
        counts(order) = counts(order)+nnz(take);
        amplitude(order,:) = amplitude(order,:)+sum(integrand(take,:),1);
    end
    take = b==3;
    if any(take)
        maxTriplePhaseResidual = max(maxTriplePhaseResidual,max(abs(opticalPhaseLength(take))));
        maxTripleDirectionError = max(maxTripleDirectionError,max(vecnorm(d(take,:)-s,2,2)));
    end
end
assert(sum(counts)==nRays);
copol = amplitude*e0';
crosspol = amplitude*ex';
sigma = 4*pi/lambda^2*abs(sum(copol))^2;
sigmaCross = 4*pi/lambda^2*abs(sum(crosspol))^2;
result = struct('method','PEC SBR plus vector exit-aperture integration',...
    'is_full_wave',false,'inner_leg_m',a,'frequency_Hz',f,...
    'azimuth_deg',az,'elevation_deg',el,'polarization',char(pol),...
    'sampling',char(options.Sampling),'N',N,'rays',nRays,...
    'rcs_m2',sigma,'rcs_dBsm',10*log10(sigma),...
    'cross_rcs_m2',sigmaCross,'cross_rcs_dBsm',10*log10(max(sigmaCross,realmin)),...
    'rays_by_order',counts','copol_amplitude_real_m2',real(copol)',...
    'copol_amplitude_imag_m2',imag(copol)',...
    'rcs_by_order_m2',(4*pi/lambda^2*abs(copol).^2)',...
    'triple_projected_area_m2',counts(3)*dAin*cosIn,...
    'max_E_dot_direction',maxTransverse,...
    'max_triple_phase_error_m',maxTriplePhaseResidual,...
    'max_triple_direction_error',maxTripleDirectionError,...
    'elapsed_seconds',toc(started));
end
