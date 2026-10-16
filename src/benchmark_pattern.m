function benchmark_pattern
a=43.6e-3; f=340e9; el0=asind(1/sqrt(3));
angles=[45,el0;22,12;80,10;5,75;45,70;1,1;89,89;45,atand(1/sqrt(2));45,atand(sqrt(2))];
for i=1:size(angles,1)
    z=sbr_pattern_point(a,f,angles(i,1),angles(i,2),192);
    hh=sbr_trihedral(a,f,angles(i,1),angles(i,2),192,'HH');
    vv=sbr_trihedral(a,f,angles(i,1),angles(i,2),192,'VV');
    dh=10*log10(z.HH)-hh.rcs_dBsm;
    dv=10*log10(z.VV)-vv.rcs_dBsm;
    fprintf('az %.6g el %.6g: HH %.9g VV %.9g differences %.3g/%.3g dB\n',...
        angles(i,:),10*log10(z.HH),10*log10(z.VV),dh,dv);
    assert(max(abs([dh,dv]))<1e-7);
end
rng(340436);
r0=rand(10000,2); mask=sum(r0,2)>1; r0(mask,:)=1-r0(mask,:);
r0=a*[r0,1-sum(r0,2)];
for i=1:40
    az=0.1+89.8*rand; el=0.1+89.8*rand;
    d=-[cosd(el)*cosd(az),cosd(el)*sind(az),sind(el)];
    [r1,d1,L1,b1,f1]=sbr_bounce(r0,d,a);
    [r2,d2,L2,b2,f2]=sbr_bounce_fast(r0,d,a);
    assert(max(abs(r1-r2),[],'all')<a*1e-9);
    assert(max(abs(L1-L2))<a*1e-9);
    assert(isequal(b1,b2) && isequal(d1,d2) && isequal(f1,f2));
end
fprintf('400000 random ray trajectories matched the original intersection tracer.\n');
tic;
for i=1:30
    sbr_pattern_point(a,f,10+2*i,10+i,384);
end
fprintf('N384 average seconds per dual-pol point = %.6g\n',toc/30);
pool=parpool('threads',4);
tic;
parfor i=1:60
    sbr_pattern_point(a,f,10+i,10+i,384);
end
fprintf('N384 thread pool seconds per dual-pol point = %.6g\n',toc/60);
delete(pool);
end
