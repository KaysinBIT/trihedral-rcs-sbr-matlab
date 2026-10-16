function validate_exit_sampling
src=fileparts(mfilename('fullpath')); outdir=pattern_output_dir;
a=43.6e-3; f=340e9;
% Independent iterative reverse tracer checks the fast visibility test and
% plane-wave phase identity, including near-grazing output directions.
rng(340436); uv=rand(3000,2); flip=sum(uv,2)>1; uv(flip,:)=1-uv(flip,:);
r0=a*[uv,1-sum(uv,2)]; normal=[1,1,1]/sqrt(3);
for angle=1:25
    az=0.5+89*rand; el=0.5+89*rand;
    s=[cosd(el)*cosd(az),cosd(el)*sind(az),sind(el)]; S=sum(s);
    for state=1:7
        bits=logical(bitget(state,1:3)); flips=1-2*double(bits);
        d=-s.*flips;
        if dot(d,normal)<=1e-13, continue; end
        [rin,reverseOut,L,~,hitFlip]=sbr_bounce(r0,-d,a);
        traced=all(hitFlip==flips,2);
        tested=max(r0(:,bits)./s(bits),[],2)<2*sum(r0(:,bits),2)/S-a*1e-11;
        assert(isequal(traced,tested),'Reverse ray visibility mismatch.');
        if any(tested)
            assert(max(abs(reverseOut(tested,:)-s),[],'all')<1e-12);
            phaseError=-rin(tested,:)*s'+L(tested)-r0(tested,:)*d';
            assert(max(abs(phaseError))<a*1e-10,'Plane-wave phase identity mismatch.');
        end
    end
end
fprintf('Reverse visibility and phase checks passed for 25 random directions x 3000 aperture points.\n');
angles=[45,asind(1/sqrt(3));26.5,24;63.5,24;45,56.5;45,60;22,12;1,45];
levels=[192,384,768,1536]; rows=[];
for i=1:size(angles,1)
    for N=levels
        old=sbr_pattern_point(a,f,angles(i,1),angles(i,2),N);
        fresh=sbr_pattern_exit(a,f,angles(i,1),angles(i,2),N);
        rows(end+1,:)=[angles(i,:),N,10*log10(old.HH),10*log10(fresh.HH),10*log10(fresh.VV)]; %#ok<AGROW>
    end
end
T=array2table(rows,'VariableNames',{'Azimuth_deg','Elevation_deg','N','EntranceHH_dBsm','ExitHH_dBsm','ExitVV_dBsm'});
disp(T); writetable(T,fullfile(outdir,'entrance_exit_sampling_comparison.csv'));
end
