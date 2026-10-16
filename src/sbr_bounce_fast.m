function [r,d,L,b,flip] = sbr_bounce_fast(r0,din,a)
% Exact event times for the SAME three orthogonal finite-plane ray problem.
% The unfolded trajectory is r_j(t)=abs(r0_j-s_j*t), s=-din>0.
% Candidate metal events are t_j=r0_j/s_j. The open aperture event is
% sum_j abs(r0_j-s_j*t)=a. Test its linear root in each of the 1-,2-,3-
% reflection intervals. This eliminates iterative intersection bookkeeping.
s=-din; S=sum(s); nr=size(r0,1); tol=a*1e-11;
times=r0./s;
[t1,first]=min(times,[],2);
[t3,reverseLast]=max(times(:,[3,2,1]),[],2);
last=4-reverseLast;
t2=sum(times,2)-t1-t3;
row=(1:nr)';
entryFirst=r0(sub2ind([nr,3],row,first));
entryLast=r0(sub2ind([nr,3],row,last));
den1=2*reshape(s(first),[],1)-S;
den2=S-2*reshape(s(last),[],1);
L1=inf(nr,1); ok=den1>1e-13; L1(ok)=2*entryFirst(ok)./den1(ok);
L2=inf(nr,1); ok=den2>1e-13; L2(ok)=2*(a-entryLast(ok))./den2(ok);
one=L1<=t2+tol;
two=~one & L2<=t3+tol;
b=3*ones(nr,1); b(one)=1; b(two)=2;
L=(2*a/S)*ones(nr,1); L(one)=L1(one); L(two)=L2(two);
flip=-ones(nr,3);
flip(one,:)=1;
flip(sub2ind([nr,3],row(one),first(one)))=-1;
flip(sub2ind([nr,3],row(two),last(two)))=1;
r=(r0-L.*s).*flip;
d=din.*flip;
assert(min(r(:))>=-tol && max(abs(sum(r,2)-a))<a*1e-9);
end
