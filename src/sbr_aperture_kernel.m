function F = sbr_aperture_kernel(E, direction, normal, observation)
% Huygens exit-aperture kernel, normalized so F=E for normal outgoing waves.
% Time convention exp(+jwt); eta*H = direction x E.
% Js = n x H, Ms = -n x E. The far field is
% E_s = jk exp(-jkr)/(2*pi*r) * integral F exp(+jk*s.r) dS.
% Thus sigma = (4*pi/lambda^2)*|integral e_rx.F exp(+jk*s.r) dS|^2.
nr = size(E,1);
n = repmat(normal,nr,1);
s = repmat(observation,nr,1);
etaJ = cross(n,cross(direction,E,2),2);
F = (cross(s,cross(s,etaJ,2),2)-cross(s,cross(n,E,2),2))/2;
end
