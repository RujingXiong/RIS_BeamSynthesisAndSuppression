% omega is the phase of vector w, i.e. omega = [0,pi,pi/2,2*pi,pi,...]
function [Fx,DeltaF] = Quadratic(h,omega)
% w = exp(1i*omega);
% hw = h'*w;
% Fx = norm(hw)^2;
% % DeltaF = 2*1i*R*w.*exp(1i.*omega);
% DeltaF = 2*imag(conj(w).*h*(hw));

w = exp(1i*omega);
Fx = real(w'*h*w);
% DeltaF = 2*1i*R*w.*exp(1i.*omega);
DeltaF = 2*imag(conj(w).*h*w);