function [u, v] = V_Kraichnan_Gauss_func(x,y,wavenum, phi, amplitude)

phase=cos(2*pi*(wavenum(:,1).*x+wavenum(:,2).*y)+phi);
u=sum(amplitude(:,1).*phase);
v=sum(amplitude(:,2).*phase);
