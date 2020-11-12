function [wavenum, phi] = Kraichnan_Exp_param_short(NMOD,ZC1,ZC2)

global state; 

gamma_1=unifrnd(0,1,NMOD,1); 
gamma_2=unifrnd(0,1,NMOD,1);
wavenum(:,1)=sqrt(1./gamma_2.^2-1).*cos(2*pi*gamma_1)...
/(2*pi*ZC1);
wavenum(:,2)=sqrt(1./gamma_2.^2-1).*sin(2*pi*gamma_1)...
/(2*pi*ZC2); 

phi = 2.*pi*rand(NMOD,1);
