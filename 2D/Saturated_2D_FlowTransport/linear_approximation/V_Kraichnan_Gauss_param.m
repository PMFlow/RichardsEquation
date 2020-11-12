function [wavenum, phi, amplitude] = V_Kraichnan_Gauss_param(NMOD,varK,ZC1,ZC2,U_MEAN,lambda)
% global state;

sq2Pi=1/sqrt(2)/pi;
dummfac=sqrt(2*varK/NMOD)*U_MEAN;

wavenum(:,1)=sq2Pi*randn(NMOD,1)/ZC1; 
wavenum(:,2)=sq2Pi*randn(NMOD,1)/ZC2; 
norm=1./sqrt(wavenum(:,1).*wavenum(:,1)+wavenum(:,2).*wavenum(:,2));

phi = 2.*pi*rand(NMOD,1);
gauss_factor=exp(-lambda*lambda./norm/2);

amplitude(:,1)=(1-wavenum(:,1).^2.*norm.^2)*dummfac.*gauss_factor;
amplitude(:,2)=-wavenum(:,1).*wavenum(:,2).*norm.^2*dummfac.*gauss_factor;
