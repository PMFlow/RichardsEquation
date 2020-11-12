%% Computes ensembles of random velocity fields
%% N=10 random modes; R=100 realizations; Tolerance in GRW flow-solver = 10^-6 
clear ; close all ;

% global state;
% initstate; state;
Nmod = 10;
varK= 0.1 ;
NRealiz = 1; % 10^2 ; % 
ZC1 = 1.0;
ZC2 = 1.0;
KMean = 15;
%% Grid Initialization
I=201;
J=101;
a=0; b=20;
c=0; d=10;
dx=(b-a)/(I-1)
x=a:dx:b;
x2=(x(1:I-1)+x(2:I))/2;
dy=(d-c)/(J-1)
y=c:dy:d;
y2=(y(1:J-1)+y(2:J))/2;
p0 = 1 - ((x-a)/(b-a))+0*y';
%% GAUSSIAN CORRELATION
tic ;
for n = 1 : NRealiz
    [wavenum, phi] = Kraichnan_Gauss_param(Nmod,ZC1,ZC2);
    [p,Vx,Vy] = realiz_Gauss_GRW(Nmod,varK,phi,wavenum,KMean,I,J,dx,dy,x,y,x2,y2,p0) ;
    fprintf('Number of realization : %d \n',n) ;
%     save(['GAUSS','\dataGAUSS_',num2str(n),'.mat'],'p','Vx','Vy') ;
end
toc ; % 4.47 h / 100 realizations; % 160.74 sec / 1 realization; 
