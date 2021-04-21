%%   Comparison GRW-Sander's solution [Watson et al., 1995, Case 1]
%%   ==============================================================
clear all;
close all
tic
%% Grid Initialization
I= 1001; 
a=0; b=10; 
dx=(b-a)/(I-1);
x=a:dx:b;
fileID = fopen('test_dt','w');
%% Soil parameters
D0=2.75862;
Ks=0.1;
v=0.85;
theta_res=0.06;
theta_sat=0.35;
Ksat=Ks/(theta_sat-theta_res);
q0=Ks/1.25; % Ks/q=1.25 in [Watson et al., 1995]
Q=q0/(theta_sat-theta_res); % [Watson et al., 1995, Eq. (12)]
%% Fujita model
DD = @(tht) D0./(1-v*tht).^2;
KK = @(tht) Ksat.*(1-v)*tht./(1-v*tht);
VV= @(tht) Ksat.*(1-v)./(1-v*tht).^2; % V=dK(tht)/d(tht)
%% Initial Conditions
tht=zeros(1,I); 
tht0=tht;
D=DD(tht);
K=KK(tht);
V=VV(tht);
Dfactor=1.2; 
%% Solution
T=0.3625;
t=0;
dt=Dfactor*(2*max(D)/dx^2); dt=1./dt;
fprintf(fileID,'t = %2.2e, dt = %2.2e\n',t,dt);
while t<=T
    t=t+dt;
    [tht]=BGRW_1D_rand(tht0,I,dx,dt,V,D,K,Q);
    tht0=tht;
    D=DD(tht0);
    K=KK(tht);
    V=VV(tht);
    dt=Dfactor*(2*max(D)/dx^2); dt=1./dt;
    fprintf(fileID,'t = %2.2e, dt = %2.2e\n',t,dt);
end
%% Results
theta=theta_res + (theta_sat-theta_res)*tht;
load thetaE
figure;
plot(theta,x-b); 
hold
plot(tE,xt,'*')
ylabel('$z$','Interpreter','latex');
xlabel('$\theta(z,t)$','Interpreter','latex'); 

fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;

tC = [theta(I) theta(I-0.2/dx) theta(I-0.4/dx) theta(I-0.6/dx) theta(I-0.8/dx) theta(I-1.0/dx) theta(I-2.0/dx)];
fprintf('tC = %0.2e %0.2e %0.2e %0.2e %0.2e %0.2e %0.2e\n', tC );
fprintf('(tC-tE)./tE = %0.2e %0.2e %0.2e %0.2e %0.2e %0.2e %0.2e\n', (tC-tE)./tE );

toc

%% tE = 0.0907    0.0861    0.0819    0.0782    0.0748    0.0719    0.0631
% xt =   0   -0.2000   -0.4000   -0.6000   -0.8000   -1.0000   -2.0000
%% The space step is : 1.00e-02 
% The time step is : 7.40e-06 
% tC = 9.29e-02 8.84e-02 8.42e-02 8.02e-02 7.66e-02 7.34e-02 6.35e-02
% (tC-tE)./tE = 2.53e-02 2.68e-02 2.70e-02 2.60e-02 2.40e-02 2.12e-02 6.40e-03
% Elapsed time is 1.749253 seconds.
