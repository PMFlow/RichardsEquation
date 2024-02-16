%%   FABRICATED Initial Conditions for the benchmark problem [List&Radu, 2016]
%%   Figs. 1 and 2 in [Song & Jiang, 2023, https://checlams.github.io/assets/pdf/paper/songescape2023.pdf]
%%   =====================================================================================================

%%   Clearing Screen and Variables
clear all; close all
tic
%%   Grid Initialization
I = 21; J=21; % ~ [Song & Jiang, 2023]
a=0; b=2;
c=0; d=2;
dx = (b-a)/(I-1);
x = a:dx:b;
dy=(d-c)/(J-1);
y=c:dy:d;
%%   Parameters
    % slit loam
    Ksat = 4.96*10^-2;
    theta_res=0.131;
    theta_sat=0.396;
    alpha=0.423;
    nGM=2.06;
    T = 3/16; past=T/3;
    S= 500; 
    dt=1/48;
t1=T/3;
maxr=0.47; 
Tolerance = 1e-5; 
%%  Parametrization functions
%% van Genuchten-Mualem
ng=(nGM-1)/nGM;
theta =  @(c)  theta_GM(theta_res,theta_sat,c,alpha,nGM,ng);
dtheta = @(c) ng*(theta_sat-theta_res)*(1./(1+(-alpha*c).^nGM)).^(ng-1).*(-alpha*(-alpha*c).^(nGM-1)./(1+(-alpha*c).^nGM).^2);
K = @(tht) Ksat.*((tht - theta_res)./(theta_sat-theta_res)).^0.5.*(1-(1-((tht - theta_res)./(theta_sat-theta_res)).^(1/ng)).^ng).^2; %.^ng).^(1/ng)).^2; 
%%   Initial Conditions 
c0 = 1 - ((y'-c)/(d-c))*3+0*x; 
i0=1/dx; j0=1/dy; 
c0(J,1:i0)=-2; % original benchmark !!!
%% FABRICATION of the Initial conditions:
c0(:,I)=0; c0(:,1)=c0(:,1)+0.1; c0(J,I)=c0(J,I)+0.1; % used to FABRICATE Figs. 1 and 2 in [Song & Jiang,2023];
c0(J,1:i0)=0; % used to obtain Figs. 1b and 2b in [Song & Jiang,2023]; when commented ---> Figs. 1a and 2a
%%
c = c0; tht = theta(c);
levels=12; 
figure;
contourf(x,y,c,levels); colormap(flipud(parula)); colorbar;
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
title('$\psi(x,z,0)$','Interpreter','latex'); 
figure;
contourf(x,y,tht,levels); colormap(flipud(parula)); colorbar;
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
title('$\theta(x,z,0)$','Interpreter','latex'); 
figure;
mesh(x,y,c); 
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
zlabel('$\psi(x,z,0)$','Interpreter','latex'); view(115,15); 
grid on ;
figure;
mesh(x,y,tht); 
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
zlabel('$\theta(x,z,0)$','Interpreter','latex'); view(115,15); 
grid on ;

toc