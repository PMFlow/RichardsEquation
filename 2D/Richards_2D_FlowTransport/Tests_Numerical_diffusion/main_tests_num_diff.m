%% test for numerical-diffusion (2D-case, constant flow velocity and diffusion coeffiicent)
%% transport step solved with either biased (BGRW) and unbiased (GRW) algorithms (lines 94-95)

close all; clear all; 
tic
%% uncomment one of the following to estimate the numerical diffusion for different Peclet numbers:
% with BGRW
%  I=21; J=31;  % TR = 2;   eps_c = 0.5300;     eps_D1 = 0.0755;     eps_D2 = 0.2596 (Pe = 3.3067)
% I=41; J=61;   % TR = 9;   eps_c = 0.1076;     eps_D1 = 1.8952e-16; eps_D2 = 1.4808e-15 (Pe = 1.6533) % optimal for the benchmark problem
% I=201; J=301; % TR = 239; eps_c = 0.0027;     eps_D1 = 4.1613e-16; eps_D2 = 1.0183e-15 (Pe = 0.3307)
% I=401; J=601; % TR = 960; eps_c = 5.7383e-04; eps_D1 = 2.9309e-15; eps_D2 = 3.6296e-15 (Pe = 0.1653)
% with unbiased GRW
% I=21; J=31;   % TR = 2;   eps_c = 1.0739; eps_D1 = 1.9423e-16; eps_D2 = 6.1421e-16 (Pe = 3.3067)
% I=41; J=61;   % TR = 4;   eps_c = 0.2038; eps_D1 = 6.5977e-17; eps_D2 = 8.0535e-16 (Pe = 1.6533)
% I=201; J=301; % TR = 19;  eps_c = 0.0538; eps_D1 = 1.9371e-16; eps_D2 = 4.7946e-16 (Pe = 0.3307)
I=401; J=601; % TR = 39;  eps_c = 0.0746; eps_D1 = 2.0965e-15; eps_D2 = 8.9192e-16 (Pe = 0.1653)
%%   Grid Initialization
x1=0; x2=2;
y1=0; y2=3;
dx = (x2-x1)/(I-1);
dy=(y2-y1)/(J-1);
x=x1:dx:x2; y=y1:dy:y2;
%% Parameters
K_sat=4.96*10^-2; % laom case in benchmark [List&Radu, 2016] ===> U_MEAN=-0.0331
theta_sat=1;
D1=0.001; D2=D1;
S=3;
L=1; 
maxr=0.8;
Tolerance = 1e-5;
%% Parametr functions
vectheta_sat=theta_sat*ones(J,I);
vectK_sat=K_sat*ones(J,I);
theta=@(p) vectheta_sat; 
K=@(tht) vectK_sat; 
%%  Initial Conditions
p0=1-((y'-y1)/(y2-y1))+0*x; 
p = p0; pa=p;
figure;
mesh(x,y,p);
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
zlabel('$\psi(x,z,t=0)$','Interpreter','latex'); view(115,15);
%% solution
eps=zeros(1,S);
tht=theta(p);
tht0=tht;
DK=K(tht); % time independent hydraulic conductivity
Dx=(DK(2:J-1,1:I-1)+DK(2:J-1,2:I))/2;
Dy=(DK(1:J-1,2:I-1)+DK(2:J,2:I-1))/2;
D=DK(2:J-1,2:I-1);
dt=2*(max(max(Dx))/(L*dx)^2+max(max(Dy))/(L*dy)^2); dt=maxr*1/dt;
rx=dt*Dx/dx^2/L; ry=dt*Dy/dy^2/L; % Eq. (14) ---> r<=1/4 such that 1-4*r>=0
rloc=1-(rx(:,1:I-2)+rx(:,2:I-1)+ry(1:J-2,:)+ry(2:J-1,:));
T=dt; % T is arbitrary for steady-state solutions
t=dt;
while t<=T
    t=t+dt;
    for s=1:S
        pp(2:J-1,2:I-1)=rloc.*p(2:J-1,2:I-1) ...
            +rx(:,1:I-2).*p(2:J-1,1:I-2)+rx(:,2:I-1).*p(2:J-1,3:I) ...
            +ry(1:J-2,:).*p(1:J-2,2:I-1) +ry(2:J-1,:).*p(3:J,2:I-1);
        %% Boundary conditions - pressure
        %%%% BCYLeft/Right
        pp(:,1)=pp(:,2); % von Neumann
        pp(:,I)=pp(:,I-1); %-qx=-D*((p(I)-p(I-1))/dx=0;
        %%%% BCXBottom/Upper
        pp(1,:)=p0(1,:); % Dirichlet
        pp(J,:)=p0(J,:);
        %% Source term - pressure
        dtht=(tht0-tht)/L; % vanishes for constant theta
        f=(ry(2:J-1,:)-ry(1:J-2,:))*dy + dtht(2:J-1,2:I-1);
        pp(2:J-1,2:I-1)=pp(2:J-1,2:I-1)+f;
        p=pp;
        %% Convergence criterion
        tol_eps=norm(p-pa)/norm(p);
        if tol_eps <= Tolerance
            fprintf('Number of iterations: %d \n',s) ;
            break
        end
        tht=theta(p);
        pa=p;
    end
    tht=theta(p);
    tht0=tht;
end
%% Velocity components
Vx=-D.*((p(2:J-1,3:I)-p(2:J-1,1:I-2))/(2*dx));
Vy=-D.*((p(3:J,2:I-1)-p(1:J-2,2:I-1))/(2*dy)+1);
U_MEAN=mean(mean(Vy));
grid_Peclet=dx*max(max(abs(Vy)))/D1;
fprintf('The constant velcity is : %0.2e \n',U_MEAN) ;
fprintf('The grid_Peclet number is : %0.2e \n',grid_Peclet) ;
%% Transport step
[c]=testBGRW_const(I,J,x,y,dx,dy,Vx,Vy,D1,D2,U_MEAN);
% [c]=testGRW_const(I,J,x,y,dx,dy,D1,D2,U_MEAN);

fprintf('The flow time step is : %0.2e \n',dt) ;
fprintf('The flow&transport space step is : %0.2e \n',dx) ;
toc
