%%   convergence test 2D_GRW-solver
clear ; close all ;

tic
global state;
initstate; state;
Nmod = 10;
KMean = 15; varK= 0.1 ;
ZC1 = 1.0; ZC2 = 1.0;
%%   Grid Initialization
I=201; J=101;
a=0; b=20;
c=0; d=10;
dx=(b-a)/(I-1);
x=a:dx:b;
x2=(x(1:I-1)+x(2:I))/2;
dy=(d-c)/(J-1);
y=c:dy:d;
y2=(y(1:J-1)+y(2:J))/2;
p0 = 1 - ((x-a)/(b-a))+0*y';
%%   Assigning Parameters
[wavenum, phi] = Kraichnan_Gauss_param(Nmod,ZC1,ZC2);
C1 = wavenum(:,1);
C2 = wavenum(:,2);
%%   Boundary and Initial Conditions
%homogeneous:
BCXL = @(v) 1 ;
BCXR = @(v) 0 ;
BCYB = @(u) 0*u ;
BCYU = @(u) 0*u ;
%%  Coefficients and Forcing Term
[X,Y] = meshgrid(x2(1:I-1),y(2:J-1));
Dx = K_repmat(X(:)',Y(:)',Nmod,KMean,varK,C1,C2,phi);
Dx = reshape(Dx,J-2,I-1);
[X,Y] = meshgrid(x(2:I-1),y2(1:J-1));
Dy = K_repmat(X(:)',Y(:)',Nmod,KMean,varK,C1,C2,phi);
Dy = reshape(Dy,J-1,I-2);
[X,Y] = meshgrid(x(2:I-1),y(2:J-1));
D = K_repmat(X(:)',Y(:)',Nmod,KMean,varK,C1,C2,phi);
D = reshape(D,J-2,I-2);
%% solution
%%%% Parameters
S=1e6; pass=1e3; % 1e5; pass=1000; % 
maxr=0.8;
Tolerance = 1e-5; % 1e-6; % 1e-7; % 1e-9; %
dt=8*(max(max(Dx))/dx^2+max(max(Dy))/dy^2); dt=maxr*1/dt
rx=dt*Dx/dx^2; ry=dt*Dy/dy^2;   % r<=1/4 such that 1-4*r>0 !!!
rloc=1-(rx(:,1:I-2)+rx(:,2:I-1)+ry(1:J-2,:)+ry(2:J-1,:));
%%%%% GRW
iS=1:S;
sump=zeros(1,S); diffL2=zeros(1,S); eps=zeros(1,S); pp=zeros(J,I);
p = p0; pa=p; sumpinit=sum(sum(p0));
for s=1:S
    pp(2:J-1,2:I-1)=rloc.*p(2:J-1,2:I-1) ...
        +rx(:,1:I-2).*p(2:J-1,1:I-2)+rx(:,2:I-1).*p(2:J-1,3:I) ...
        +ry(1:J-2,:).*p(1:J-2,2:I-1)+ry(2:J-1,:).*p(3:J,2:I-1);
    %%%% boundary conditions
    %%%% BCXL/R
    pp(:,1)=p0(:,1); pp(:,I)=p0(:,I);
    %%%% BCYB/U
    derB=BCYB(x(2:I-1));
    pp(1,2:I-1)=pp(2,2:I-1)-derB*dy;
    derU=BCYU(x(2:I-1));
    pp(J,2:I-1)=pp(J-1,2:I-1)+derU*dy;
    p=pp;
    eps(s)=dx*norm(p-pa)+norm(p-pa)/norm(p);
    pa=p;
    sump(s)=sum(sum(p));
    diffL2(s) = ( dx * dy )^(1/2) * norm(p-p0) ;
    if mod(s,pass)==0
        fprintf('s= %d\n',s);
    end
    if eps(s) <= Tolerance
        fprintf('Number of iterations: %d \n',s) ;
        iS=1:s;
        break
    end
end
%%%% Velocities
Vx=-(p(2:J-1,3:I)-p(2:J-1,1:I-2)).*D/(2*dx);
Vy=-(p(3:J,2:I-1)-p(1:J-2,2:I-1)).*D/(2*dy);

%% Plots
figure
plot(iS,eps(1:max(iS)),'r');
set(gca,'yscale','log'); 
xlabel('$s$','Interpreter','latex'); ylabel('$\|h^s - h^{s-1}\|$','Interpreter','latex');
% Constant nr. of particles indicates that the stationary state of h(t,x,y)is attained
figure;
plot(iS,sump(1:max(iS)),'o'); hold; plot(0,sumpinit,'pr')
set(gca,'yscale','log'); 
xlabel('t'); ylabel('$ \sum n $','Interpreter','latex');
figure;
plot(iS,diffL2(1:max(iS)),'b'); 
set(gca,'yscale','log'); 
xlabel('t'); ylabel('$\| h=h_0 \|_{L_2}$','Interpreter','latex')
%Velocity components
figure;
plot(x(2:I-1),Vx,x(2:I-1),Vy);
xlabel('x','Interpreter','latex'); ylabel('Velocity components')
%  Stationary numerical solution h(x,y), V_x(x,y), V_y(x,y)
figure;
mesh(x,y,p-p0); view(20,50);
xlabel('x','Interpreter','latex'); ylabel('y','Interpreter','latex'); zlabel('$h(x,y)-h_0(x,y)$','Interpreter','latex')
grid on ;
figure;
mesh(x(2:I-1),y(2:J-1),Vx); view(20,50);
xlabel('x','Interpreter','latex'); ylabel('y','Interpreter','latex'); zlabel('$V_x(x,y)$','Interpreter','latex')
figure;
mesh(x(2:I-1),y(2:J-1),Vy); view(20,50);
xlabel('x','Interpreter','latex'); ylabel('y','Interpreter','latex'); zlabel('$V_y(x,y)$','Interpreter','latex')

toc