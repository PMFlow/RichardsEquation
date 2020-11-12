function [p,Vx,Vy] = realiz_Gauss_GRW(Nmod,varK,phi,wavenum,KMean,I,J,dx,dy,x,y,x2,y2,p0)
%%   2D flow-solver for realizations of random velocity fields
close all

%%   Assigning Parameters
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
S=1e5; 
maxr=0.8;
Tolerance=1e-6; 

dt=8*(max(max(Dx))/dx^2+max(max(Dy))/dy^2); dt=maxr*1/dt
rx=dt*Dx/dx^2; ry=dt*Dy/dy^2;   % r<=1/4 such that 1-4*r>0 !!!
rloc=1-(rx(:,1:I-2)+rx(:,2:I-1)+ry(1:J-2,:)+ry(2:J-1,:));
%%%%% GRW
p = p0; pa=p;
pp=zeros(J,I); eps=zeros(1,S); 
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
    if eps(s) <= Tolerance
        fprintf('Number of iterations = %d Error = %e\n',s,eps(s)) ;
        break
    end
    pa=p;   
end
%%%% Velocities
Vx=-(p(2:J-1,3:I)-p(2:J-1,1:I-2)).*D/(2*dx);
Vy=-(p(3:J,2:I-1)-p(1:J-2,2:I-1)).*D/(2*dy);
