%% Coupled flow & transport 2D: benchmark problems, laom and clay soils
%% Deterministic-GRW flow solver; random-BGRW transport solver

clear all; close all
tic
global state;
initstate; state;
itest=0; % 0=loam soil; 1=clay soil;
irand=1; % 0=constant Ksat; 1=random Ksat;
icorr=0; % 0=Gaussian correlation; 1=exponential correlation; 
%%   Grid Initialization
I=41; J=61;
a=0; b=2;
c=0; d=3;
dx = (b-a)/(I-1);
x = a:dx:b;
x2 = (x(1:I-1)+x(2:I))/2;
dy=(d-c)/(J-1);
y=c:dy:d;
y2=(y(1:J-1)+y(2:J))/2;
%%   Parameters
aGM=0.44;
bGM=0.0046;
if itest==1
    disp('Beit Netofa clay')
    Ksat = 8.2*10^-4;
    theta_res=0.0;
    theta_sat=0.446;
    alpha=0.152;
    nGM=1.17;   
    T = 3;
    past=T/3;
    dt=1/3;
    Lp=100; Lc=Lp; 
else
    disp('slit loam')
    Ksat = 4.96*10^-2;
    theta_res=0.131;
    theta_sat=0.396;
    alpha=0.423;
    nGM=2.06;
    T = 3;
    past=T/3;
    dt=1/48;
    Lp=20; Lc=Lp;
end
% dt<min(dtp,dtc); dtp, dtc satisfy r<=1 in flow and transport solvers
t1=T/3;
maxr=0.8;
D1=0.001; D2=D1;
Tolerance = 1e-5;
S= 50000;
%% Initialization of random/determinist hydraulic conductivity K
if irand==1
    Nmod = 100;
    ZC1 = 0.1; 
    ZC2 = 0.01; 
    varK = 0.5;
    [X,Y] = meshgrid(x,y);
    if icorr==1
        [wavenum, phi] = Kraichnan_Exp_param(Nmod,ZC1,ZC2);
    else
        [wavenum, phi] = Kraichnan_Gauss_param(Nmod,ZC1,ZC2);
    end
    C1 = wavenum(:,1);
    C2 = wavenum(:,2);
    Ks= K_r(X(:)',Y(:)',Nmod,Ksat,varK,C1,C2,phi);
    Ks=reshape(Ks,J,I);
else
    Ks=Ksat*ones(J,I);
end
%% van Genuchten-Mualem parameter functions
ng=(nGM-1)/nGM;
theta =  @(p,c)  theta_GM(theta_res,theta_sat,p,c,alpha,nGM,ng,aGM,bGM);
dtheta = @(p) ng*(theta_sat-theta_res)*(1./(1+(-alpha*p).^nGM)).^(ng-1).*(-alpha*(-alpha*p).^(nGM-1)./(1+(-alpha*p).^nGM).^2);
K = @(tht) Ks.*((tht - theta_res)./(theta_sat-theta_res)).^0.5.*(1-(1-((tht - theta_res)./(theta_sat-theta_res)).^(1/ng)).^ng).^2; % 8,7 min. 
%%   Initial Conditions 
p0 = 1 - ((y'-c)/(d-c))*3+0*x;
c0 = (y'-c)/(d-c)/1.2+0*x;
i0=1/dx; j0=1/dy;
p0(J,1:i0)=-2; 
p = p0; pBC=p0;
c0(J,1:i0)=1; 
c0(1:j0,I)=zeros(j0,1);
c=c0; cBC=c0;  
figure;
mesh(x,y,p);
xlabel('x','Interpreter','latex'); ylabel('z','Interpreter','latex');
zlabel('$\psi(x,z,t=0)$','Interpreter','latex'); view(115,15);
figure;
mesh(x,y,c);
xlabel('x','Interpreter','latex'); ylabel('z','Interpreter','latex');
zlabel('$c(x,z,t=0)$','Interpreter','latex'); view(115,15);
%% Solution
tgraf=0; t=0; kt=1; 
Vx=zeros(J,I); Vy=zeros(J,I);
convf=zeros(3,S); convt=zeros(3,S); pp=zeros(J,I);
tht = theta(p,c);
tht0=tht; thtc0=tht.*c0; 
pa=p; ca=c;
% Dfactor=1.2;
% dtc=Dfactor*(2*D1/Lc/dx^2+2*D2/Lc/dy^2); dtc=1./dtc % BGRW
% dt=min(dtc,dt)
iS=1:S;
while t<=T
    t=t+dt;
    eps=zeros(1,S); epsc=zeros(1,S);     
    for s=1:S
%% Flow step        
        DK=K(tht);
        Dx=(DK(2:J-1,1:I-1)+DK(2:J-1,2:I))/2;
        Dy=(DK(1:J-1,2:I-1)+DK(2:J,2:I-1))/2;
        D=DK(2:J-1,2:I-1);
        rx=dt*Dx/dx^2/Lp; ry=dt*Dy/dy^2/Lp; % Eq. (14) ---> r<=1/4 such that 1-4*r>=0
        rloc=1-(rx(:,1:I-2)+rx(:,2:I-1)+ry(1:J-2,:)+ry(2:J-1,:));
        pp(2:J-1,2:I-1)=rloc.*p(2:J-1,2:I-1) ...
            +rx(:,1:I-2).*p(2:J-1,1:I-2)+rx(:,2:I-1).*p(2:J-1,3:I) ...
            +ry(1:J-2,:).*p(1:J-2,2:I-1) +ry(2:J-1,:).*p(3:J,2:I-1);       
        %% Boundary conditions - pressure
        %%%% BCYLeft/Right
        pp(:,1)=pp(:,2); % no flux
        pp(1:j0,I)=pBC(1:j0,I); % \Gamma_{D2}
        pp(j0+1:J,I)=pp(j0+1:J,I-1); % no flux       
        %%%% BCXBottom/Upper
        pp(1,2:I-1)=pp(2,2:I-1)+dy; % no flux
        if t<=t1 % \Gamma_{D1}
            pp(J,1:i0)=pBC(J,1:i0)+2.2*t/t1;
        else
            pp(J,1:i0)=0.2;
        end
        pp(J,i0+1:I-1)=pp(J-1,i0+1:I-1)-dy;  % no flux
        %% Source term - pressure
        dtht=(tht0-tht)/Lp;
        f=(ry(2:J-1,:)-ry(1:J-2,:))*dy + dtht(2:J-1,2:I-1); 
        pp(2:J-1,2:I-1)=pp(2:J-1,2:I-1)+f;
        p=pp;
        tht = theta(p,c);
%         pa=p; 
        thtc=tht.*c;  p0=p; %% Alt-scheme
%% Transport step
        dthtc=(thtc0-thtc)/Lc;
%       V_interior of \Omega        
        Vx(2:J-1,2:I-1)=-D.*((p0(2:J-1,3:I)-p0(2:J-1,1:I-2))/(2*dx)); 
        Vy(2:J-1,2:I-1)=-D.*((p0(3:J,2:I-1)-p0(1:J-2,2:I-1))/(2*dy)+1); 
%       V_normal to boundary approximated by finie-differences
        Vx(:,1)=-DK(:,1).*(p0(:,2)-p0(:,1))/dx; Vx(:,I)=-DK(:,I).*(p0(:,I)-p0(:,I-1))/dx; 
        Vx(1,2:I-1)=-DK(1,2:I-1).*(p0(1,3:I)-p0(1,1:I-2))/(2*dy); 
        Vx(1,1)=-DK(1,1).*(p0(1,2)-p0(1,1))/dy; 
        Vx(1,I)=-DK(1,I-1).*(p0(1,I)-p0(1,I-1))/dy; 
        Vx(J,2:I-1)=-DK(J,2:I-1).*(p0(J,3:I)-p0(J,1:I-2))/(2*dy); 
        Vx(J,1)=-DK(J,1).*(p0(J,2)-p0(J,1))/dy; 
        Vx(J,I)=-DK(J,I-1).*(p0(J,I)-p0(J,I-1))/dy; 
        Vy(1,:)=-DK(1,:).*((p0(2,:)-p0(1,:))/dy+1); Vy(J,:)=-DK(J,:).*((p0(J,:)-p0(J-1,:))/dy+1);
        Vy(2:J-1,1)=-DK(2:J-1,1).*((p0(3:J,1)-p0(1:J-2,1))/(2*dy)+1);
        Vy(1,1)=-DK(1,1).*((p0(2,1)-p0(1,1))/dy+1);
        Vy(J,1)=-DK(J-1,1).*((p0(J,1)-p0(J-1,1))/dy+1);
        Vy(2:J-1,I)=-DK(2:J-1,I).*((p0(3:J,I)-p0(1:J-2,I))/(2*dy)+1);
        Vy(1,I)=-DK(1,I).*((p0(2,I)-p0(1,I))/dy+1);
        Vy(J,I)=-DK(J-1,I).*((p0(J,I)-p0(J-1,I))/dy+1);
%       V_normal to boundary extended from interior of \Omega
%         Vx(:,1)=Vx(:,2); Vx(:,I)=Vx(:,I-1); 
%         Vx(1,:)=Vx(2,:); Vx(J,:)=Vx(J-1,:); 
%         Vy(:,1)=Vy(:,2); Vy(:,I)=Vy(:,I-1);
%         Vy(1,:)=Vy(2,:); Vy(J,:)=Vy(J-1,:); 
%       Transport solver        
        [c]=BGRW_Alt_rand_benchmark(c0,cBC,dthtc,I,J,dx,dy,i0,j0,dt,Lc,Vx,Vy,D1,D2);
        c0=c; 
        tht = theta(p,c);
        thtc=tht.*c;  
%% Convergence criterion
        tol_eps=dx*norm(p-pa)+norm(p-pa)/norm(p);
        tol_epsc=dx*norm(c-ca)+norm(c-ca)/norm(c);
        if kt*past>=t && kt*past<t+dt && t<=T
            eps(s)=tol_eps;
            epsc(s)=tol_epsc;
        end
        if max(tol_eps,tol_epsc) <= Tolerance
            break
        end
        pa=p; ca=c;       
    end
    if  kt*past>=t && kt*past<t+dt && t<=T
        tgraf=tgraf+1;
        rndt=kt*past; % time in days
        str=['t=',num2str(rndt)];
        strvect(tgraf,1:length(str))=str;
        convf(kt,:)=eps;
        convt(kt,:)=epsc;
        kt=kt+1;
    end
    tht = theta(p,c);
    tht0=tht; thtc0=tht.*c0; 
end
if  kt==4
    if itest == 0
        save('convf_loam');
        save('convt_loam');
    else
        save('convf_clay');
        save('convt_clay');
    end
end

toc % random Ksat, Gaussian corr., loam: cca 7.6 min; clay: cca 2.4 min; (Laptop 16 GB, 2.11 GHz)

%% Results
plot_fig(x,y,p,c,tht,Vx,Vy)
if irand==0
    if itest==0
        save('solution_loam_det','x','y','p','c','tht','Vx','Vy','convf','convt','DK');
    else
        save('solution_clay_det','x','y','p','c','tht','Vx','Vy','convf','convt','DK');
    end
else
    if itest==0
        save('solution_loam','x','y','p','c','tht','Vx','Vy','convf','convt','DK');
    else
        save('solution_clay','x','y','p','c','tht','Vx','Vy','convf','convt','DK');
    end
end
fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;

%% mean and max Peclet numbers
V=mean(mean(sqrt(Vx.^2+Vy.^2)))
Pe=V*dx/D1
maxV=max(max(sqrt(Vx.^2+Vy.^2)))
maxPe=maxV*dx/D1
