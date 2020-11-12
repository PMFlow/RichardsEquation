%% Flow-benchmark [Schneid, 2000; List&Radu, 2016]

clear all; close all
tic
itest=0; % 0=loam soil; 1=clay soil;
%%   Grid Initialization
I=21; J=31; 
a=0; b=2;
c=0; d=3;
dx=(b-a)/(I-1);
x=a:dx:b;
x2=(x(1:I-1)+x(2:I))/2;
dy=(d-c)/(J-1);
y=c:dy:d;
y2=(y(1:J-1)+y(2:J))/2;
%% Parameters
if itest==1
    disp('Beit Netofa clay')
    Ksat = 8.2*10^-4;
    theta_res=0.0;
    theta_sat=0.446;
    alpha=0.152;
    nGM=1.17;   
    T = 3; past=T/3;
    S= 500;
    dt=1/3
else
    disp('slit loam')
    Ksat = 4.96*10^-2;
    theta_res=0.131;
    theta_sat=0.396;
    alpha=0.423;
    nGM=2.06;
    T = 3/16; past=T/3;
    S= 500; 
    dt=1/48;
end
% dt satisfies r<=1 for both clay and loam problems
t1=T/3;
maxr=0.8;
Tolerance = 1e-5;
%% van Genuchten-Mualem parameter functions
ng=(nGM-1)/nGM;
theta =  @(p)  theta_GM(theta_res,theta_sat,p,alpha,nGM,ng);
dtheta = @(p) ng*(theta_sat-theta_res)*(1./(1+(-alpha*c).^nGM)).^(ng-1).*(-alpha*(-alpha*c).^(nGM-1)./(1+(-alpha*c).^nGM).^2);
K = @(tht) Ksat.*((tht - theta_res)./(theta_sat-theta_res)).^0.5.*(1-(1-((tht - theta_res)./(theta_sat-theta_res)).^(1/ng)).^ng).^2; %.^ng).^(1/ng)).^2; 
%% Initial Conditions 
p0 = 1 - ((y'-c)/(d-c))*3+0*x;
i0=1/dx; j0=1/dy;
p0(J,1:i0)=-2; 
p = p0;
%% Solution
tgraf=0; kt=1;
pp=zeros(J,I);
tht = theta(p);
tht0=tht;
Ltht_neg = dtheta(p);
Ltht_pos = zeros(J,I);
neg = p <0; 
Ltht = Ltht_neg.*neg + Ltht_pos.*~neg;
% L=max(max(abs(Ltht)))/2; % ~[List&Radu, 2016] =0.0032 (clay); =0.0109 (loam)
if itest==1
    L=0.12; % convergent for itest=1 (clay soil)
else
    L=0.5;  % convergent for itest=0 (loam soil)
end
pa=p;
iS=1:S;
t=0;
while t<=T
    t=t+dt;
    eps=zeros(1,S); 
    for s=1:S
        DK=K(tht);
        Dx=(DK(2:J-1,1:I-1)+DK(2:J-1,2:I))/2;
        Dy=(DK(1:J-1,2:I-1)+DK(2:J,2:I-1))/2;
        D=DK(2:J-1,2:I-1);
        rx=dt*Dx/dx^2/L; ry=dt*Dy/dy^2/L;
        rloc=1-(rx(:,1:I-2)+rx(:,2:I-1)+ry(1:J-2,:)+ry(2:J-1,:));
        pp(2:J-1,2:I-1)=rloc.*p(2:J-1,2:I-1) ...
            +rx(:,1:I-2).*p(2:J-1,1:I-2)+rx(:,2:I-1).*p(2:J-1,3:I) ...
            +ry(1:J-2,:).*p(1:J-2,2:I-1) +ry(2:J-1,:).*p(3:J,2:I-1);
        %% Boundary conditions
        %%%% BCYLeft/Right
        pp(:,1)=pp(:,2); % no flux
        pp(1:j0,I)=p0(1:j0,I); % \Gamma_{D2}
        pp(j0+1:J,I)=pp(j0+1:J,I-1); % no flux        
        %%%% BCXBottom/Upper
        pp(1,2:I-1)=pp(2,2:I-1)+dy; % no flux
        if t<=t1 % \Gamma_{D1}
            pp(J,1:i0)=p0(J,1:i0)+2.2*t/t1;
        else
            pp(J,1:i0)=0.2;
        end
        pp(J,i0+1:I-1)=pp(J-1,i0+1:I-1)-dy;  % no flux
        dtht=(tht0-tht)/L;
        %% Source term
        f=(ry(2:J-1,:)-ry(1:J-2,:))*dy + dtht(2:J-1,2:I-1); 
        pp(2:J-1,2:I-1)=pp(2:J-1,2:I-1)+f;
        p=pp;
        %% Convergence criterion
        tol_eps=dx*norm(p-pa)+norm(p-pa)/norm(p);
        if kt*past>=t && kt*past<t+dt && t<=T 
            eps(s)=tol_eps;
        end
        if tol_eps <= Tolerance
            break
        end        
        tht = theta(p);
        pa=p;
    end
    if  kt*past>=t && kt*past<t+dt && t<=T 
        tgraf=tgraf+1; 
        if itest==1
            rndt=kt*past; % days
        else
            rndt=kt*past*24; % hours
        end
        fprintf('kt*past= %d\n',rndt);
        str=['t=',num2str(rndt)];
        strvect(tgraf,1:length(str))=str;
        figure(3); box; hold all;
        P(tgraf)=plot(iS(1:2:end),eps(1:2:end));
        kt=kt+1;
    end
    tht = theta(p);
    tht0=tht;
end
%% Velocity components
Vx=-D.*((p(2:J-1,3:I)-p(2:J-1,1:I-2))/(2*dx));
Vy=-D.*((p(3:J,2:I-1)-p(1:J-2,2:I-1))/(2*dy)+1);

toc % loam: cca 0.65 sec; clay: cca. 1.15 sec (Laptop 16 GB, 2.11 GHz)

%% Plots
NameArray = {'Marker'}; ValueArray = {'o','+','x'}';
set(P,NameArray,ValueArray);
set(gca,'yscale','log'); box on;
xlabel('$s$','Interpreter','latex'); 
ylabel('$\|\psi^s - \psi^{s-1}\|$','Interpreter','latex');
legend(strvect); legend('boxoff'); ylim([0.99*Tolerance inf]);
levels=12;
plot_contours(I,J,x,y,p,tht,Vx,Vy,levels);
plot_fig(I,J,x,y,p,tht,Vx,Vy)
fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;
if itest==1
    save('sol_GRW_2D_clay','x','y','p','tht','Vx','Vy') 
else
    save('sol_GRW_2D_loam','x','y','p','tht','Vx','Vy') 
end
