%% Coupled flow & transport 1D: benchmark problems, laom and clay soils
%% Random GRW/BGRW flow/transport solvers

clear all; close all
tic
global state;
initstate; state;
itest=0; % 0=loam soil; 1=clay soil;
irand=1; % 0=constant Ksat; 1=random Ksat; 
icorr=0; % 0=Gaussian correlation; 1=exponential correlation;
%%   Grid Initialization
I=61; 
a=0; b=3;
dx=(b-a)/(I-1);
x=a:dx:b;
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
    dt=1/3;
    Lp=0.5; Lc=Lp; 
else
    disp('slit loam')
    Ksat = 4.96*10^-2;
    theta_res=0.131;
    theta_sat=0.396;
    alpha=0.423;
    nGM=2.06;
    dt=1/48;
    Lp=2; Lc=Lp;
end
% dt<min(dtp,dtc); dtp, dtc satisfy r<=1 in flow and transport solvers;
% see 'Richards_1D_code_verification' for determination of dtp and dtc
T = 3;
past=T/3;
t1=T/3;
D1=0.001; 
Tolerance = 1e-5;
S= 50000;
%% Initialization of random/deterministic hydraulic conductivity K
if irand==1
    Nmod = 100;
    ZC1 = 0.1; 
    ZC2 = 0.01; 
    varK = 0.5;
    if icorr==1
        [wavenum, phi] = Kraichnan_Exp_param(Nmod,ZC1,ZC2);
    else
        [wavenum, phi] = Kraichnan_Gauss_param(Nmod,ZC1,ZC2);
    end
    C1 = wavenum(:,1);
    C2 = wavenum(:,2);
    Ks= K_r(x(:)',Nmod,Ksat,varK,C1,C2,phi);
else
    Ks=Ksat*ones(1,I);
end
%% van Genuchten-Mualem parameter functions
ng=(nGM-1)/nGM;
theta =  @(p,c)  theta_GM(theta_res,theta_sat,p,c,alpha,nGM,ng,aGM,bGM);
K = @(tht) Ks.*((tht - theta_res)./(theta_sat-theta_res)).^0.5.*(1-(1-((tht - theta_res)./(theta_sat-theta_res)).^(1/ng)).^ng).^2;
%%   Initial Conditions
p0=1-((x-a)/(b-a))*3;
c0=(x-a)/(b-a)/1.2;
p = p0; pBC=p0;
c0(I)=1; 
c0(1)=0;
c=c0; cBC=c0;  
figure;
plot(p,x); 
ylabel('$z$','Interpreter','latex');
xlabel('$\psi(z,t=0)$','Interpreter','latex'); 
figure;
plot(c,x); 
ylabel('$z$','Interpreter','latex');
xlabel('$c(z,t=0)$','Interpreter','latex'); 
%% Solution
tgraf=0; t=0; kt=1; 
convf=zeros(3,S); convt=zeros(3,S); nn=zeros(I);
tht=theta(p,c); thtc=tht.*c;
tht0=tht; thtc0=thtc;
pa=p; ca=c;
dN=10^24; 
n0 = floor(dN*p); n0E=round(pBC*dN);
n=n0; 
restr=0; restsar1=0; restsarI=0; rest1=0; rest2=0; restf=0;
iS=1:S;
while t<=T
    DK=K(tht); 
    D=(DK(1:I-1)+DK(2:I))/2;
    Dq=DK(2:I-1);
    t=t+dt;
    eps=zeros(1,S); epsc=zeros(1,S);     
    for s=1:S
%% Flow step        
        DK=K(tht); 
        D=(DK(1:I-1)+DK(2:I))/2; 
        Dq=DK(2:I-1);
        r=dt*D/dx^2/Lp; 
        rloc=[1-2*r(1),1-(r(1:I-2)+r(2:I-1)),1-2*r(I-1)];
        rapr=[0.5,r(1:I-2)./(r(1:I-2)+r(2:I-1))];
        restr=rloc.*n+restr; nn=floor(restr); restr=restr-nn;  nsar=n-nn;
        restsar1=r(1)*n(1)+restsar1; nsar1=floor(restsar1); restsar1=restsar1-nsar1;
        nn(2)=nn(2)+nsar1;
        rest1=r(1:I-2).*n(2:I-1)+rest1; nsarleft=floor(rest1); rest1=rest1-nsarleft;
        nn(1:I-2)=nn(1:I-2)+nsarleft;
        nn(3:I)=nn(3:I)+nsar(2:I-1)-nsarleft; 
        restsarI=r(I-1)*n(I)+restsarI; nsarI=floor(restsarI); restsarI=restsarI-nsarI;
        nn(I-1)=nn(I-1)+nsarI;    
        %% Boundary conditions - pressure
        %%%% BC_Bottom/Top
        nn(1)=n0E(1); 
        if t<=t1 
            nn(I)=n0E(I)+2.2*t/t1; 
        else
            nn(I)=0.2; 
        end
        %% Source term - pressure
        dtht=(tht0-tht)/Lp;
        f=diff(r)*dx+dtht(2:I-1);
        restf=dN*f+restf; nf=floor(restf); restf=restf-nf;
        nn(2:I-1)=nn(2:I-1)+nf;
        n=nn; p=n/dN;
        thtc=theta(p,c).*c;  p0=p; % Alternating L-scheme      
%% Transport step
        dthtc=(thtc0-thtc)/Lc;
        q(2:I-1)=-Dq.*((p0(3:I)-p0(1:I-2))/(2*dx)+1); % flow-interior of \Omega
%         q(1)=q(2); q(I)=q(I-1); % flow BCs extended from interor of \Omega
        q(1)=-DK(1).*((p0(2)-p0(1))/dx+1); % flow left_BC approximated by finie-difference 
        q(I)=-DK(I).*((p0(I)-p0(I-1))/dx+1); % flow right_BC approximated by finie-difference 
        mean_q=mean(q)*ones(1,I);
        [c]=BGRW_1D_Alt_rand_benchmark(c0,cBC,dthtc,I,dx,dt,Lc,q,D1); % transport solver
        c0=c; tht=theta(p,c);
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
        test_kt=1;
    else
        test_kt=0;
    end
    tht0=tht; thtc0=tht.*c0; 
end
if  kt==4
    save('convf');
    save('convt');
end

toc % random Ksat, Gaussian corr., loam: 5.15 sec; clay: 0.65 sec; (Laptop 16 GB, 2.11 GHz)

%% Plots
kt_plot=kt-1;
plot_conv_Alt(kt_plot,S,strvect)
max_grid_Peclet_x=dx*max(max(abs(q)))/D1
mean_grid_Peclet_x=dx*mean(mean(abs(q)))/D1
plot_fig(x,p,c,tht,q)
fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;
