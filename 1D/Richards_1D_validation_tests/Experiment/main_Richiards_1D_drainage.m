%% GRW solution for the free drainage problem (Hydrus 1D example)
%% ==============================================================
clear all; 
close all
tic
%%   Grid Initialization
I= 61; 
a=0; b=6; 
dx=(b-a)/(I-1);
x=a:dx:b;
%%   Parameters
Ksat = 0.25; 
theta_res=0.0; 
theta_sat=0.331;
alpha=1.43; 
nGM=1.5;
dt=1/3;
Lp= 0.5; 
Tolerance = 1e-5; 
S= 50000;
maxr=1; 
%% van Genuchten-Mualem parameter functions
ng=(nGM-1)/nGM;
theta =  @(p) theta_GM(theta_res,theta_sat,p,alpha,nGM,ng);
K = @(tht) Ksat.*((tht - theta_res)./(theta_sat-theta_res)).^0.5.*(1-(1-((tht - theta_res)./(theta_sat-theta_res)).^(1/ng)).^ng).^2;
%%   Initial Conditions
p0=zeros(1,I); 
p = p0; pBC=p0;

%% Solutions for T = 1, 4, 20, and 100 days
jT = 1;
pT=zeros(4,I); thT=zeros(4,I); qT=zeros(4,I);
fileID = fopen('test_dt','w');
for it = [1 4 20 100]
T = it; 
past=T/3;
fprintf(fileID,'T= %2.2e (days)\n',it);
%% Solution for t\in[0,T]
tgraf=0; t=0; kt=1; 
convf=zeros(3,S); convt=zeros(3,S); nn=zeros(I);
tht=theta(p);
tht0=tht;
pa=p;
dN=10^24; 
n0 = floor(dN*p); n0E=round(pBC*dN);
n=n0; 
restr=0; restsar1=0; restsarI=0; rest1=0; rest2=0; restf=0;
iS=1:S;
while t<=T
    DK=K(tht); 
    D=(DK(1:I-1)+DK(2:I))/2;
    Dq=DK(2:I-1);
    dt=Lp*dx^2*maxr/max(D)/2; % r<=1/2 such that 1-2*r>0;
    t=t+dt;
    eps=zeros(1,S); epsc=zeros(1,S);     
    fprintf(fileID,'t = %2.2e, dt = %2.2e\n',t,dt);
    for s=1:S
%% Flow step        
        DK=K(tht); 
        D=(DK(1:I-1)+DK(2:I))/2; 
        Dq=DK(2:I-1);
        dt=Lp*dx^2*maxr/max(D)/2; % r<=1/2 such that 1-2*r>0;
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
        %%%% BC_Top
        nn(I)=nn(I-1)-dx*dN; % free drainage (zero flux on top)
        %% Source term - pressure
        dtht=(tht0-tht)/Lp;
        f=diff(r)*dx+dtht(2:I-1);
        restf=dN*f+restf; nf=floor(restf); restf=restf-nf;
        nn(2:I-1)=nn(2:I-1)+nf;
        n=nn; p=n/dN;
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
        rndt=kt*past; % time in days
        str=['t=',num2str(rndt)];
        strvect(tgraf,1:length(str))=str;
        convf(kt,:)=eps;
        kt=kt+1;
        test_kt=1;
    else
        test_kt=0;
    end
    tht0=tht;
end
if  kt==4
    save('convf');
end
kt_plot=kt-1; 
plot_conv(kt_plot,S,strvect)

%% Velocity
        p0=p; 
        q(2:I-1)=-Dq.*((p0(3:I)-p0(1:I-2))/(2*dx)+1); % flow-interior of \Omega
%         q(1)=q(2); q(I)=q(I-1); % flow BCs extended from interor of \Omega
        q(1)=-DK(1).*((p0(2)-p0(1))/dx+1); % flow left_BC approximated by finie-difference 
        q(I)=-DK(I).*((p0(I)-p0(I-1))/dx+1); % flow right_BC approximated by finie-difference 
        mean_q=mean(q)*ones(1,I);

%%
strT=['T=',num2str(it)];
strvecT(jT,1:length(str))=strT;
pT(jT,:)=p; thT(jT,:)=tht; qT(jT,:)=q;
if  kt==4
    save('pT','pT'); save('thT','thT'); save('qT','qT');
end
jT = jT+1;
end

%% Plots
kt_plot=kt-1; jT_plot=jT-1; xb=(x-b)*100;
plot_T(jT_plot,xb,strvecT)

fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;
toc 
