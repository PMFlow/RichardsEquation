%%   GRW-Warrick's solution [Warrick et al. 1985]
%%   ============================================
clear all;
close all
tic
%% Grid Initialization
I= 101; 
a=0; b=1;
dx=(b-a)/(I-1);
x=a:dx:b;
%%   Parameters
Ksat = 0.5184; % m/d = 6e-4 cm/s [Warrick et al. 1985]
theta_res=0.1;
theta_sat=0.45;
alpha=1; % m^-1 = 0.01 cm^-1 [Warrick et al. 1985]
nGM=1.5;
Lp= 0.2; 
Tolerance = 1e-5; 
S= 50000;
maxr=1; 
%% van Genuchten-Mualem parameter functions
ng=(nGM-1)/nGM;
theta =  @(p) theta_GM(theta_res,theta_sat,p,alpha,nGM,ng);
K = @(tht) Ksat.*((tht - theta_res)./(theta_sat-theta_res)).^0.5.*(1-(1-((tht - theta_res)./(theta_sat-theta_res)).^(1/ng)).^ng).^2;
%%   Initial Conditions
% Tht=0.2; % Tht=(tht-theta_res)/(theta_sat-theta_res)
% p0W = ((Tht^(-1/ng)-1)^(1/nGM)); = theta(p)^-1 (van Genuchten-Mualem model)
p0W = -24.8665;
p0=p0W*ones(1,I); p0(I)=0;
p = p0; pBC=p0;
%% Solutions for T = 0.5, 1, 1.5, 2 hours
jT = 1; nt=4;
pT=zeros(nt,I); thT=zeros(nt,I); qT=zeros(nt,I);
fileID = fopen('test_dt','w');
for it = [0.0208    0.0417    0.0625    0.0833] % (days)
T=it;
past=T/3;
fprintf(fileID,'T= %2.2e (hours)\n',it*24);
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
        r=dt*D/dx^2/Lp; rq=dt*Dq/dx^2/Lp; 
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
        %%%% BC_Top/Bottom
        nn(I)=0; % infiltration
        nn(1)=n0E(1); 
        %% Source term - pressure
        dtht=(tht0-tht)/Lp;
        f=-diff(r)*dx+dtht(2:I-1); % -diff(r) because z is positive downward 
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
q=zeros(1,I);
        p0=p; 
        q(2:I-1)=-Dq.*((p0(3:I)-p0(1:I-2))/(2*dx)+1); % flow-interior of \Omega
%         q(1)=q(2); q(I)=q(I-1); % flow BCs extended from interor of \Omega
        q(1)=-DK(1).*((p0(2)-p0(1))/dx+1); % flow left_BC approximated by finie-difference 
        q(I)=-DK(I).*((p0(I)-p0(I-1))/dx+1); % flow right_BC approximated by finie-difference 
        mean_q=mean(q)*ones(1,I);
%% save solutions at t=T
rit=round(it*24*2)/2;
strT=['T=',num2str(rit)];
strvecT(jT,1:length(strT))=strT;
pT(jT,:)=p; thT(jT,:)=tht; qT(jT,:)=q;
% if  kt==4
%     save('pT','pT'); save('thT','thT'); save('qT','qT');
% end
jT = jT+1;
end

%% Plots
jT_plot=jT-1; xb=(x-b);
plot_T(jT_plot,xb,strvecT)
fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;

toc
