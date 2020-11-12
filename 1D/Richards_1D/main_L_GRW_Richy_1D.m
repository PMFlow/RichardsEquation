%%   Richards-1D: two scenarios with transition unsaturate/saturated flow
close all
tic
itest=0; % 0=scenario1; 1=scenario2;
%%   Grid Initialization
I = 101;
a = 0 ;
b = 2 ;
dx = (b-a)/(I-1)
x = [a:dx:b];
x2 = (x(1:I-1)+x(2:I))/2;
T = 1*10^4; past=T/5; 
S = 20000; 
Tolerance=1e-9; 
%%   Parameters
Ksat = 2.77*10^-6;
theta_res=0.06;
theta_sat=0.36;
alpha=10;
q0=2.77*10^-7;
q1=2.50*10^-6;
t1=T/100;
maxr=0.8;
%%  Coefficients
vectKsat=Ksat*ones(1,I);
if itest == 1
    vectKsat((I-1)/2+1:I)=500*Ksat;
end
theta = @(p) theta_exp(theta_res,theta_sat,alpha,p);
K = @(theta) vectKsat.*(theta-theta_res)/(theta_sat-theta_res); 
%% Initial conditions
if itest == 0
    load IC_GRW_scenario1; 
else
    load IC_GRW_scenario2; 
end
Ltht=alpha*(theta_sat-theta_res)*exp(alpha*max(p));
% L=Ltht; % = 445.2395 (for both itest=1 and itest=0)
L=2;% for L=1.96  both scenario1 and scenario2 are already convergent
dN=10^10;
n0 = floor(dN*p);
p0=n0/dN;
n=n0; nA=n(1);
pinit=p;
thtinit=tht;
figure(1); plot(pinit,x,'*');
ylabel('$z$','Interpreter','latex'); xlabel('$\psi(z,t=0)$','Interpreter','latex');
figure(2); plot(thtinit,x,'*');
ylabel('$z$','Interpreter','latex'); xlabel('$\theta(z,t=0)$','Interpreter','latex');
%% Solution
tgraf=0; tconv=0; t=0; kt=1;
restr=0; restsar1=0; restsarI=0; rest1=0; rest2=0; restf=0;
tht=theta(p);
tht0=tht;
pa=p;
iS=1:S;
if itest == 0
    fileID = fopen('dt_scenario1','w');
else
    fileID = fopen('dt_scenario2','w');
end
while t<=T
    D=K(tht); D=D(1:I-1);
    dt=dx^2*maxr/max(D)/2; % r<=1/2 such that 1-2*r>=0
    t=t+dt;
    fprintf(fileID,'%2.2e  %2.2e\n',t, dt ) ;
    eps=zeros(1,S); 
    for s=1:S
        D=K(tht); D=D(1:I-1); 
        r=dt*D/dx^2/L;
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
        %% Dirichlet BC
        nn(1)=nA;
        %% Neuman BC
        if t<=t1
            qR=q0+t*(q1-q0)/t1; % -qR=-D(I-1)*((c(I)-c(I-1))/dx+1);
        else
            qR=q1;
        end
        nn(I)=nn(I-1)+(qR/D(I-1)-1)*dx*dN;
        dtht=(tht0-tht)/L;
        %% Source term
        f=diff(r)*dx+dtht(2:I-1); 
        restf=dN*f+restf; nf=floor(restf); restf=restf-nf;
        nn(2:I-1)=nn(2:I-1)+nf;
        n=nn; p=n/dN;
        na=n;
        %% Convergence criterion
        tol_eps=norm(p-pa)/norm(p);
        if kt*past>=t && kt*past<t+dt 
            eps(s)=tol_eps;
        end
        if tol_eps <= Tolerance
            break
        end        
        tht=theta(p);       
        pa=p;
    end
    %%
    if  kt*past>=t && kt*past<t+dt 
        rndt=kt*past;
        fprintf('t= %d\n',rndt);
        str=['t=',num2str(rndt)];
        strvect(kt,1:length(str))=str;
        figure(3); box; hold all;
        P(kt)=plot(iS,eps); 
        kt=kt+1;
    end
    tht=theta(p);    
    tht0=tht;
end
q=-D.*((p(2:I)-p(1:I-1))/dx+1);
%% Plots
NameArray = {'Marker'}; ValueArray = {'o','+','x','*','.'}';
set(P,NameArray,ValueArray); 
set(gca,'yscale','log'); %box on;
xlabel('$s$','Interpreter','latex');
ylabel('$\|\psi^s - \psi^{s-1}\| \;/\; \|\psi^s\|$','Interpreter','latex');
legend(strvect); legend('boxoff'); ylim([Tolerance inf]);
figure(4);
plot(p,x,'mo--','LineWidth',1.2,'MarkerFaceColor','yellow','MarkerSize',7.2) ;
ylabel('$z$','Interpreter','latex'); xlabel('$\psi(z)$','Interpreter','latex')
grid on ;
figure(5);
plot(tht,x,'mo--','LineWidth',1.2,'MarkerFaceColor','yellow','MarkerSize',7.2) ;
ylabel('$z$','Interpreter','latex'); xlabel('$\theta(z)$','Interpreter','latex')
grid on ;
figure(6)
plot(q,x2,'.g'); 
ylabel('$z$','Interpreter','latex'); xlabel('$q(z)$','Interpreter','latex')
%% Results
fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;
% if itest==0
%     save('solution_scenario1','x','x2','p','tht','q')
% else
%     save('solution_scenario2','x','x2','p','tht','q')
% end
toc % scenario 1: 1.3 sec; scenario 2: cca. 3 min, on a Laptop 16 GB, 2.11 GHz