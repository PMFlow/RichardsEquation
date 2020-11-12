%%   IC for Richards-1D as steady state solutions for time independent BC
close all
tic
itest=0; % 0=scenario1; 1=scenario2; 
%%   Grid Initialization
I = 101;
a = 0;
b = 2;
dx = (b-a)/(I-1)
x = [a:dx:b];
x2 = (x(1:I-1)+x(2:I))/2;
T=1;
S=1e7; 
if itest == 0
    pass=1e3; % for scenario1
else
    pass=1e5; % for scenario2
end
Tolerance=1e-9; 
%%   Parameters
Ksat = 2.77*10^-6;
theta_res=0.06;
theta_sat=0.36;
alpha=10;
q0=2.77*10^-7;
q1=2.50*10^-6;
maxr=0.8; 
%%  Coefficients and Forcing Term
vectKsat=Ksat*ones(1,I);
if itest == 1
    vectKsat((I-1)/2+1:I)=500*Ksat;
end
theta = @(p) theta_exp(theta_res,theta_sat,alpha,p);
K = @(theta) vectKsat.*(theta-theta_res)/(theta_sat-theta_res); 
%% solution
cmax=0.5;
L=1;
dN=1e24;
n0 = zeros(1,I); n0(1)=0.5*dN; n=n0; 
sumninit=sum(n0); 
pinit=n/dN; p=pinit;
figure(1); plot(pinit,x,'*');
ylabel('$z$','Interpreter','latex'); xlabel('$\psi(z,t=0)$','Interpreter','latex');
ks=1;
sumn=zeros(1,floor(S/pass)); tvect=zeros(1,floor(S/pass));
flux1=zeros(1,floor(S/pass)); fluxI=zeros(1,floor(S/pass));
tevol=zeros(1,floor(S/pass)); eps=zeros(1,floor(S/pass)); 
restr=0; restsar1=0; restsarI=0; rest1=0; rest2=0; restf=0;

for t=1:T
    for s=1:S
        pa=p;   
        tht=theta(p);
        D=K(tht); D=D(1:I-1);
        dt=dx^2*maxr/max(D)/2; % r<=1/2 such that 1-2*r>=0 
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
        nn(1)=n0(1);
        %% Neuman BC
        qR=q0; % water flux discretized as: -qR=-D(I-1)*((p(I)-p(I-1))/dx+1)
        nn(I)=nn(I-1)+(qR/D(I-1)-1)*dx*dN;
        %% Source term
        f=diff(r)*dx; % for steady-state problems the increments of theta vanish
        restf=dN*f+restf; nf=floor(restf); restf=restf-nf;
        nn(2:I-1)=nn(2:I-1)+nf;
        n=nn; p=n/dN; 
        na=n;
        %%
        if mod(s,pass)==0
            fprintf('s= %d\n',s);
            tvect(ks)=s;
            %%            
            sumn(ks)=sum(n);
            flux1(ks)=nsar(1)/2-nsar(2)*rapr(2);
            fluxI(ks)=nsar(I-1)*(1-rapr(I-1))-nsar(I)/2;
            %% Convergence criterion
            eps(ks)=norm(p-pa)/norm(p); 
            if eps(ks) <= Tolerance
                break
            end                       
            ks=ks+1;
        end        
    end
end
q=-D.*((p(2:I)-p(1:I-1))/dx+1);
%% plots
figure(2); plot(tvect,sumn,'o'); hold; plot(0,sumninit,'pr')
xlabel('$s$','Interpreter','latex'); 
title('Total number of particles')
figure(3); plot(tvect,flux1,'o'); hold all; plot(tvect,fluxI,'*')
xlabel('$s$','Interpreter','latex'); 
legend('In-flux (x=-1)','Out-flux (x=1)','Location','best')
title('In- and out-particle flux')
figure(4); plot(tvect,flux1-fluxI,'o')
xlabel('$s$','Interpreter','latex'); 
title('Total particle flux')
figure(5); semilogy(tvect,eps);
xlabel('$s$','Interpreter','latex'); ylabel('$\|\psi^s-\psi^{s-1}\| \;/\; \|\psi^s\|$','Interpreter','latex'); 
title('Convergence of the L-scheme');
sumnfin=round(sum(n))
figure(6);
plot(p,x,'mo--','LineWidth',1.2,'MarkerFaceColor','yellow','MarkerSize',7.2) ;
title('Initial pressure head');
ylabel('$z$','Interpreter','latex'); xlabel('$\psi(z)$','Interpreter','latex')
grid on ;
figure(7);
plot(tht,x,'mo--','LineWidth',1.2,'MarkerFaceColor','yellow','MarkerSize',7.2) ;
title('Initial water content');
ylabel('$z$','Interpreter','latex'); xlabel('$\theta(z)$','Interpreter','latex')
grid on ;
figure(8)
plot(q,x2,'.g'); 
title('Initial water flux');
ylabel('$z$','Interpreter','latex'); xlabel('$q(z)$','Interpreter','latex')
%% Results
fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;
if itest == 0
    save('IC_GRW_scenario1','x','x2','p','tht','q')
else
    save('IC_GRW_scenario2','x','x2','p','tht','q')
end
toc