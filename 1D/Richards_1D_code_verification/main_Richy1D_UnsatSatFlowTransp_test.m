%%   code validation: coupled flow & transport 1D: degenerate problem
%%   ================================================================

%%   Clearing Screen and Variables
clear all; close all
tic
Itest = 1;
L2_p = zeros(4,1); eoc_p=zeros(3,1);
L2_c = zeros(4,1); eoc_c=zeros(3,1);
for it = [10 20 40 80]
    %% Grid Initialization
    I = it+1;
    a=0; b=1;
    dx = (b-a)/(I-1);
    x = [a:dx:b];
    x2 = (x(1:I-1)+x(2:I))/2;
    %% Parameters
    D1=1.0;
    T = 1;
    past=T/3;
    S= 10000;
    maxr=1; 
    Tolerance = 1e-6;
    Lp=100; Lc=Lp; 
    Dfactor=2; 
    dtc=Dfactor*(2*D1/Lc/dx^2); dtc=1./dtc;
    %% Parameter functions and sources
    solEp = @(t,x)  -t.*x.*(1-x)+x/4;
    solEc = @(t,x)   t.*x.*(1-x)+1; 
    sol_q = @(t,x) - t*x - t*(x - 1) - 5/4;
    K = @(p) ones(1,length(p));
    %% Initial Conditions
    q=zeros(1,I);
    p0=solEp(0,x); 
    c0=solEc(0,x); 
    p=p0; c=c0; p0E=p0; c0E=c0;
    %% Solution
    tgraf=0; t=0; kt=1; test_kt=0;
    eps=zeros(1,S); epsc=zeros(1,S); convf=zeros(3,S); convt=zeros(3,S); nn=zeros(I);
    tht=theta(p,c); thtc=tht.*c;
    tht0=tht; thtc0=thtc;
    pa=p; ca=c;
    iS=1:S;   
    dN=10^24;
    n0 = floor(dN*p); n0E=round(p0E*dN);
    sumninit=sum(n0);
    n=n0;
    pinit=p;
    restr=0; restsar1=0; restsarI=0; rest1=0; rest2=0; restf=0;   
    fileID = fopen('test_dt','w');
    while t<=T
        DK=K(p); 
        D=(DK(1:I-1)+DK(2:I))/2;
        Dq=DK(2:I-1);
        [fs,fg]=Fsg(t,x,p); 
        dtp=Lp*dx^2*maxr/max(D)/2; 
        dt=min(dtc,dtp);
        t=t+dt;
        fprintf(fileID,'%2.2e  %2.2e\n',t, dt ) ;
        eps=zeros(1,S);
        for s=1:S
            DK=K(p); 
            D=(DK(1:I-1)+DK(2:I))/2;
            Dq=DK(2:I-1);            
            dtp=Lp*dx^2*maxr/max(D)/2; 
            dt=min(dtc,dtp);
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
            %% boundary conditions
            %%%% BC_Left/Right
            nn(1)=n0E(1);
            nn(I)=n0E(I);
            %% Source term
            dtht=(tht0-tht)/Lp;
            f=diff(r)*dx+dtht(2:I-1)+fs(2:I-1)*dt/Lp;
            restf=dN*f+restf; nf=floor(restf); restf=restf-nf;
            nn(2:I-1)=nn(2:I-1)+nf;
            n=nn; p=n/dN;
            thtc=theta(p,c).*c;  p0=p; %% Alt-scheme
            %% Transport step
            dthtc=(thtc0-thtc)/Lc;
            q(2:I-1)=-Dq.*((p0(3:I)-p0(1:I-2))/(2*dx)+1); % q - interior of \Omega
            q(1)=q(2); q(I)=q(I-1);                       % q - extended from interior of \Omega
%             q(1)=-DK(1).*((p0(2)-p0(1))/dx+1);            % q(1)- approximated by finite differences
%             q(I)=-DK(I).*((p0(I)-p0(I-1))/dx+1);          % q(I)- approximated by finite differences
%             q=sol_q(t,x);
            mean_q=mean(q)*ones(1,I);
            [c]=BGRW_1D_Alt_rand(c0,c0E,fg,dthtc,I,dx,dt,Lc,q,D1);
            c0=c; tht=theta(p,c);
            %% Convergence criterion
            tol_eps=norm(p-pa);
            tol_epsc=norm(c-ca);
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
            rndt=kt*past;
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
        tE=t;
    end
    if  kt==4
        save('convf');
        save('convt');
    end
    %% Results
    kt_plot=kt-1;
    plot_conv_Alt(kt_plot,S,strvect)
    fprintf('The space step is : %0.2e \n',dx) ;
    fprintf('The time step is : %0.2e \n',dt) ;
    fprintf('The total time is : %0.2e \n',tE) ;    
    cE_p=solEp(tE,x);
    L2_p(Itest) = ( dx )^(1/2) *norm(p-cE_p); 
    cE_c=solEc(tE,x);
    L2_c(Itest) = ( dx )^(1/2) *norm(c-cE_c); 
    if Itest >1
        eoc_p(Itest-1)=log10(L2_p(Itest-1)/L2_p(Itest))/log10(2);
        eoc_c(Itest-1)=log10(L2_c(Itest-1)/L2_c(Itest))/log10(2);
    end
    Itest = Itest+1;
end
fprintf('L2_p  : %0.2e \n',L2_p)
fprintf('EOC_p : %0.2e \n',eoc_p)
fprintf('L2_c  : %0.2e \n',L2_c)
fprintf('EOC_c : %0.2e \n',eoc_c)
max_grid_Peclet_x=dx*max(max(abs(q)))/D1
mean_grid_Peclet_x=dx*mean(mean(abs(q)))/D1
p1=solEp(1,x); c1=solEc(1,x);
figure;
plot(x,p1,'-*',x,p,'-+');
xlabel('$z$','Interpreter','latex');
ylabel('$\psi(z,t)$','Interpreter','latex'); 
legend('analytical','numerical','Location','best'); legend('boxoff')
figure;
plot(x,c1,'-*',x,c,'-+');
xlabel('$z$','Interpreter','latex');
ylabel('$c(z,t)$','Interpreter','latex'); 
legend('analytical','numerical','Location','best'); legend('boxoff')

toc
%% extended-q;
%% 10<=I<=80, Lp=Lc=100
% main_Richy1D_UnsatSatFlowTransp_test0
% The space step is : 1.00e-01 
% The time step is : 2.50e-01 
% The total time is : 1.00e+00 
% The space step is : 5.00e-02 
% The time step is : 6.25e-02 
% The total time is : 1.00e+00 
% The space step is : 2.50e-02 
% The time step is : 1.56e-02 
% The total time is : 1.00e+00 
% The space step is : 1.25e-02 
% The time step is : 3.91e-03 
% The total time is : 1.00e+00 
% L2_p  : 4.59e-02 
% L2_p  : 1.14e-02 
% L2_p  : 2.94e-03 
% L2_p  : 9.99e-04 
                    % EOC_p : 2.01e+00 
                    % EOC_p : 1.96e+00 
                    % EOC_p : 1.56e+00 
% L2_c  : 5.38e-02 
% L2_c  : 1.43e-02 
% L2_c  : 3.84e-03 
% L2_c  : 1.44e-03 
                    % EOC_c : 1.91e+00 
                    % EOC_c : 1.90e+00 
                    % EOC_c : 1.42e+00 
% max_grid_Peclet_x = 0.0277
% mean_grid_Peclet_x = 0.0156
% Elapsed time is 23.865306 seconds.
%% approximated-q
%% 10<=I<=80, Lp=Lc=100
% main_Richy1D_UnsatSatFlowTransp_test0
% The space step is : 1.00e-01 
% The time step is : 2.50e-01 
% The total time is : 1.00e+00 
% The space step is : 5.00e-02 
% The time step is : 6.25e-02 
% The total time is : 1.00e+00 
% The space step is : 2.50e-02 
% The time step is : 1.56e-02 
% The total time is : 1.00e+00 
% The space step is : 1.25e-02 
% The time step is : 3.91e-03 
% The total time is : 1.00e+00 
% L2_p  : 4.59e-02 
% L2_p  : 1.14e-02 
% L2_p  : 2.94e-03 
% L2_p  : 9.99e-04 
                    % EOC_p : 2.01e+00 
                    % EOC_p : 1.96e+00 
                    % EOC_p : 1.56e+00 
% L2_c  : 5.02e-02 
% L2_c  : 1.31e-02 
% L2_c  : 3.51e-03 
% L2_c  : 1.36e-03 
                    % EOC_c : 1.94e+00 
                    % EOC_c : 1.90e+00 
                    % EOC_c : 1.37e+00 
% max_grid_Peclet_x = 0.0279
% mean_grid_Peclet_x = 0.0156
% Elapsed time is 23.751488 seconds.
%% exact-q;
%% 10<=I<=80, Lp=Lc=100
% main_Richy1D_UnsatSatFlowTransp_test0
% The space step is : 1.00e-01 
% The time step is : 2.50e-01 
% The total time is : 1.00e+00 
% The space step is : 5.00e-02 
% The time step is : 6.25e-02 
% The total time is : 1.00e+00 
% The space step is : 2.50e-02 
% The time step is : 1.56e-02 
% The total time is : 1.00e+00 
% The space step is : 1.25e-02 
% The time step is : 3.91e-03 
% The total time is : 1.00e+00 
% L2_p  : 4.59e-02 
% L2_p  : 1.14e-02 
% L2_p  : 2.95e-03 
% L2_p  : 9.97e-04 
                    % EOC_p : 2.00e+00 
                    % EOC_p : 1.95e+00 
                    % EOC_p : 1.57e+00 
% L2_c  : 7.64e-03 
% L2_c  : 1.95e-03 
% L2_c  : 4.38e-04 
% L2_c  : 8.99e-05 
                    % EOC_c : 1.97e+00 
                    % EOC_c : 2.16e+00 
                    % EOC_c : 2.28e+00 
% max_grid_Peclet_x = 0.0281
% mean_grid_Peclet_x = 0.0156
% Elapsed time is 23.163856 seconds.

