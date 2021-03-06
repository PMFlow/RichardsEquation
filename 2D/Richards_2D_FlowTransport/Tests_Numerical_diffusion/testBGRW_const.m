function [c]=testBGRW_const(I,J,xx,yy,dx,dy,Vx,Vy,D1,D2,U_MEAN)
%% 2-dim BGRW algorithm for advection-diffusion processes

%% Domain / discretization
X=xx(2:I-1); Y=yy(2:J-1);
Lx=I-2; Ly=J-2; T=3; S=1;
if I==201
    T=2.99;
end
N=10^24; % number of particles
%% Time-step size & number of time steps
Dfactor=2;
dt=Dfactor*(2*D1/dx^2+2*D2/dy^2); dt=1./dt;
TR=floor(T/dt); dTR=round(TR/3);
tsel=dTR:dTR:TR;
%% Initial condition
y0=(Ly+1)*0.7*dy; x0=(Lx+1)/2*dx; ti=1;
gss=Gauss_IC(ti,dx,dy,x0,y0,Lx,Ly,U_MEAN,D1,D2);
c0=gss/sum(sum(gss)); % IC given by a thin Gaussian distribution
n0=round(c0*N);
figure(2); hold all;
mesh(X,Y,c0')
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
zlabel('$c(x,z,t=0)$','Interpreter','latex'); grid on; view(115,15);
%% Initial moments
xr=(X-x0); yr=(Y-y0);
Mx0=sum(xr*n0)/N; Varx0=sum(xr.^2*n0)/N-Mx0*Mx0;
My0=sum(yr*n0')/N; Vary0=sum(yr.^2*n0')/N-My0*My0;
%% Initialization spatial moments of the solution
Dx=zeros(1,TR); Mx=zeros(1,TR); Varx=zeros(1,TR);
Dy=zeros(1,TR); My=zeros(1,TR); Vary=zeros(1,TR);
%% BGRW solution
v=Vy'*dt/dx;
u=Vx'*dt/dy;
ru=2*D1*dt/(dx*dx)*ones(Lx,Ly);
rv=2*D2*dt/(dy*dy)*ones(Lx,Ly);
n=n0; tr=dt; t=0; nn=zeros(Lx,Ly);
restr=0; restjump=0; restjumpx=0; restjumpy=0;
while tr<=T
    t=t+1;
    ntot=0;
    for s=1:S
        for x=1:(Lx-1)
            for y=1:(Ly-1)
                if n(x,y) > 0
                    rx=ru(x,y); ry=rv(x,y); r=rx+ry;
                    restr=n(x,y)*(1-r)+restr; nsta=floor(restr);
                    restr=restr-nsta; njump=n(x,y)-nsta;
                    nn(x,y)=nn(x,y)+nsta;
                    restjump=njump*rx/r+restjump;
                    njumpx=floor(restjump); restjump=restjump-njumpx;
                    njumpy=njump-njumpx;
                    if(njumpx)>0
                        restjumpx=njumpx*0.5*(1-u(x,y)/rx)+restjumpx;
                        nj(1)=floor(restjumpx); restjumpx=restjumpx-nj(1);
                        nj(2)=njumpx-nj(1);
                        for i=1:2
                            if nj(i)>0
                                xd=x+(2*i-3);
                                nn(xd,y)=nn(xd,y)+nj(i);
                            end
                        end
                    end
                    if(njumpy)>0
                        restjumpy=njumpy*0.5*(1-v(x,y)/ry)+restjumpy;
                        nj(1)=floor(restjumpy); restjumpy=restjumpy-nj(1);
                        nj(2)=njumpy-nj(1);
                        for i=1:2
                            if nj(i)>0
                                yd=y+(2*i-3);
                                nn(x,yd)=nn(x,yd)+nj(i);
                            end
                        end
                    end
                end
            end
        end
        n=nn; nn=zeros(Lx,Ly); ntot=sum(sum(n));
        Mx(t)=sum(xr*n); Varx(t)=sum(xr.^2*n)+sum(sum(u.^2.*n*dx^2*t)); %for Biased_GRW only !
        My(t)=sum(yr*n'); Vary(t)=sum(yr.^2*n')+sum(sum(v.^2.*n*dy^2*t));
        if s==1 && sum(ismember(t,tsel))
            [lin,col]=find(n>0);
            [nl,nc]=size(n);
            linind=sub2ind([nl,nc],lin,col);
            figure(3); plot3(lin*dx,col*dy,n(linind)/N);
            if t==tsel(1)
                hold;
            end
            grid on; view(115,15);
            xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
            zlabel('$c(x,z,t)$','Interpreter','latex');
        end
        tr=tr+dt;
    end
end
c=n/N;
tt=T+ti;
gss=Gauss_IC(tt,dx,dy,x0,y0,Lx,Ly,U_MEAN,D1,D2);
c_gss=gss/sum(sum(gss));
eps_c=norm(c-c_gss)/norm(c_gss);
figure(4)
mesh(X,Y,c'-c_gss'); view(115,15);
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
zlabel('$c(x,z,T)-c_{_{Gauss}}(x,z,T)$','Interpreter','latex');
for t=1:TR
    Mx(t)=Mx(t)/ntot; My(t)=My(t)/ntot;
    Varx(t)=Varx(t)/ntot-Mx(t)*Mx(t)-Varx0; Vary(t)=Vary(t)/ntot-My(t)*My(t)-Vary0;
    Dx(t)=Varx(t)/(2*t*dt); Dy(t)=Vary(t)/(2*t*dt);
    Mx(t)=(Mx(t)-Mx0)/(t*dt); My(t)=(My(t)-My0)/(t*dt);
end
t=1:TR;
figure(5);
subplot(1,2,1)
plot(t*dt,Mx,t*dt,My);
legend('$V_x(t)$','$V_z(t)$','Interpreter','latex','Location','best');
xlabel('$t$','Interpreter','latex');
subplot(1,2,2)
plot(t*dt,Dx,t*dt,Dy);
legend('$D_x(t)$','$D_z(t)$','Interpreter','latex','Location','best');
xlabel('$t$','Interpreter','latex');
figure(6); hold all;
subplot(1,2,1)
semilogy(t*dt,abs(Mx-0)/abs(U_MEAN),'.',t*dt,abs(My-U_MEAN)/abs(U_MEAN),'.');
legend('$|V_{x}(t) - V_{0x}| /|V_{0z}|$','$|V_{z}(t) - V_{0z}| /|V_{0z}|$','Interpreter','latex','Location','best');
xlabel('$t$','Interpreter','latex');
subplot(1,2,2)
semilogy(t*dt,abs(Dx-D1)/D1,'.',t*dt,abs(Dy-D2)/D2,'.');
legend('$|D_{x}(t) - D_{0x}| / D_{0x}$','$|D_{z}(t) - D_{0z}| / D_{0z}$','Interpreter','latex','Location','best');
xlabel('$t$','Interpreter','latex');

fprintf('eps_c = %0.2e \n',eps_c);
fprintf('eps_D1 = %0.2e \n',sqrt(dt)*norm(Dx-D1)/D1/T);
fprintf('eps_D2 = %0.2e \n',sqrt(dt)*norm(Dy-D2)/D2/T);
fprintf('The number of transport time steps is : %d \n',TR) ;
fprintf('The transport time step is : %0.2e \n',dt) ;

