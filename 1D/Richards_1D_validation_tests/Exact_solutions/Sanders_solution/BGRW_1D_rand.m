function [tht]=BGRW_1D_rand(tht0,I,dx,dt,V,D,K,Q)
% BGRW function (random) - transport step in GRW scheme for Richards 1D/Theta-formulation
%% Initial condition
N=10^24; n0=round(tht0*N);
%% BGRW solution
u=V*dt/dx;
ru=2*D*dt/dx^2;
if max(abs(u))>=ru
   max(u) 
   max(ru) 
   return 
end
n=n0; nn=zeros(1,I);
restr=0; restjump=0; 

for x=1:I
    if n(x) > 0
        r=ru(x);
        restr=n(x)*(1-r)+restr; nsta=floor(restr);
        restr=restr-nsta; njump=n(x)-nsta;
        nn(x)=nn(x)+nsta;
        if(njump)>0
            restjump=njump*0.5*(1-u(x)/r)+restjump;
            nj(1)=floor(restjump); restjump=njump-nj(1);
            nj(2)=floor(restjump); restjump=restjump-nj(2);
            if x==1
                nn(2)=nn(2)+nj(2);
            elseif x==I
                nn(I-1)=nn(I-1)+nj(1);
            else
                for i=1:2
                    xd=x+(2*i-3);
                    nn(xd)=nn(xd)+nj(i);
                end
            end
        end
    end
end
%% boundary conditions
%%%% BC_Top Boundary: Q=constant
nn(I)=nn(I-1)+(Q-K(I))/D(I)*dx*N; 
n=nn;
tht=n/N;