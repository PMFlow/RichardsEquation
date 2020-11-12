function [c]=BGRW_1D_Alt_rand_benchmark(c0,cBC,dthtc,I,dx,dtt,Lc,q,D)
%% BGRW function (random) - transport step in coupled flow & trnasport 1D

%% Initial condition
N=10^24; n0=round(c0*N); nBC=round(cBC*N);
%% BGRW solution
u=q*dtt/Lc/dx;
ru=2*D*dtt/Lc/dx^2*ones(1,I);
n=n0; nn=zeros(1,I);
restr=0; restjump=0; restf=0;
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
%% Boundary conditions - concentration
%%%% BC_Bottom/Top
nn(1)=nBC(1);
nn(I)=nBC(I);
%% Source term - concentration
restf=dthtc(2:I-1)*N + 1e-3*nn(2:I-1)./(1+nn(2:I-1))*N*dtt/Lc+restf; nf=floor(restf);
nn(2:I-1)=nn(2:I-1)+nf;
n=nn;
c=n/N;
